"""1-D finite-volume shallow-water solver for canal reaches separated by weirs.

Unknowns per cell: depth h (m) and unit discharge q = h u (m^2/s).

    h_t + q_x = 0
    q_t + (q^2/h + G h^2/2)_x = -G h z_x - G h S_f

Scheme (sources in README, Sources; all [VERIFY]):
- HLL approximate Riemann solver (Harten, Lax and van Leer 1983).
- Hydrostatic reconstruction of the bed slope term (Audusse et al. 2004),
  which keeps a lake at rest exactly at rest.
- order=1: piecewise-constant states, forward Euler.
  order=2: minmod-limited linear states (MUSCL), Heun's method.
- Manning friction treated semi-implicitly inside each stage.

Reaches are joined in series: the weir outflow of one reach is the inflow
of the next, exchanged inside every stage, so volume is conserved exactly.

Boundary conditions, given as tuples:
    ("wall",)                 no flow
    ("inflow", q)             prescribed inflow, q a number or a function of t
    ("depth", h)              prescribed depth at the face (subcritical outflow)
    ("weir", crest, cd)       free overflow, weir_unit_discharge(level, crest, cd)
    ("upstream",)             fed by the weir of the reach upstream
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .hydraulics import G, hydraulic_radius, weir_unit_discharge

DRY = 1e-8  # m; below this depth a cell is dry and its velocity is zero


def minmod(a, b):
    return np.where(a * b > 0.0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def limited_slopes(v):
    """Minmod slope of v in every cell; the two end cells use one-sided slopes."""
    d = np.diff(v)
    s = np.empty_like(v)
    s[1:-1] = minmod(d[:-1], d[1:])
    s[0], s[-1] = d[0], d[-1]
    return s


def velocity(h, q):
    return np.where(h > DRY, q / np.maximum(h, DRY), 0.0)


def hll(h_l, q_l, h_r, q_r):
    """HLL flux of the shallow-water equations between a left and a right state."""
    u_l, u_r = velocity(h_l, q_l), velocity(h_r, q_r)
    c_l, c_r = np.sqrt(G * h_l), np.sqrt(G * h_r)
    s_l = np.minimum(u_l - c_l, u_r - c_r)
    s_r = np.maximum(u_l + c_l, u_r + c_r)
    fh_l, fq_l = q_l, q_l * u_l + 0.5 * G * h_l**2
    fh_r, fq_r = q_r, q_r * u_r + 0.5 * G * h_r**2
    span = np.where(s_r - s_l > 0.0, s_r - s_l, 1.0)
    mid_h = (s_r * fh_l - s_l * fh_r + s_l * s_r * (h_r - h_l)) / span
    mid_q = (s_r * fq_l - s_l * fq_r + s_l * s_r * (q_r - q_l)) / span
    fh = np.where(s_l >= 0.0, fh_l, np.where(s_r <= 0.0, fh_r, mid_h))
    fq = np.where(s_l >= 0.0, fq_l, np.where(s_r <= 0.0, fq_r, mid_q))
    return fh, fq


@dataclass
class Reach:
    """A straight rectangular canal reach divided into equal cells."""

    length: float                 # m
    cells: int
    bed: np.ndarray               # bed elevation at the cell centres, m
    manning_n: float = 0.0        # s/m^(1/3); 0 means frictionless
    width: float | None = None    # m; None means a wide canal (R = h)
    west: tuple = ("wall",)
    east: tuple = ("wall",)

    def __post_init__(self):
        self.dx = self.length / self.cells
        self.x = (np.arange(self.cells) + 0.5) * self.dx
        self.bed = np.asarray(self.bed, dtype=float)
        if self.bed.shape != (self.cells,):
            raise ValueError("bed needs one value per cell")


class Model:
    """Reaches in series, advanced together with one time step."""

    def __init__(self, reaches, depths, discharges, order=2, cfl=0.45):
        if order not in (1, 2):
            raise ValueError("order must be 1 or 2")
        self.reaches = list(reaches)
        self.h = [np.array(d, dtype=float) for d in depths]
        self.q = [np.array(d, dtype=float) for d in discharges]
        self.order, self.cfl = order, cfl
        self.t = 0.0
        self.inflow_volume = 0.0    # m^3 (m^2 per metre of width if width is None)
        self.outflow_volume = 0.0
        self.outflows = [0.0] * len(self.reaches)   # east-end discharge, last step

    # ---------- fluxes ----------

    @staticmethod
    def _boundary(bc, h, u, level, t, q_feed):
        """Mass and momentum flux through a boundary face (positive along +x)."""
        kind = bc[0]
        if kind == "wall":
            return 0.0, 0.5 * G * h * h
        if kind in ("inflow", "upstream"):
            qb = q_feed if kind == "upstream" else (bc[1](t) if callable(bc[1]) else bc[1])
            return qb, (qb * qb / h if h > DRY else 0.0) + 0.5 * G * h * h
        if kind == "depth":
            hb, qb = bc[1], h * u
            return qb, qb * qb / hb + 0.5 * G * hb * hb
        if kind == "weir":
            qb = float(weir_unit_discharge(level, bc[1], bc[2] if len(bc) > 2 else 1.0))
            return qb, (qb * qb / h if h > DRY else 0.0) + 0.5 * G * h * h
        raise ValueError(f"unknown boundary condition {kind!r}")

    def _fluxes(self, r, h, q, t, q_feed):
        u = velocity(h, q)
        level = h + r.bed
        if self.order == 2:
            sh, su, sl = limited_slopes(h), limited_slopes(u), limited_slopes(level)
        else:
            sh = su = sl = np.zeros_like(h)
        # states at the west (w) and east (e) face of every cell
        h_w, h_e = np.maximum(h - 0.5 * sh, 0.0), np.maximum(h + 0.5 * sh, 0.0)
        u_w, u_e = u - 0.5 * su, u + 0.5 * su
        z_w, z_e = level - 0.5 * sl - h_w, level + 0.5 * sl - h_e

        # hydrostatic reconstruction at the interior faces
        z_star = np.maximum(z_e[:-1], z_w[1:])
        hl = np.maximum(h_e[:-1] + z_e[:-1] - z_star, 0.0)
        hr = np.maximum(h_w[1:] + z_w[1:] - z_star, 0.0)
        fh, fq = hll(hl, hl * u_e[:-1], hr, hr * u_w[1:])

        n = r.cells
        west_h, west_q = np.empty(n), np.empty(n)
        east_h, east_q = np.empty(n), np.empty(n)
        east_h[:-1], east_q[:-1] = fh, fq + 0.5 * G * (h_e[:-1] ** 2 - hl**2)
        west_h[1:], west_q[1:] = fh, fq + 0.5 * G * (h_w[1:] ** 2 - hr**2)
        west_h[0], west_q[0] = self._boundary(r.west, h_w[0], u_w[0],
                                              h_w[0] + z_w[0], t, q_feed)
        east_h[-1], east_q[-1] = self._boundary(r.east, h_e[-1], u_e[-1],
                                                h_e[-1] + z_e[-1], t, None)
        source = -G * 0.5 * (h_w + h_e) * (z_e - z_w)   # zero when order == 1
        return west_h, west_q, east_h, east_q, source

    # ---------- time stepping ----------

    def _stage(self, hs, qs, t, dt):
        """One forward-Euler stage for every reach, upstream to downstream."""
        new_h, new_q, ends = [], [], []
        q_feed = 0.0
        for r, h, q in zip(self.reaches, hs, qs):
            west_h, west_q, east_h, east_q, source = self._fluxes(r, h, q, t, q_feed)
            h1 = h - dt / r.dx * (east_h - west_h)
            q1 = q - dt / r.dx * (east_q - west_q - source)
            if r.manning_n > 0.0:
                hp = np.maximum(h1, DRY)
                radius = hydraulic_radius(hp, r.width)
                k = G * r.manning_n**2 * np.abs(q) / (hp * radius ** (4.0 / 3.0))
                q1 = q1 / (1.0 + dt * k)
            q1 = np.where(h1 > DRY, q1, 0.0)
            new_h.append(np.maximum(h1, 0.0))
            new_q.append(q1)
            ends.append((west_h[0], east_h[-1]))
            q_feed = east_h[-1]
        return new_h, new_q, ends

    def max_wave_speed(self):
        speeds = [np.max((np.abs(velocity(h, q)) + np.sqrt(G * h)) / r.dx)
                  for r, h, q in zip(self.reaches, self.h, self.q)]
        return max(speeds)

    def stable_dt(self):
        return self.cfl / self.max_wave_speed()

    def _width_factor(self, r):
        return 1.0 if r.width is None else r.width

    def step(self, dt):
        """Advance by dt; returns the largest change rate of h (m/s)."""
        h0, q0 = self.h, self.q
        h1, q1, ends1 = self._stage(h0, q0, self.t, dt)
        if self.order == 1:
            h_new, q_new, w = h1, q1, [(1.0, ends1)]
        else:
            h2, q2, ends2 = self._stage(h1, q1, self.t + dt, dt)
            h_new = [0.5 * (a + b) for a, b in zip(h0, h2)]
            q_new = [0.5 * (a + b) for a, b in zip(q0, q2)]
            w = [(0.5, ends1), (0.5, ends2)]
        widths = [self._width_factor(r) for r in self.reaches]
        self.outflows = [sum(weight * ends[k][1] for weight, ends in w) * widths[k]
                         for k in range(len(self.reaches))]
        inflow = sum(weight * ends[0][0] for weight, ends in w) * widths[0]
        self.inflow_volume += dt * inflow
        self.outflow_volume += dt * self.outflows[-1]
        rate = max(np.max(np.abs(a - b)) for a, b in zip(h_new, h0)) / dt
        self.h, self.q = h_new, q_new
        self.t += dt
        return rate

    def run(self, t_end, callback=None, every=None):
        """Advance to t_end. callback(model) runs every `every` seconds."""
        next_out = self.t if every is not None else None
        while self.t < t_end - 1e-12:
            if callback is not None and next_out is not None and self.t >= next_out - 1e-9:
                callback(self)
                next_out += every
            self.step(min(self.stable_dt(), t_end - self.t))
        if callback is not None:
            callback(self)

    def run_to_steady(self, tol=1e-9, t_max=1e7):
        """Advance until every |dh/dt| is below tol (m/s). Returns the time."""
        while self.t < t_max:
            if self.step(self.stable_dt()) < tol:
                return self.t
        raise RuntimeError(f"no steady state before t = {t_max} s")

    # ---------- diagnostics ----------

    def volume(self):
        """Stored water volume, m^3 (m^2 per metre of width for wide canals)."""
        return sum(np.sum(h) * r.dx * self._width_factor(r)
                   for r, h in zip(self.reaches, self.h))