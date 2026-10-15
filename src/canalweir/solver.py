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