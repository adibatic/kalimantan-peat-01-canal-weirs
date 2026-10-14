"""Exact and reference answers the solver is checked against.

Each one is derived from the 1-D shallow-water equations (see README, Method)
and each is checked against those equations by tests/test_exact.py, so none
depends on a value recalled from memory. Sources to compare with are given
per function and in the README, Sources.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp

from .hydraulics import G, WEIR_COEFF, hydraulic_radius, normal_depth, weir_head


def ritter(x, t, h_left, x_dam):
    """Dam break onto a dry, flat, frictionless bed (Ritter 1892) [VERIFY].

    Water of depth h_left at rest for x < x_dam, dry for x > x_dam, at t = 0.
    Inside the rarefaction the left Riemann invariant u + 2c = 2 c0 holds and
    the characteristic through the dam gives x - x_dam = (u - c) t.
    Returns depth h (m) and velocity u (m/s).
    """
    x = np.asarray(x, dtype=float)
    c0 = np.sqrt(G * h_left)
    xi = (x - x_dam) / t
    fan = (xi > -c0) & (xi < 2.0 * c0)
    h = np.where(xi <= -c0, h_left, 0.0)
    h = np.where(fan, (2.0 * c0 - xi) ** 2 / (9.0 * G), h)
    u = np.where(fan, 2.0 / 3.0 * (c0 + xi), 0.0)
    return h, u


def friction_slope(h, q, n, width=None):
    """Manning friction slope S_f = n^2 q|q| / (h^2 R^(4/3))."""
    h = np.asarray(h, dtype=float)
    return n * n * q * np.abs(q) / (h * h * hydraulic_radius(h, width) ** (4.0 / 3.0))


class SteadyBenchmark:
    """A steady flow with friction whose depth is known exactly.

    Choose the depth h(x) and the discharge q; the steady momentum equation
    then gives the bed slope that makes h(x) exact (the approach of
    MacDonald et al. 1997, collected in Delestre et al. 2013) [VERIFY]:

        dz/dx = (q^2 / (G h^3) - 1) dh/dx - S_f(h, q)

    Here h(x) = h0 + dh (x / L)^2 on 0 <= x <= L, deepening downstream as in
    the pool behind a weir.
    """

    def __init__(self, length=100.0, q=0.2, n=0.03, h0=0.6, dh=0.2, width=None):
        self.length, self.q, self.n = length, q, n
        self.h0, self.dh, self.width = h0, dh, width

    def depth(self, x):
        return self.h0 + self.dh * (np.asarray(x, dtype=float) / self.length) ** 2

    def depth_gradient(self, x):
        return 2.0 * self.dh * np.asarray(x, dtype=float) / self.length**2

    def bed_slope(self, x):
        h = self.depth(x)
        return ((self.q**2 / (G * h**3) - 1.0) * self.depth_gradient(x)
                - friction_slope(h, self.q, self.n, self.width))

    def bed(self, x):
        """Bed elevation z(x) with z(0) = 0, by adaptive quadrature."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        order = np.argsort(x)
        z = np.empty_like(x)
        start, total = 0.0, 0.0
        for k in order:
            total += quad(lambda s: float(self.bed_slope(s)), start, x[k],
                          epsabs=1e-14, epsrel=1e-13)[0]
            z[k], start = total, x[k]
        return z


def backwater_profile(q, slope, n, crest_height, length, width=None, cd=1.0,
                      points=2001):
    """Steady depth upstream of a weir, by integrating the gradually varied
    flow equation from the weir towards upstream (Chow 1959) [VERIFY]:

        dh/dx = (S0 - S_f) / (1 - q^2 / (G h^3))

    The weir sits at x = length on a bed of constant slope S0; the depth there
    is crest_height + weir_head(q). Returns x (m, increasing downstream) and h.
    """
    h_weir = crest_height + weir_head(q, cd)

    def rhs(x, h):
        return (slope - friction_slope(h, q, n, width)) / (1.0 - q * q / (G * h**3))

    x = np.linspace(length, 0.0, points)
    sol = solve_ivp(rhs, (length, 0.0), [h_weir], t_eval=x, method="DOP853",
                    rtol=1e-11, atol=1e-12)
    return sol.t[::-1], sol.y[0][::-1]


def backwater_length(x, h, h_normal, x_weir, threshold=0.05):
    """Distance upstream of the weir over which h exceeds h_normal + threshold.

    x increases downstream and h is assumed to rise towards the weir (an M1
    profile), so the raised part is one stretch ending at the weir.
    """
    raised = np.asarray(h) - h_normal > threshold
    if not raised.any():
        return 0.0
    if raised[0]:
        raise ValueError("water is raised over the whole domain; make it longer")
    first = np.argmax(raised)              # first raised point from upstream
    x0, x1 = x[first - 1], x[first]
    e0, e1 = h[first - 1] - h_normal - threshold, h[first] - h_normal - threshold
    return x_weir - (x0 - e0 * (x1 - x0) / (e1 - e0))


def level_pool_length(q, slope, n, crest_height, width=None, cd=1.0):
    """Quick estimate: a flat pool at the weir's upstream level, cut off where
    it meets the normal-depth surface, (crest + head - h_n) / S0."""
    h_n = normal_depth(q, slope, n, width)
    return (crest_height + weir_head(q, cd) - h_n) / slope