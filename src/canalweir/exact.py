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