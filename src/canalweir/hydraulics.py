"""Closed-form open-channel relations for a rectangular canal.

All depths in m, unit discharges q in m^2/s (discharge per metre of width),
slopes dimensionless. Every function is checked by tests/test_hydraulics.py.

Sources (see README, Sources):
- Manning's equation and critical depth: Chow (1959) [VERIFY].
- Weir law: derived here from the energy equation with critical flow on a
  broad crest; compare with Chow (1959) [VERIFY].
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

G = 9.81  # gravitational acceleration, m/s^2

# Ideal broad-crested weir, per metre of crest: q = WEIR_COEFF * H^(3/2).
# On the crest the flow is critical, so the depth there is 2/3 of the upstream
# head H, and q = sqrt(G * (2/3 H)^3) = (2/3)^(3/2) sqrt(G) H^(3/2).
WEIR_COEFF = (2.0 / 3.0) ** 1.5 * np.sqrt(G)


def hydraulic_radius(h, width=None):
    """Flow area over wetted perimeter. width=None means a wide canal (R = h)."""
    h = np.asarray(h, dtype=float)
    if width is None:
        return h
    return width * h / (width + 2.0 * h)


def manning_unit_discharge(h, slope, n, width=None):
    """Uniform-flow discharge per metre of width, q = h R^(2/3) S^(1/2) / n."""
    h = np.asarray(h, dtype=float)
    return h * hydraulic_radius(h, width) ** (2.0 / 3.0) * np.sqrt(slope) / n


def normal_depth(q, slope, n, width=None):
    """Depth of uniform flow carrying q on a bed of the given slope."""
    def residual(h):
        return manning_unit_discharge(h, slope, n, width) - q

    upper = 1.0
    while residual(upper) < 0.0:
        upper *= 2.0
    return brentq(residual, 1e-12, upper, xtol=1e-14, rtol=1e-14)


def critical_depth(q):
    """Depth at which the Froude number is 1 in a rectangular canal."""
    return (q * q / G) ** (1.0 / 3.0)


def froude(h, q):
    """Froude number u / sqrt(G h) with u = q / h."""
    h = np.asarray(h, dtype=float)
    return q / (h * np.sqrt(G * h))


def weir_unit_discharge(level, crest, cd=1.0):
    """Free overflow per metre of crest for a water level above the crest.

    level and crest are elevations (m); the head is level - crest, zero when
    the water is below the crest. cd scales the ideal law (cd = 1).
    The approach velocity is neglected and the weir is never drowned.
    """
    head = np.maximum(np.asarray(level, dtype=float) - crest, 0.0)
    return cd * WEIR_COEFF * head**1.5


def weir_head(q, cd=1.0):
    """Head over the crest that passes q; the inverse of weir_unit_discharge."""
    return (q / (cd * WEIR_COEFF)) ** (2.0 / 3.0)