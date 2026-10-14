"""The reference answers are checked against the equations they claim to solve,
so an error in a derivation shows up here, not in the solver tests."""
import numpy as np
import pytest

from canalweir.exact import (SteadyBenchmark, backwater_length, backwater_profile,
                             drawdown_head, friction_slope, level_pool_length, ritter)
from canalweir.hydraulics import G, WEIR_COEFF, normal_depth, weir_head


def test_ritter_satisfies_the_shallow_water_equations():
    h_left, x_dam, t, eps = 1.0, 0.0, 10.0, 1e-4
    c0 = np.sqrt(G * h_left)
    x = np.linspace(-0.9 * c0 * t, 1.9 * c0 * t, 50)      # inside the rarefaction

    def state(xx, tt):
        h, u = ritter(xx, tt, h_left, x_dam)
        return h, h * u

    h_t = (state(x, t + eps)[0] - state(x, t - eps)[0]) / (2 * eps)
    q_t = (state(x, t + eps)[1] - state(x, t - eps)[1]) / (2 * eps)
    hp, qp = state(x + eps, t)
    hm, qm = state(x - eps, t)
    mass = h_t + (qp - qm) / (2 * eps)
    flux_p = qp**2 / hp + 0.5 * G * hp**2
    flux_m = qm**2 / hm + 0.5 * G * hm**2
    momentum = q_t + (flux_p - flux_m) / (2 * eps)
    assert np.max(np.abs(mass)) < 1e-7
    assert np.max(np.abs(momentum)) < 1e-6


def test_ritter_conserves_volume():
    h_left, x_dam, t = 1.0, 500.0, 30.0
    x = np.linspace(0.0, 1000.0, 200001)
    h, _ = ritter(x, t, h_left, x_dam)
    volume = np.sum(0.5 * (h[1:] + h[:-1]) * np.diff(x))      # trapezoidal rule
    assert volume == pytest.approx(h_left * x_dam, rel=1e-6)


def test_ritter_front_positions():
    h_left, t = 2.0, 5.0
    c0 = np.sqrt(G * h_left)
    h, u = ritter(np.array([-c0 * t - 1e-9, 2 * c0 * t + 1e-9]), t, h_left, 0.0)
    assert h[0] == h_left and u[0] == 0.0          # undisturbed pool
    assert h[1] == 0.0 and u[1] == 0.0             # still dry


def test_steady_benchmark_satisfies_steady_momentum():
    b = SteadyBenchmark()
    x = np.linspace(1.0, 99.0, 40)
    eps = 1e-3

    def momentum_flux(xx):
        h = b.depth(xx)
        return b.q**2 / h + 0.5 * G * h**2

    dflux = (momentum_flux(x + eps) - momentum_flux(x - eps)) / (2 * eps)
    dz = (b.bed(x + eps) - b.bed(x - eps)) / (2 * eps)
    h = b.depth(x)
    residual = dflux + G * h * dz + G * h * friction_slope(h, b.q, b.n)
    assert np.max(np.abs(residual)) < 1e-6


def test_backwater_profile_returns_to_normal_depth():
    q, slope, n, width = 0.2, 2e-4, 0.03, 3.0
    x, h = backwater_profile(q, slope, n, 1.0, 40_000.0, width)
    assert h[-1] == pytest.approx(1.0 + weir_head(q), rel=1e-12)
    assert h[0] == pytest.approx(normal_depth(q, slope, n, width), abs=1e-6)
    assert np.all(np.diff(h) > -1e-12)              # M1 curve: never shallower downstream


def test_backwater_length_of_a_straight_line():
    x = np.linspace(0.0, 100.0, 101)
    h = 1.0 + 0.01 * x                              # 0.05 above 1.0 at x = 5
    assert backwater_length(x, h, 1.0, 100.0, 0.05) == pytest.approx(95.0, rel=1e-12)


def test_backwater_length_is_zero_when_nothing_is_raised():
    x = np.linspace(0.0, 10.0, 11)
    assert backwater_length(x, np.ones(11), 1.0, 10.0) == 0.0


def test_level_pool_length():
    q, slope, n = 0.2, 2e-4, 0.03
    expected = (1.0 + weir_head(q) - normal_depth(q, slope, n)) / slope
    assert level_pool_length(q, slope, n, 1.0) == pytest.approx(expected, rel=1e-14)


def test_drawdown_head_satisfies_the_mass_balance():
    head0, length, eps = 0.2, 100.0, 1e-3
    t = np.linspace(0.0, 500.0, 11)
    later = drawdown_head(t + eps, head0, length)
    earlier = drawdown_head(t - eps, head0, length)
    dHdt = (later - earlier) / (2 * eps)
    residual = length * dHdt + WEIR_COEFF * drawdown_head(t, head0, length) ** 1.5
    assert drawdown_head(0.0, head0, length) == head0
    assert np.max(np.abs(residual)) < 1e-8