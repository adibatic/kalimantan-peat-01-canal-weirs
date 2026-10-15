"""The solver against exact answers. Tolerances are fixed in CLAUDE.md."""
import numpy as np
import pytest

from canalweir.exact import SteadyBenchmark, backwater_length, backwater_profile, ritter
from canalweir.hydraulics import normal_depth, weir_head
from canalweir.solver import Model, Reach


def bumpy_reach(cells=60, **kw):
    reach = Reach(100.0, cells, np.zeros(cells), **kw)
    reach.bed = 0.3 * np.exp(-(((reach.x - 50.0) / 10.0) ** 2))
    return reach


@pytest.mark.parametrize("order", [1, 2])
def test_lake_at_rest_stays_at_rest(order):
    reach = bumpy_reach(manning_n=0.03)
    model = Model([reach], [1.0 - reach.bed], [np.zeros(reach.cells)], order=order)
    model.run(300.0)
    assert np.max(np.abs(model.q[0])) < 1e-12
    assert np.max(np.abs(model.h[0] + reach.bed - 1.0)) < 1e-12


@pytest.mark.parametrize("order", [1, 2])
def test_closed_canal_conserves_volume(order):
    reach = bumpy_reach(manning_n=0.03)
    depth = 1.0 - reach.bed + 0.2 * np.exp(-(((reach.x - 20.0) / 5.0) ** 2))
    model = Model([reach], [depth], [np.zeros(reach.cells)], order=order)
    v0 = model.volume()
    model.run(120.0)
    assert abs(model.volume() - v0) / v0 < 1e-13


@pytest.mark.parametrize("order", [1, 2])
def test_volume_balance_with_inflow_and_weirs(order):
    reaches = []
    for k in range(2):
        r = Reach(500.0, 25, np.zeros(25), manning_n=0.03, width=3.0,
                  west=("inflow", lambda t: 0.1 + 0.05 * np.sin(t / 50.0)) if k == 0
                  else ("upstream",), east=("weir", 0.8, 1.0))
        reaches.append(r)
    model = Model(reaches, [np.full(25, 0.9)] * 2, [np.zeros(25)] * 2, order=order)
    v0 = model.volume()
    model.run(900.0)
    balance = model.volume() - v0 - model.inflow_volume + model.outflow_volume
    assert abs(balance) / v0 < 1e-13


@pytest.mark.parametrize("order, max_error", [(1, 6e-3), (2, 3e-3)])
def test_dam_break_matches_ritter(order, max_error):
    reach = Reach(1000.0, 200, np.zeros(200))
    model = Model([reach], [np.where(reach.x < 500.0, 1.0, 0.0)], [np.zeros(200)],
                  order=order)
    model.run(30.0)
    exact, _ = ritter(reach.x, 30.0, 1.0, 500.0)
    error = np.sum(np.abs(model.h[0] - exact)) * reach.dx / reach.length
    assert error < max_error
    assert np.all(model.h[0] >= 0.0)


def steady_error(cells, order):
    b = SteadyBenchmark()
    x = (np.arange(cells) + 0.5) * b.length / cells
    reach = Reach(b.length, cells, b.bed(x), manning_n=b.n, west=("inflow", b.q),
                  east=("depth", float(b.depth(b.length))))
    model = Model([reach], [b.depth(x)], [np.full(cells, b.q)], order=order)
    model.run_to_steady(tol=1e-12)
    return np.mean(np.abs(model.h[0] - b.depth(x)))


@pytest.mark.parametrize("order", [1, 2])
def test_observed_order_of_accuracy(order):
    e = [steady_error(cells, order) for cells in (25, 50, 100)]
    observed = np.log2(e[1] / e[2])
    assert abs(observed - order) < 0.1


def test_steady_weir_passes_the_inflow_at_the_rating_head():
    q, crest, cells = 0.2, 1.0, 40
    reach = Reach(400.0, cells, np.zeros(cells), manning_n=0.03, width=3.0,
                  west=("inflow", q), east=("weir", crest, 1.0))
    model = Model([reach], [np.full(cells, 1.3)], [np.zeros(cells)])
    model.run_to_steady(tol=1e-11)
    assert model.outflows[0] / 3.0 == pytest.approx(q, rel=1e-8)
    # flat bed, slow flow: the level at the weir is the crest plus the rating head
    assert model.h[0][-1] == pytest.approx(crest + weir_head(q), abs=2e-3)


def test_backwater_length_agrees_with_gvf():
    q, slope, n, width, crest, length, cells = 0.1, 2e-4, 0.03, 3.0, 0.5, 3000.0, 75
    reach = Reach(length, cells, np.zeros(cells), manning_n=n, width=width,
                  west=("inflow", q), east=("weir", crest, 1.0))
    reach.bed = slope * (length - reach.x)
    h_n = normal_depth(q, slope, n, width)
    model = Model([reach], [np.maximum(h_n, crest + weir_head(q) - reach.bed)],
                  [np.full(cells, q)])
    model.run_to_steady(tol=1e-10)
    x, h = backwater_profile(q, slope, n, crest, length, width)
    L_gvf = backwater_length(x, h, h_n, length)
    L_fv = backwater_length(reach.x, model.h[0], h_n, length)
    assert abs(L_fv - L_gvf) / L_gvf < 0.01