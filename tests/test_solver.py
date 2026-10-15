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