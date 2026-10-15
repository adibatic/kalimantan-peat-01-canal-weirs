"""Verification: observed order of accuracy on a steady flow with friction.

The SteadyBenchmark has an exact depth h(x); the bed is built so that h(x)
satisfies the steady equations. Halving dx should divide the error by 2
(order 1) or 4 (order 2).
"""
from __future__ import annotations

import csv

import numpy as np
from common import BLUE, MUTED, ORANGE, TABLES, plt, save

from canalweir.exact import SteadyBenchmark
from canalweir.solver import Model, Reach

CELLS = (25, 50, 100, 200, 400)


def run(bench, cells, order):
    x = (np.arange(cells) + 0.5) * bench.length / cells
    reach = Reach(bench.length, cells, bench.bed(x), manning_n=bench.n,
                  west=("inflow", bench.q),
                  east=("depth", float(bench.depth(bench.length))))
    model = Model([reach], [bench.depth(x)], [np.full(cells, bench.q)], order=order)
    model.run_to_steady(tol=1e-12)
    return np.mean(np.abs(model.h[0] - bench.depth(x)))


def main():
    bench = SteadyBenchmark()
    rows, errors = [], {1: [], 2: []}
    for order in (1, 2):
        previous = None
        for cells in CELLS:
            error = run(bench, cells, order)
            p = np.log2(previous / error) if previous else float("nan")
            rows.append({"order": order, "cells": cells, "dx_m": bench.length / cells,
                         "L1_error_m": error, "observed_order": p})
            errors[order].append(error)
            print(f"order {order}  cells {cells:4d}  L1 {error:.3e} m  p {p:.3f}")
            previous = error
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "steady_convergence.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    dx = bench.length / np.array(CELLS)
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for order, color in ((1, ORANGE), (2, BLUE)):
        e = np.array(errors[order])
        ax.loglog(dx, e, "o-", color=color, ms=6, label=f"order {order} scheme")
        ref = 0.4 * e[0] * (dx / dx[0]) ** order       # offset below the data
        ax.loglog(dx, ref, ls=":", color=MUTED, lw=1.2)
        ax.annotate(f"slope {order}", (dx[-1], ref[-1]), textcoords="offset points",
                    xytext=(4, -12), color=MUTED, fontsize=9)
    ax.set_ylim(bottom=1.5e-7)                       # room for the slope labels
    ax.set_xticks(dx, [f"{d:g}" for d in dx])
    ax.minorticks_off()
    ax.set(xlabel="cell size, dx (m)", ylabel="mean depth error (m)",
           title="Order of accuracy, steady flow with friction")
    ax.legend(loc="upper left")
    save(fig, "steady_convergence.png")


if __name__ == "__main__":
    main()