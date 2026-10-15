"""Verification: a weir fails and its pool spills onto a dry canal bed.

Compares the solver with Ritter's exact solution at t = 30 s and measures the
observed order of convergence of the L1 depth error for both schemes.
"""
from __future__ import annotations

import csv

import numpy as np
from common import BLUE, INK, MUTED, ORANGE, TABLES, plt, save

from canalweir.exact import ritter
from canalweir.solver import Model, Reach

LENGTH, X_DAM, H_POOL, T_END = 1000.0, 500.0, 1.0, 30.0
CELLS = (100, 200, 400, 800, 1600)


def run(cells, order):
    reach = Reach(LENGTH, cells, np.zeros(cells))
    depth = np.where(reach.x < X_DAM, H_POOL, 0.0)
    model = Model([reach], [depth], [np.zeros(cells)], order=order)
    volume0 = model.volume()
    model.run(T_END)
    exact, _ = ritter(reach.x, T_END, H_POOL, X_DAM)
    error = np.sum(np.abs(model.h[0] - exact)) * reach.dx / LENGTH
    return reach, model, error, abs(model.volume() - volume0) / volume0


def main():
    rows = []
    for order in (1, 2):
        previous = None
        for cells in CELLS:
            _, _, error, drift = run(cells, order)
            p = np.log2(previous / error) if previous else float("nan")
            rows.append({"order": order, "cells": cells, "dx_m": LENGTH / cells,
                         "L1_error_m": error, "observed_order": p,
                         "volume_drift": drift})
            print(f"order {order}  cells {cells:5d}  L1 {error:.3e} m  p {p:.2f}  "
                  f"volume drift {drift:.1e}")
            previous = error
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "dam_break.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    x = np.linspace(380, 720, 2000)
    exact, _ = ritter(x, T_END, H_POOL, X_DAM)
    ax.plot(x, exact, color=INK, lw=1.2, label="Ritter, exact")
    for order, color in ((1, ORANGE), (2, BLUE)):
        reach, model, *_ = run(200, order)
        keep = (reach.x > 380) & (reach.x < 720)
        ax.plot(reach.x[keep], model.h[0][keep], "o", ms=3.2, color=color,
                label=f"solver, order {order}, 200 cells")
    ax.axvline(X_DAM, color=MUTED, lw=1, ls="--")
    ax.text(X_DAM + 4, 0.93, "failed weir", color=MUTED, fontsize=9)
    ax.set(xlabel="distance along canal, x (m)", ylabel="depth, h (m)",
           title=f"Dam break onto a dry bed, t = {T_END:.0f} s")
    ax.legend(loc="upper right")
    save(fig, "dam_break.png")


if __name__ == "__main__":
    main()