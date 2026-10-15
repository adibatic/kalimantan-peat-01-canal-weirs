"""Shared settings for the experiments: the illustrative canal and the figure style.

The canal values are chosen for illustration, not measured on a peat canal.
Project 03 replaces them with surveyed values; every result here scales with them.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from canalweir.hydraulics import normal_depth, weir_head  # noqa: E402
from canalweir.solver import Model, Reach  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "results" / "figures"
TABLES = ROOT / "results" / "tables"
ANIMATIONS = ROOT / "results" / "animations"
MEDIA = ROOT / "media"

# ---------- the illustrative canal (chosen, not measured) ----------
WIDTH = 3.0                  # m, rectangular canal
MANNING_N = 0.03             # s/m^(1/3)
SLOPE = 2e-4                 # bed slope, m per m (0.2 m per km)
DISCHARGES = (0.1, 0.3, 0.6, 1.2)   # m^3/s, dry season to wet season
CREST_HEIGHTS = (0.5, 1.0, 1.5)     # m above the bed at the weir
THRESHOLD = 0.05             # m; "raised" means at least this much above normal depth

# ---------- figure style ----------
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#8f8d87"
SURFACE, GRID = "#fcfcfb", "#e4e2dc"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
BLUE_RAMP = ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")   # light to dark, ordinal
WATER, PEAT = "#9ec5f4", "#6b5444"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "figure.dpi": 110, "savefig.dpi": 150,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "lines.linewidth": 2.0, "legend.frameon": False,
})


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path.relative_to(ROOT).as_posix())


def sloped_reach(length, cells, crest_height, inflow, slope=SLOPE, n=MANNING_N,
                 width=WIDTH):
    """A reach on a constant slope, bed 0 at the weir (east end)."""
    reach = Reach(length, cells, np.zeros(cells), manning_n=n, width=width,
                  west=("inflow", inflow), east=("weir", crest_height, 1.0))
    reach.bed = slope * (length - reach.x)
    return reach


def steady_weir_run(q, crest_height, length, cells, order=2, tol=1e-8):
    """Steady flow upstream of a weir on the illustrative canal, by the solver.

    tol = 1e-8 m/s: in deep, slow pools the second-order residual stops falling
    near 6e-9 m/s (the limiter switches on an almost flat surface), long after
    the backwater length has settled to 0.1 m.
    """
    reach = sloped_reach(length, cells, crest_height, q)
    h_n = normal_depth(q, SLOPE, MANNING_N, WIDTH)
    start = np.maximum(h_n, crest_height + weir_head(q) - reach.bed)
    model = Model([reach], [start], [np.full(cells, q)], order=order)
    model.run_to_steady(tol=tol)
    return reach, model