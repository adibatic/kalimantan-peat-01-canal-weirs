"""The answer, part 1: how far upstream does a weir raise the water?

For each discharge and crest height:
- the gradually varied flow (GVF) profile gives the backwater length (reference),
- the finite-volume solver run to steady state gives it independently,
- the level-pool estimate (crest + head - normal depth) / slope is the quick rule.
"""
from __future__ import annotations

import csv

import numpy as np
from common import (BLUE, BLUE_RAMP, CREST_HEIGHTS, DISCHARGES, INK, INK_2, MANNING_N,
                    MUTED, ORANGE, PEAT, SLOPE, TABLES, THRESHOLD, WATER, WIDTH, plt,
                    save, steady_weir_run)

from canalweir.exact import backwater_length, backwater_profile, level_pool_length
from canalweir.hydraulics import normal_depth, weir_head

GVF_LENGTH = 40_000.0   # m, long enough for every case to return to normal depth
CELL_SIZE = 40.0        # m, finite-volume cell size for the cross-check


def gvf_length(q, crest):
    x, h = backwater_profile(q, SLOPE, MANNING_N, crest, GVF_LENGTH, WIDTH)
    return backwater_length(x, h, normal_depth(q, SLOPE, MANNING_N, WIDTH),
                            GVF_LENGTH, THRESHOLD)


def main():
    rows = []
    for Q in DISCHARGES:
        q = Q / WIDTH
        h_n = normal_depth(q, SLOPE, MANNING_N, WIDTH)
        for crest in CREST_HEIGHTS:
            head = weir_head(q)
            L_gvf = gvf_length(q, crest)
            L_pool = max(level_pool_length(q, SLOPE, MANNING_N, crest, WIDTH), 0.0)
            L_fv = float("nan")
            if L_gvf > 0.0:
                length = 1.3 * L_gvf + 1000.0
                reach, model = steady_weir_run(q, crest, length, int(length / CELL_SIZE))
                L_fv = backwater_length(reach.x, model.h[0], h_n, length, THRESHOLD)
            rows.append({"Q_m3s": Q, "crest_m": crest, "normal_depth_m": h_n,
                         "weir_head_m": head, "rise_at_weir_m": crest + head - h_n,
                         "L_gvf_m": L_gvf, "L_solver_m": L_fv,
                         "solver_minus_gvf_m": L_fv - L_gvf, "L_level_pool_m": L_pool,
                         "gvf_over_pool": L_gvf / L_pool if L_pool > 0 else float("nan")})
            r = rows[-1]
            print(f"Q {Q:.1f}  crest {crest:.1f}  h_n {h_n:.3f}"
                  f"  rise {r['rise_at_weir_m']:+.3f}"
                  f"  L_gvf {L_gvf:7.0f}  L_solver {L_fv:7.0f}  L_pool {L_pool:6.0f}")
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "backwater.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    profile_figure()
    length_figure()


def profile_figure(Q=0.3, crest=1.0):
    q = Q / WIDTH
    h_n = normal_depth(q, SLOPE, MANNING_N, WIDTH)
    L_gvf = gvf_length(q, crest)
    length = 1.3 * L_gvf + 1000.0
    reach, model = steady_weir_run(q, crest, length, int(length / CELL_SIZE))
    x, h = backwater_profile(q, SLOPE, MANNING_N, crest, length, WIDTH)
    up = (length - x) / 1000.0                     # km upstream of the weir
    bed = SLOPE * (length - x)
    up_fv = (length - reach.x) / 1000.0
    pool_level = crest + weir_head(q)
    L_pool = level_pool_length(q, SLOPE, MANNING_N, crest, WIDTH)

    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.fill_between(up, bed - 0.3, bed, color=PEAT, alpha=0.9, lw=0, label="canal bed")
    ax.fill_between(up, bed, bed + h, color=WATER, lw=0, alpha=0.8)
    ax.plot(up, bed + h_n, color=MUTED, lw=1.4, ls="--", label="no weir (normal depth)")
    ax.plot(up, bed + h, color=BLUE, label="with weir, GVF profile")
    ax.plot(up_fv[::4], (reach.bed + model.h[0])[::4], "o", ms=3.5, color=ORANGE,
            label="with weir, solver")
    ax.plot([0, L_pool / 1000], [pool_level, pool_level], color=INK_2, lw=1.2, ls=":",
            label="level-pool estimate")
    ax.bar([0], [crest], width=0.12, bottom=0, color=INK, align="center")
    ax.annotate("", xy=(L_gvf / 1000, pool_level + 0.25), xytext=(0, pool_level + 0.25),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
    ax.text(L_gvf / 2000, pool_level + 0.3,
            f"raised more than {THRESHOLD * 100:.0f} cm over {L_gvf / 1000:.2f} km",
            ha="center", color=INK, fontsize=9)
    ax.set_xlim(-0.15, min(length, 1.15 * L_gvf) / 1000)
    ax.set_ylim(-0.3, pool_level + 0.55)
    ax.set(xlabel="distance upstream of the weir (km)", ylabel="elevation (m)",
           title=f"Backwater of a {crest:.1f} m weir at Q = {Q} m$^3$/s")
    ax.legend(loc="lower right", fontsize=8.5, ncol=2, frameon=True, facecolor="white",
              edgecolor="none", framealpha=0.92)
    save(fig, "backwater_profile.png")


def length_figure():
    crests = np.linspace(0.0, 2.0, 81)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for Q, color in zip(DISCHARGES, BLUE_RAMP):
        q = Q / WIDTH
        L = np.array([gvf_length(q, c) for c in crests]) / 1000.0
        pool = np.array([max(level_pool_length(q, SLOPE, MANNING_N, c, WIDTH), 0.0)
                         for c in crests]) / 1000.0
        ax.plot(crests, L, color=color, label=f"Q = {Q} m$^3$/s")
        ax.plot(crests, pool, color=color, lw=1.0, ls="--")
    ax.plot([], [], color=MUTED, lw=1.0, ls="--", label="level-pool estimate")
    ax.set(xlabel="crest height above the bed (m)",
           ylabel=f"length raised > {THRESHOLD * 100:.0f} cm (km)",
           title="How far upstream a weir raises the water")
    ax.set_xlim(0, 2.0)
    ax.legend(loc="upper left", fontsize=8.5)
    save(fig, "backwater_length.png")


if __name__ == "__main__":
    main()