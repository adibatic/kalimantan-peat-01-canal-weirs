"""The answer, part 2: how does the outflow change as the level falls?

1. Rating curves: outflow against water depth at the canal end, with and
   without a weir (no weir = free overfall at bed level, crest 0).
2. Recession: inflow stops; compare the volume a canal keeps with and
   without a weir over three days.
3. Check: a short flat pool emptying over a weir against the exact lumped
   solution H(t) = H0 / (1 + C sqrt(H0) t / (2 L))^2.
"""
from __future__ import annotations

import csv

import numpy as np
from common import (BLUE, INK, INK_2, MANNING_N, MUTED, ORANGE, SLOPE, TABLES, WIDTH, plt,
                    save, sloped_reach, steady_weir_run)

from canalweir.exact import drawdown_head
from canalweir.hydraulics import WEIR_COEFF, normal_depth, weir_unit_discharge
from canalweir.solver import Model, Reach

Q_WET, CREST, LENGTH, CELLS, DAYS = 0.6, 1.0, 4000.0, 100, 3.0


def rating_figure():
    depth = np.linspace(0.0, 1.5, 400)
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.plot(depth, weir_unit_discharge(depth, 0.0) * WIDTH, color=ORANGE,
            label="no weir (critical flow at the end)")
    ax.plot(depth, weir_unit_discharge(depth, CREST) * WIDTH, color=BLUE,
            label=f"{CREST:.1f} m weir")
    ax.axvline(CREST, color=MUTED, lw=1, ls="--")
    ax.text(CREST + 0.02, 4.6, "crest", color=MUTED, fontsize=9)
    ax.set(xlabel="water depth at the canal end (m)", ylabel="outflow (m$^3$/s)",
           title="Outflow as the level falls", xlim=(0, 1.5), ylim=(0, 5))
    ax.legend(loc="upper left")
    save(fig, "rating_curves.png")


def pool_volume_below_crest(crest):
    """Volume a weir keeps when the level has fallen to its crest:
    width times the integral of max(crest - z, 0) along the reach, z = S0 (L - x)."""
    wet = min(LENGTH, crest / SLOPE)                 # length of bed below the crest
    return WIDTH * (crest * wet - SLOPE * wet**2 / 2.0)


def recession(crest):
    """Steady at Q_WET, then the inflow stops. Returns times (s), volumes, outflows."""
    q = Q_WET / WIDTH
    if crest > 0:
        reach, model = steady_weir_run(q, crest, LENGTH, CELLS)
        h, qq = model.h[0], model.q[0]
    else:
        reach = sloped_reach(LENGTH, CELLS, 0.0, q)
        h_n = normal_depth(q, SLOPE, MANNING_N, WIDTH)
        model = Model([reach], [np.full(CELLS, h_n)], [np.full(CELLS, q)])
        model.run_to_steady(tol=1e-10)
        h, qq = model.h[0], model.q[0]
    reach.west = ("inflow", 0.0)
    model = Model([reach], [h], [qq])
    times, volumes, outflows = [], [], []

    def record(m):
        times.append(m.t)
        volumes.append(m.volume())
        outflows.append(m.outflows[0])

    model.run(DAYS * 86400.0, callback=record, every=600.0)
    return np.array(times), np.array(volumes), np.array(outflows)


def recession_figure():
    rows = []
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.6))
    cases = ((0.0, ORANGE, "no weir"), (CREST, BLUE, f"{CREST:.1f} m weir"))
    for crest, color, label in cases:
        t, v, out = recession(crest)
        days = t / 86400.0
        ax1.plot(days[1:], out[1:], color=color, label=label)
        ax2.plot(days, v, color=color, label=label)
        for d in (0.25, 1.0, 3.0):
            k = np.argmin(np.abs(days - d))
            rows.append({"case": label, "day": d, "outflow_m3s": out[k],
                         "volume_m3": v[k], "volume_percent": 100 * v[k] / v[0]})
            print(f"{label:12s} day {d:4.2f}  outflow {out[k]:.2e} m3/s  "
                  f"volume {v[k]:7.0f} m3 ({100 * v[k] / v[0]:5.1f} %)")
    pool = pool_volume_below_crest(CREST)
    ax2.axhline(pool, color=INK_2, lw=1.0, ls=":")
    ax2.text(DAYS, pool * 1.03, "pool below the crest, exact", ha="right", va="bottom",
             color=INK_2, fontsize=8.5)
    print(f"pool below the crest, exact: {pool:.0f} m3")
    ax1.set(xlabel="days after the inflow stops", ylabel="outflow (m$^3$/s)",
            title="Outflow during a dry spell", yscale="log", ylim=(1e-4, 1))
    ax2.set(xlabel="days after the inflow stops", ylabel="water in the canal (m$^3$)",
            title="Storage during a dry spell", ylim=(0, None))
    ax1.legend(loc="upper right")
    ax2.legend(loc="upper right")
    save(fig, "recession.png")
    with open(TABLES / "recession.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def pool_check():
    """Flat pool, wall upstream, weir downstream, no friction, no inflow."""
    length, cells, crest, head0 = 100.0, 50, 1.5, 0.1
    reach = Reach(length, cells, np.zeros(cells), east=("weir", crest, 1.0))
    model = Model([reach], [np.full(cells, crest + head0)], [np.zeros(cells)])
    t_end = 6.0 * length / (WEIR_COEFF * np.sqrt(head0))
    times, heads = [], []
    model.run(t_end, callback=lambda m: (times.append(m.t),
                                         heads.append(m.volume() / length - crest)),
              every=t_end / 60)
    times, heads = np.array(times), np.array(heads)
    exact = drawdown_head(times, head0, length)
    err = np.max(np.abs(heads - exact)) / head0
    balance = model.volume() + model.outflow_volume - (crest + head0) * length
    print(f"pool check: max |H - H_exact| / H0 = {err:.4f}, "
          f"volume balance {balance:.1e} m2")
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    tt = np.linspace(0, t_end, 400)
    ax.plot(tt / 60, drawdown_head(tt, head0, length) * 100, color=INK, lw=1.2,
            label="exact lumped solution")
    ax.plot(times / 60, heads * 100, "o", ms=4, color=BLUE, label="solver, mean level")
    ax.set(xlabel="time (min)", ylabel="head over the crest (cm)",
           title="A flat pool emptying over a weir")
    ax.legend(loc="upper right")
    ax.text(0.98, 0.55, f"max difference {100 * err:.1f} % of H$_0$",
            transform=ax.transAxes, ha="right", color=INK_2, fontsize=9)
    save(fig, "pool_drawdown.png")
    return err


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    rating_figure()
    pool_check()
    recession_figure()


if __name__ == "__main__":
    main()