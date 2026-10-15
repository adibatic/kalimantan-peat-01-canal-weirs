"""Animation: a storm, then a dry spell, on a canal with three weirs and without.

A synthetic inflow (not measured) passes a 9 km canal. With weirs, each
reach keeps a pool when the inflow stops; without, the canal drains.
Writes results/animations/weir_cascade.gif, or media/weir_cascade.mp4 with --reel.
"""
from __future__ import annotations

import argparse

import numpy as np
from common import (ANIMATIONS, BLUE, INK, INK_2, MANNING_N, MEDIA, MUTED, ORANGE, PEAT,
                    ROOT, SLOPE, WATER, WIDTH, plt)
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

from canalweir.hydraulics import normal_depth, weir_head
from canalweir.solver import Model, Reach

REACHES, REACH_LENGTH, CREST, CELL = 3, 3000.0, 0.8, 50.0
DAYS, FRAMES = 3.0, 90
DAY = 86400.0


def inflow(t):
    """Synthetic hydrograph (m^3/s): a fading base flow and one storm on day 1."""
    d = t / DAY
    return 0.15 * np.exp(-d / 1.5) + 1.0 * np.exp(-(((d - 1.0) / 0.25) ** 2))


def bed_of(x_global):
    return SLOPE * (REACHES * REACH_LENGTH - x_global)


def build(with_weirs):
    """Weired canal: three reaches in series. Plain canal: one reach, critical outlet."""
    cells = int(REACH_LENGTH / CELL)
    q0 = inflow(0.0) / WIDTH
    h_n = normal_depth(q0, SLOPE, MANNING_N, WIDTH)
    reaches, depths = [], []
    count = REACHES if with_weirs else 1
    length = REACH_LENGTH if with_weirs else REACH_LENGTH * REACHES
    for k in range(count):
        n = cells if with_weirs else cells * REACHES
        west = ("inflow", lambda t: inflow(t) / WIDTH) if k == 0 else ("upstream",)
        r = Reach(length, n, np.zeros(n), manning_n=MANNING_N, width=WIDTH, west=west,
                  east=("weir", 0.0, 1.0))
        r.offset = k * length
        r.bed = bed_of(r.offset + r.x)
        crest = r.bed[-1] - 0.5 * r.dx * SLOPE + (CREST if with_weirs else 0.0)
        r.east = ("weir", crest, 1.0)
        reaches.append(r)
        depths.append(np.maximum(h_n, crest + weir_head(q0) - r.bed))
    model = Model(reaches, depths, [np.full(r.cells, q0) for r in reaches])
    # settle to the steady state of the starting inflow, then run the hydrograph
    for r in reaches[:1]:
        r.west = ("inflow", q0)
    model.run_to_steady(tol=1e-9)
    reaches[0].west = ("inflow", lambda t: inflow(t) / WIDTH)
    model = Model(reaches, model.h, model.q)
    return model


def simulate(with_weirs):
    model = build(with_weirs)
    frames = {"t": [], "level": [], "outflow": [], "volume": []}

    def record(m):
        frames["t"].append(m.t)
        frames["level"].append(np.concatenate([r.bed + h for r, h in zip(m.reaches, m.h)]))
        frames["outflow"].append(m.outflows[-1])
        frames["volume"].append(m.volume())

    model.run(DAYS * DAY, callback=record, every=DAYS * DAY / FRAMES)
    x = np.concatenate([r.offset + r.x for r in model.reaches]) / 1000.0
    crests = [(r.offset + r.length, r.east[1]) for r in model.reaches]
    return x, {k: np.array(v) for k, v in frames.items()}, crests


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reel", action="store_true", help="write an MP4 to media/")
    args = parser.parse_args()

    xw, weired, crests = simulate(True)
    xp, plain, _ = simulate(False)
    print(f"frames: {len(weired['t'])} with weirs, {len(plain['t'])} without")
    print(f"water kept after {DAYS:.0f} days: with weirs {weired['volume'][-1]:.0f} m3, "
          f"without {plain['volume'][-1]:.0f} m3")

    x_line = np.linspace(0, REACHES * REACH_LENGTH, 400)
    fig = plt.figure(figsize=(8.0, 4.8))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.5, 1.0), hspace=0.55, wspace=0.28)
    ax = fig.add_subplot(grid[0, :])
    ax_q = fig.add_subplot(grid[1, 0])
    ax_v = fig.add_subplot(grid[1, 1])

    bed = bed_of(x_line)
    ax.fill_between(x_line / 1000, bed - 0.4, bed, color=PEAT, lw=0)
    for xe, crest in crests:
        ax.plot([xe / 1000, xe / 1000], [bed_of(xe), crest], color=INK, lw=4,
                solid_capstyle="butt")
    fill = [ax.fill_between(xw, bed_of(xw * 1000), weired["level"][0], color=WATER, lw=0)]
    (line_w,) = ax.plot(xw, weired["level"][0], color=BLUE, lw=2, label="with 3 weirs")
    (line_p,) = ax.plot(xp, plain["level"][0], color=ORANGE, lw=1.6, ls="--",
                        label="no weirs")
    clock = ax.text(0.01, 0.06, "", transform=ax.transAxes, color=INK, fontsize=10,
                    fontweight="bold")
    ax.set(xlim=(0, REACHES * REACH_LENGTH / 1000), ylim=(-0.4, bed[0] + 1.6),
           xlabel="distance along the canal (km), flow to the right",
           ylabel="elevation (m)", title="Storm, then a dry spell: what the weirs keep")
    ax.legend(loc="upper right", ncol=2, fontsize=9)

    days = weired["t"] / DAY
    ax_q.plot(days, inflow(weired["t"]), color=MUTED, lw=1.4, label="inflow")
    ax_q.plot(days, weired["outflow"], color=BLUE, lw=1.8, label="out, weirs")
    ax_q.plot(plain["t"] / DAY, plain["outflow"], color=ORANGE, lw=1.6, ls="--",
              label="out, no weirs")
    ax_q.set(xlabel="day", ylabel="discharge (m$^3$/s)", xlim=(0, DAYS))
    ax_q.legend(loc="upper right", fontsize=7.5)
    ax_v.plot(days, weired["volume"] / 1000, color=BLUE, lw=1.8)
    ax_v.plot(plain["t"] / DAY, plain["volume"] / 1000, color=ORANGE, lw=1.6, ls="--")
    ax_v.set(xlabel="day", ylabel="water stored (1000 m$^3$)", xlim=(0, DAYS),
             ylim=(0, None))
    cursors = [a.axvline(0, color=INK_2, lw=1) for a in (ax_q, ax_v)]

    def draw(i):
        fill[0].remove()
        fill[0] = ax.fill_between(xw, bed_of(xw * 1000), weired["level"][i], color=WATER,
                                  lw=0)
        line_w.set_ydata(weired["level"][i])
        line_p.set_ydata(plain["level"][min(i, len(plain["level"]) - 1)])
        clock.set_text(f"day {days[i]:.1f}   inflow {inflow(weired['t'][i]):.2f} m$^3$/s")
        for c in cursors:
            c.set_xdata([days[i], days[i]])
        return line_w, line_p, clock, *cursors

    anim = FuncAnimation(fig, draw, frames=len(days), interval=1000 / 15)
    if args.reel:
        MEDIA.mkdir(exist_ok=True)
        path = MEDIA / "weir_cascade.mp4"
        anim.save(path, writer=FFMpegWriter(fps=30, bitrate=4000), dpi=200)
    else:
        ANIMATIONS.mkdir(parents=True, exist_ok=True)
        path = ANIMATIONS / "weir_cascade.gif"
        anim.save(path, writer=PillowWriter(fps=15), dpi=90)
    print("wrote", path.relative_to(ROOT).as_posix())


if __name__ == "__main__":
    main()