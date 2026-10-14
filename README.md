# Canal weirs

How far upstream does an overflow weir raise a peatland canal, and how does its outflow change as the level falls? A 1-D shallow-water solver, written from scratch and verified against exact answers.

Kalimantan peat · 01 of 11

## Setup

    git clone https://github.com/adibatic/kalimantan-peat-01-canal-weirs.git
    cd kalimantan-peat-01-canal-weirs
    conda env create -f environment.yml
    conda activate peat01
    python -m pip install --no-deps -e .
    pytest

## Environment

| Item | Value |
|---|---|
| OS | Microsoft Windows 11 Home Single Language |
| CPU | AMD Ryzen 7 6800HS Creator Edition |
| RAM | 28 GB |
| conda | 26.7.2 |
| Python | 3.14.8 |
| NumPy | 2.5.4 |
| SciPy | 1.18.1 |
| Matplotlib | 3.11.2 |
| pytest | 9.1.1 |

## Agent log

| Date | What the agent did wrong | How I caught it | Rule added to CLAUDE.md |
|---|---|---|---|

## Problem

A rectangular canal drains a peat dome. A fixed-crest overflow weir holds water back.

1. How far upstream does the weir raise the water? Measured as the length over which the depth is more than 5 cm above normal depth.
2. How does the outflow change as the level falls? With and without the weir, when the inflow stops.

The canal is chosen for illustration, not measured: width 3 m, bed slope 2e-4, Manning n = 0.03, discharges 0.1 to 1.2 m³/s, crests 0.5 to 1.5 m above the bed. Results are properties of this model canal, not of any real canal.

## Sources

☐ = equation not yet checked against my own copy; ☑ = checked.

| | Source | Used for |
|---|---|---|
| ☐ | Chow, V. T. (1959). *Open-Channel Hydraulics.* McGraw-Hill. | Manning's equation, critical flow, gradually varied flow, broad-crested weir |
| ☐ | Chaudhry, M. H. (2008). *Open-Channel Flow* (2nd ed.). Springer. | Saint-Venant equations, dam break |
| ☐ | Ritter, A. (1892). Die Fortpflanzung der Wasserwellen. *Zeitschrift des Vereines Deutscher Ingenieure* 36. | Dam break on a dry bed |
| ☐ | Delestre, O., et al. (2013). SWASHES. *Int. J. Numer. Meth. Fluids* 72(3), 269–300. doi:10.1002/fld.3741 | Ritter and MacDonald-type benchmarks |
| ☐ | MacDonald, I., et al. (1997). Analytic benchmark solutions for open-channel flows. *J. Hydraul. Eng.* 123(11). | Manufactured steady flow with friction |
| ☐ | Harten, A., Lax, P. D., van Leer, B. (1983). *SIAM Review* 25(1). | HLL flux |
| ☐ | Audusse, E., et al. (2004). *SIAM J. Sci. Comput.* 25(6). | Hydrostatic reconstruction |
| ☐ | Toro, E. F. (2001). *Shock-Capturing Methods for Free-Surface Shallow Flows.* Wiley. | MUSCL, minmod, wave speeds |

## Progress

Self-paced. Target: the whole project in about a week.

- [x] Milestone 1: Foundation. The repository installs as a package, pytest runs, the rules are set.
- [ ] Milestone 2: Problem defined. Question, sources, exact answers and checks written before any code.
- [ ] Milestone 3: Hydraulics. Manning, critical depth and the weir law, tested.
- [ ] Milestone 4: Exact answers. Ritter, a manufactured steady flow, the GVF profile and the pool drawdown, each tested against its equations.
- [ ] Milestone 5: Solver. Well-balanced finite volumes, two schemes, weirs as boundaries, volume conserved.
- [ ] Milestone 6: Verification. Observed orders 1 and 2, the dam break, the weir rating, the backwater length against GVF.
- [ ] Milestone 7: The answer. Backwater lengths, outflow as the level falls, the dry-spell recession.
- [ ] Milestone 8: Package the result. Animation, CI badge, front page, v1.0 release.