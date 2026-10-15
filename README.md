<div align="center">

# Canal weirs

**How far upstream does an overflow weir raise a peatland canal,<br>and how does its outflow change as the level falls?**

A 1-D shallow-water solver, written from scratch and checked against exact answers.

[![tests](https://github.com/adibatic/kalimantan-peat-01-canal-weirs/actions/workflows/tests.yml/badge.svg)](https://github.com/adibatic/kalimantan-peat-01-canal-weirs/actions/workflows/tests.yml)
![python](https://img.shields.io/badge/python-3.10%2B-3776ab)
![licence](https://img.shields.io/badge/licence-MIT-2a78d6)
![series](https://img.shields.io/badge/kalimantan--peat-01%20of%2011-6b5444)

Kalimantan peat · **01 of 11** · next: 02 groundwater (coming)

<img src="results/animations/weir_cascade.gif" alt="A 9 km canal with three weirs and without, through a storm and a dry spell. With weirs the canal keeps three pools; without, it drains." width="720">

<sub>A synthetic storm, then six dry days, on the same 9 km canal with three 0.8 m weirs (blue) and without (orange, dashed).<br>After six days the weired canal holds 13,706 m³; the open canal holds 839 m³.</sub>

</div>

---

## Key finding

In the illustrative canal (3 m wide, slope 0.2 m per km, Manning n = 0.03, all chosen, not measured):

1. **A weir only works if its crest plus its overflow head stands above the normal depth.** A 0.5 m weir raises the water by more than 5 cm at 0.1 and 0.3 m³/s, but not at 0.6 m³/s or above: there the canal's own depth already drowns it.
2. **Above that, every metre of crest raises roughly 5 km of canal** (about 1/slope), from 4.0 to 4.7 km for a 1.0 m weir across all four flows.
3. **The level-pool estimate is good in the dry season and short in the wet season.** At 0.1 m³/s it is within 6 % of the computed backwater length. At 0.3 to 1.2 m³/s the water is raised 1.2 to 3.2 times farther than it predicts, because a flowing canal's surface is not flat.
4. **As the level falls, the weir's outflow falls to zero at the crest; the open canal's does not.** After the inflow stops, a 4 km reach with a 1.0 m weir keeps 7,204 m³ after three days (the exact pool below the crest is 7,200 m³). Without the weir, the same reach keeps 14 m³.

These are properties of a model canal. They say nothing about a real canal until it is surveyed (project 03) and the model is validated against field loggers (project 11).

## Results at a glance

Length of canal raised more than 5 cm above normal depth, upstream of one weir. Gradually varied flow (GVF) reference; the level-pool estimate in brackets.

| Discharge (m³/s) | Normal depth (m) | 0.5 m crest | 1.0 m crest | 1.5 m crest |
|---:|---:|---:|---:|---:|
| 0.1 | 0.215 | 1.90 km (1.79) | 4.45 km (4.29) | 6.96 km (6.79) |
| 0.3 | 0.437 | 1.64 km (1.07) | 4.56 km (3.57) | 7.15 km (6.07) |
| 0.6 | 0.696 | **0** (0.22) | 4.70 km (2.72) | 7.56 km (5.22) |
| 1.2 | 1.135 | **0** (0) | 3.95 km (1.23) | 8.15 km (3.73) |

The finite-volume solver, run to steady state with 40 m cells, gives the same lengths within 5 m (0.06 %) in all ten cases where water is raised ([`results/tables/backwater.csv`](results/tables/backwater.csv)).

<p align="center">
<img src="results/figures/backwater_profile.png" width="49%" alt="Water surface upstream of a 1.0 m weir at 0.3 m3/s: GVF profile, solver, normal depth and level-pool estimate">
<img src="results/figures/backwater_length.png" width="49%" alt="Raised length against crest height for four discharges">
</p>

## Verification

Every check compares the solver with an exact answer, never with itself. The tests in [`tests/`](tests) run all of them on every push.

| Check | Exact answer | Result | Criterion |
|---|---|---|---|
| Lake at rest over a bump | velocity stays 0 | below 1e-12 m²/s, both schemes | < 1e-12 |
| Volume, closed canal and two reaches with weirs | stored change = inflow - outflow | relative error below 1e-13 | < 1e-13 |
| Order of accuracy, steady flow with friction | manufactured depth h(x) | observed order **0.999** and **2.000** | 1 ± 0.1, 2 ± 0.1 |
| Weir fails onto a dry bed | Ritter's solution | L1 error 2.5 mm (order 2), 5.2 mm (order 1), 200 cells | < 3 mm, < 6 mm |
| Steady flow over a weir | outflow = inflow, head = rating curve | outflow within 1e-8 | 1e-8 |
| Backwater length | GVF integration (independent method) | within 5 m of up to 8.1 km | < 1 % |
| Pool left when the inflow stops | width × ∫ max(crest - z, 0) dx = 7,200 m³ | 7,204 m³ after 3 days | (converging) |
| Flat pool emptying over a weir | lumped H(t) = H₀ / (1 + C√H₀ t / 2L)² | within 1.0 % of H₀ | model-to-model |

The dam break converges at order 1.00 (second-order scheme) and 0.71 (first order), as expected for a solution with a kink and a dry front. Its volume drifts by 1e-16.

<p align="center">
<img src="results/figures/steady_convergence.png" width="32%" alt="Error against cell size, slopes 1 and 2">
<img src="results/figures/dam_break.png" width="35%" alt="Dam break: solver against Ritter">
<img src="results/figures/pool_drawdown.png" width="31%" alt="Flat pool draining over a weir: solver against lumped exact solution">
</p>

## Outflow as the level falls

<p align="center">
<img src="results/figures/rating_curves.png" width="36%" alt="Outflow against depth with and without a weir">
<img src="results/figures/recession.png" width="62%" alt="Outflow and storage after the inflow stops, with and without a weir">
</p>

The open canal ends in critical flow, so its outflow falls with depth but only reaches zero when the canal is empty. The weir's outflow reaches zero at the crest. That is the feedback a drying peat surface provides on its own, rebuilt in the canal: the drier it gets, the less it loses.

## How it works

**Equations.** Depth h and discharge per metre of width q = hu along the canal, bed elevation z, Manning friction (sources below):

$$h_t + q_x = 0, \qquad q_t + \left(\frac{q^2}{h} + \frac{g h^2}{2}\right)_x = -g h z_x - g h S_f, \qquad S_f = \frac{n^2 q|q|}{h^2 R^{4/3}}$$

**Weir law.** On a broad crest the flow passes through critical depth, two thirds of the head H above the crest. Then

$$q = \sqrt{g\left(\tfrac{2}{3}H\right)^3} = \left(\tfrac{2}{3}\right)^{3/2}\sqrt{g}\, H^{3/2} \approx 1.70\, H^{3/2}\ \text{m}^2/\text{s}$$

ideal, approach velocity neglected, never drowned. "No weir" is the same law with the crest at the bed: critical flow at the canal end.

**Reference profile.** Steady flow turns the momentum equation into the gradually varied flow equation, integrated upstream from the weir:

$$\frac{dh}{dx} = \frac{S_0 - S_f}{1 - q^2/(g h^3)}$$

**Solver** ([`src/canalweir/solver.py`](src/canalweir/solver.py)). Finite volumes with the HLL flux, hydrostatic reconstruction of the bed (a lake at rest stays exactly at rest), and two schemes: first order with forward Euler, or second order with minmod-limited reconstruction and Heun's method. Manning friction is semi-implicit inside each stage, so a steady state does not depend on the time step. Reaches are chained by their weirs: the outflow of one is the inflow of the next, inside every stage, so volume is conserved exactly.

**Exact answers** ([`src/canalweir/exact.py`](src/canalweir/exact.py)). Each one is derived in its docstring from the equations above and tested against them in [`tests/test_exact.py`](tests/test_exact.py), so a slip in a derivation fails a test before it can mislead the solver tests.

## Scope and limits

- 1-D, rectangular, prismatic canal. No seepage to or from the peat (project 02), no rain or evaporation on the canal.
- Weir: ideal broad-crested law, discharge coefficient 1, approach velocity neglected, free overflow only (a drowned weir is not modelled).
- The canal, flows, crests and storm are chosen for illustration. Peat canals vary; every length here scales roughly with 1/slope.
- Steady runs stop when |dh/dt| < 1e-8 m/s. In deep, slow pools the second-order residual stalls near 6e-9 m/s (the limiter switches on an almost flat surface); the backwater length has settled to 0.1 m long before.
- Planned but not done: a code-to-code comparison with HEC-RAS.

## Sources

The equations in the code name their source. ☐ means I have not yet checked the equation against my own copy of the source; ☑ means I have. Only the SWASHES entry was checked against a search result; the other bibliographic details (pages, volume) must be checked too.

| | Source | Used for |
|---|---|---|
| ☐ | Chow, V. T. (1959). *Open-Channel Hydraulics.* McGraw-Hill. | Manning's equation, critical flow, gradually varied flow, broad-crested weir |
| ☐ | Chaudhry, M. H. (2008). *Open-Channel Flow* (2nd ed.). Springer. | Saint-Venant equations, dam break |
| ☐ | Ritter, A. (1892). Die Fortpflanzung der Wasserwellen. *Zeitschrift des Vereines Deutscher Ingenieure* 36. | Dam break on a dry bed |
| ☐ | Delestre, O., Lucas, C., Ksinant, P.-A., Darboux, F., Laguerre, C., Vo, T.-N.-T., James, F., Cordier, S. (2013). SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. *Int. J. Numer. Meth. Fluids* 72(3), 269–300. [doi:10.1002/fld.3741](https://doi.org/10.1002/fld.3741), free at [arXiv:1110.0288](https://arxiv.org/abs/1110.0288). | Ritter and MacDonald-type benchmarks |
| ☐ | MacDonald, I., Baines, M. J., Nichols, N. K., Samuels, P. G. (1997). Analytic benchmark solutions for open-channel flows. *J. Hydraul. Eng.* 123(11). | Manufactured steady flow with friction |
| ☐ | Harten, A., Lax, P. D., van Leer, B. (1983). On upstream differencing and Godunov-type schemes for hyperbolic conservation laws. *SIAM Review* 25(1). | HLL flux |
| ☐ | Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004). A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows. *SIAM J. Sci. Comput.* 25(6). | Hydrostatic reconstruction, first and second order |
| ☐ | Toro, E. F. (2001). *Shock-Capturing Methods for Free-Surface Shallow Flows.* Wiley. | MUSCL reconstruction, minmod limiter, wave-speed estimates |
| ☐ | Urzainki, I., et al. (2020). Canal blocking optimization in restoration of drained peatlands. *Biogeosciences* 17. | Context: block placement is already studied |

## Reproduce

With conda (Windows, Linux or macOS; the commands are the same in PowerShell and in bash):

```
git clone https://github.com/adibatic/kalimantan-peat-01-canal-weirs.git
cd kalimantan-peat-01-canal-weirs
conda env create -f environment.yml
conda activate peat01
python -m pip install --no-deps -e .
pytest                                        # 29 tests, about 15 s
python experiments/dam_break.py               # seconds
python experiments/steady_convergence.py      # about 45 s
python experiments/backwater.py               # about 45 s
python experiments/drawdown.py                # about 20 s
python experiments/weir_cascade.py            # about 2 min, the GIF above
python experiments/weir_cascade.py --reel     # MP4 in media/ (needs ffmpeg)
```

Without conda, use `pip install -r requirements.txt` in any Python 3.10 to 3.14 environment instead of the two `conda` lines. The tests run on Linux and Windows with Python 3.10 and 3.14 on every push.

## Repository

```
kalimantan-peat-01-canal-weirs/
├── src/canalweir/
│   ├── hydraulics.py       Manning, critical depth, weir law
│   ├── exact.py            Ritter, manufactured steady flow, GVF profile, pool drawdown
│   └── solver.py           finite-volume solver, reaches chained by weirs
├── tests/                  29 checks against exact answers
├── experiments/            one script per figure, shared settings in common.py
├── results/                figures, tables (CSV), the animation
├── environment.yml         conda environment (Python 3.14)
├── requirements.txt        the same libraries for pip
├── CLAUDE.md               rules for the AI agent
└── .github/workflows/      pytest on Linux and Windows, Python 3.10 and 3.14
```

## Environment

The machine that produced `results/`:

| Item | Value |
|---|---|
| Python | 3.14.8 |
| NumPy | 2.5.4 |
| SciPy | 1.18.1 |
| Matplotlib | 3.11.2 |
| pytest | 9.1.1 |

## Agent log

| Date | What the agent did wrong | How I caught it | Rule added to CLAUDE.md |
|---|---|---|---|

## Progress

Self-paced. Target: the whole project in about a week.

- [x] Milestone 1: Foundation. The repository installs as a package, pytest runs, the rules are set.
- [x] Milestone 2: Problem defined. Question, sources, exact answers and checks written before any code.
- [x] Milestone 3: Hydraulics. Manning, critical depth and the weir law, tested.
- [x] Milestone 4: Exact answers. Ritter, a manufactured steady flow, the GVF profile and the pool drawdown, each tested against its equations.
- [x] Milestone 5: Solver. Well-balanced finite volumes, two schemes, weirs as boundaries, volume conserved.
- [x] Milestone 6: Verification. Observed orders 1 and 2, the dam break, the weir rating, the backwater length against GVF.
- [x] Milestone 7: The answer. Backwater lengths, outflow as the level falls, the dry-spell recession.
- [x] Milestone 8: Package the result. Animation, CI badge, front page, v1.0 release.

## Licence and credits

Code under the [MIT licence](LICENSE). Written by Adriel Imaran Santoso, with Claude Code as a reviewing agent (rules in [CLAUDE.md](CLAUDE.md)). Exact solutions are the work of the authors cited above.