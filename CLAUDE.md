# Project rules

Canal weirs (Kalimantan peat · 01 of 11): a 1-D shallow-water solver for a canal
with overflow weirs, Python. Package in src/canalweir, tests in tests/,
experiments in experiments/, results in results/.

## Units and numbers
- SI units everywhere: m, s, m^3/s. g = 9.81 m/s^2.
- The canal in experiments/common.py (width 3 m, n = 0.03, slope 2e-4, discharges
  0.1 to 1.2 m^3/s, crests 0.5 to 1.5 m) is chosen for illustration, not measured.
  Never present a result from it as a property of a real peat canal.
- Never change a test tolerance, a cell size, the CFL number, a steady-state tolerance,
  or a canal parameter without saying so in the commit message.
- Never loosen a test to make it pass. Find the cause.
- Equations, model constants and reference values must come from a source I have read
  (README, Sources). Never invent or "recall" one, or an equation or page number.
- No placement, forecast or verification claim about a real canal before validation
  on held-out field data.

## Commands
Windows 11, PowerShell 7, conda (Miniforge). Commands are PowerShell.
- Environment (once): conda env create -f environment.yml
- Activate: conda activate peat01
- Install the package (once): python -m pip install --no-deps -e .
- Tests: pytest
- Experiment: python experiments/<name>.py (dam_break, steady_convergence, backwater, drawdown, weir_cascade)
- Video: python experiments/weir_cascade.py --reel

## Checks that must pass before a result is used
- pytest passes.
- Lake at rest over a bump: |q| and level error below 1e-12.
- Volume balance (stored change = inflow - outflow): relative error below 1e-13.
- Observed order on the steady benchmark: 1 +/- 0.1 (order=1), 2 +/- 0.1 (order=2).
- Steady runs stop at |dh/dt| below 1e-8 m/s; solver and GVF backwater lengths agree within 1 %.

## Workflow
- I build this project step by step and type every file myself. Explain, review, and find bugs. Do not create or edit files unless I ask for it in that message.
- Do not commit or push unless I ask. If a request is ambiguous ("proceed", "go ahead"), ask what I mean before changing anything.
- Commits follow Conventional Commits 1.0.0 (feat, fix, perf, style, test, docs, build, ci, chore, refactor).
- media/ is never committed.
- Explain every change you make and why.