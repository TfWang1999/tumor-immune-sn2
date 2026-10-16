# Antigen-driven tumour–immune multistability: minimal reproduction code

Code and source data for the five main figures and Table 1 of
**Antigen-driven restructuring of multistability in a tumour–immune–dendritic-cell model**.

This standalone collection copies the existing plotting, scalar-equilibrium,
shooting and integration routines, with local paths adapted. It contains only
the displayed parameter planes, selected trajectories and mechanism windows.
Plotting uses the supplied results and does not launch simulations or fitting.

## Replot the five main figures

Use Python 3.11 or newer and MATLAB (tested with MATLAB R2025a).
AUTO07p is needed only for optional continuation, not for plotting.

```sh
python -m pip install -r requirements.txt
matlab -batch "figure_01_fits"
python plot_figure_02.py
python plot_figure_03.py
python plot_figure_04.py
python plot_figure_05.py
python build_table_01.py
```

Run these commands from this directory. Every plot is saved to `outputs/` as
PDF, SVG and PNG. Table 1 is written as CSV and Markdown. The plotting scripts
preserve stability-based branch splitting: separate branches and opposite
stability intervals are not joined by a single line.

| Main item | Content | Supplied data | Data-generation entry |
|---|---|---|---|
| Figure 1 | Four conditional model fits, RMSE and R² | `data/fits/`, `observations.csv`, `parameters.json` | `generate_figure_01.m` |
| Figure 2 | Global equilibria and matched high-tumour local enlargements | `data/equilibria/` | `generate_figure_02_03.py` |
| Figure 3 | SN2 inverse-recruitment law, constant term and M4–M3 shift | `data/sn2/sn2_weak_recruitment.csv` | `generate_figure_02_03.py` |
| Figure 4 | Two equilibrium/cycle slices, periodic companions, two parameter planes and Floquet multipliers | `data/continuation/`, `periodic/`, `planes/`, `trajectories/` | `generate_figure_04.py`, optional `--auto` |
| Figure 5 | Four slices, two initial states per attracting outcome | `data/timeseries/` | `generate_figure_05.py` |
| Table 1 | M1–M4 window widths at three DC decay rates | `data/sn2/decay_interaction_*.csv` | `generate_table_01.py` |

## Regenerate the numerical data

These commands overwrite **only this collection's `data/` files**. Use a fresh
clone if you want to retain the supplied publication data. Sampling and AUTO
branch point spacing can vary slightly across platforms and solver versions.

```sh
matlab -batch "generate_figure_01"
python generate_figure_02_03.py
python generate_figure_04.py
python generate_figure_05.py
python generate_table_01.py
```

Figure 1 reintegrates the model at the saved fitted estimates and recomputes
RMSE and R²; it does not repeat parameter optimisation. The original fit is
conditional on fixed immune coefficients, as described in the manuscript.
Figure 2–3 generation uses the scalar equilibrium and fold equations. The
default Figure 4 generator independently refines the two periodic companions,
computes transverse Floquet multipliers along the saved branch seeds and
refines the cycle fold. Table 1 uses compact fold seeds to solve the periodic
closure, phase and unit-transverse-multiplier equations for twelve slices;
it also regenerates the four recruitment-path widths discussed in the text.

Figure 5 retains the original long integration horizons, up to
`tau = 220000`, to distinguish slow damping from sustained cycles. It tests
all sixteen trajectories and records terminal outcomes in JSON, while the
CSV stores only the visible early windows. B1's long-time maxima/minima are
retained in the JSON used by its relaxation inset. These computations can
take several minutes; plotting the supplied results takes seconds.

### Optional AUTO07p continuation

Install AUTO07p with its normal Python interface and a Fortran compiler.
Only `auto/model.f90` and `auto/continue.auto` are supplied: they regenerate
the two Figure 4 one-parameter slices and the `(sigma, alpha)` and
`(beta, alpha)` planes. No pre-existing AUTO restart files are required.

On Linux/macOS, activate the AUTO environment so `auto` is on `PATH`, then run:

```sh
python generate_figure_04.py --auto
```

For Windows with MSYS2, set `AUTO_MSYS_BASH` to the MSYS `bash.exe` and
`AUTO_DIR` to the **POSIX path** of your AUTO installation, then run the same
command. For example in PowerShell (replace the paths):

```powershell
$env:AUTO_MSYS_BASH = 'D:/msys64/usr/bin/bash.exe'
$env:AUTO_DIR = '/home/yourname/auto-07p'
python generate_figure_04.py --auto
```

AUTO writes raw branches and its log to `_work/auto/`; only the small CSV
exports are copied into `data/`. Remove that task-specific working directory
when changing parameters to start a fresh continuation. It is excluded from
Git, together with generated plots. Both requested parameter planes and the
independent Figure 4 shooting/Floquet calculation were tested locally.

## Coordinates

All dynamical figures use the manuscript's growth-time dimensionless
coordinates `tau=r*t`, `x=C/K`, `y=eta*T/r`, `z=D/K_D` and parameters
`sigma, alpha, delta, beta, a, gamma, kappa`. Figure 1 shows the observed
day/volume coordinates.

The copied AUTO and shooting routines retain the original fixed reference
clock `tau0=r0*t`, with `rho=r/r0` and reference effector state `y0=rho*y`.
Their internal coefficients are dimensionless too. To obtain the plotted
growth-time parameters: `alpha=alpha0/rho`, `sigma=sigma0/rho**2`,
`beta=beta0/rho`, `a=a0/rho`, `gamma=gamma0/rho`, and `tau=rho*tau0`.
Some historical CSV column names include `_per_day`; the plotting scripts
explicitly convert them before plotting. They are not additional model
parameters. `model.py` contains the common vector field and Jacobian.

## Data provenance and scope

See [data/README.md](data/README.md) for panel-level inputs and experimental
source citations, and [TESTED.md](TESTED.md) for the local test record.
This release covers the five main figures and Table 1. Historical exploratory
plots, unused parameter planes, robustness sweeps, fitting optimisation and
independent CAPD benchmark certificate machinery are outside this minimal
figure-reproduction collection.

For citation metadata, see `CITATION.cff`. The repository carries no claim of
a software DOI or archived dataset DOI.
