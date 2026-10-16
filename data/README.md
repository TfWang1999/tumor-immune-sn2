# Source data and panel mapping

The numerical tables are copies of the verified manuscript results. Only the
rows and branches needed by the five retained figures and Table 1 are included.
Compact numerical seeds support optional regeneration; they are starting
guesses, not substitutes for evaluating closure and stability residuals.

| Panel/item | Files relative to `data/` |
|---|---|
| Figure 1a–d | `fits/fits.mat`, the corresponding `fits/*_curves.csv`; `observations.csv` and `parameters.json` for reintegration |
| Figure 2a–d | `equilibria/global_equilibrium_branches.csv`, `ablation_summary.csv`, `ablation_landmarks.csv`, `upper_fold_checks.csv` |
| Figure 3a–c | `sn2/sn2_weak_recruitment.csv`; asymptotic constants are evaluated by the plotting script from `parameters.json` |
| Figure 4a | `continuation/branches.csv`, `landmarks.csv`, `trajectories/cycle_alpha_4_metrics.csv`, `cycle_alpha_4p2_metrics.csv` |
| Figure 4b | `periodic/branches.csv`, `equilibrium_overlay.csv`, `landmarks.csv`, `fold_verified.json`, `hopf_reference.json` |
| Figure 4c | `periodic/stable_0p96925_corrected.csv`, `saddle_0p96925_corrected.csv` and companion JSONs; `equilibria/common_parameter_equilibria.csv` |
| Figure 4d,e | `planes/comparison.csv`: only the recruitment–expansion and activation–expansion planes, split by family and direction; `valid` controls display |
| Figure 4f | `periodic/branch_metrics.csv`, `fold_verified.json`, `hopf_reference.json` |
| Figure 5a–d | `timeseries/Figure_05_timeseries.csv`, `Figure_05_parameters.json`; CSV contains visible early windows, JSON contains terminal diagnostics and B1 relaxation envelope |
| Table 1 | `sn2/decay_interaction_summary.csv`, `decay_interaction_effects.csv` |
| Recruitment-path widths mentioned in Results | `equilibria/verified_windows.csv` (four recruitment slices only) |

`cycle_seeds.json`, `branch_seeds.json`, `table_seeds.json` and
`recruitment_seeds.json` store small scalar/state initial guesses copied from
the existing verified periodic solutions. They remove the need to distribute
raw AUTO profiles and historical continuation directories. The default
Figure 4 generator recomputes orbit closures and transverse multipliers from
these seeds; `--auto` starts the displayed continuation branches from the
model equilibria without archived restart solutions.

## Four experimental control datasets

The observations are the control tumour-volume series used in the manuscript,
not new experimental measurements. The third observation column retains the
historical name `sd_mm3`; the MC38 source reports **SEM**, which is preserved
without conversion. Other supplied error bars follow their source studies.

- **LLC:** Chen Z, Ji W, Feng W, Cui J, Wang Y, Li F, et al. PTPRT loss
  enhances anti-PD-1 therapy efficacy by regulation of STING pathway in
  non-small cell lung cancer. *Science Translational Medicine*
  2024;16(763):eadl3598. https://doi.org/10.1126/scitranslmed.adl3598
- **MC38:** Zhang SL, Mao YQ, Zhang ZY, Li ZM, Kong CY, Chen HL, et al.
  Pectin supplement significantly enhanced the anti-PD-1 efficacy in
  tumor-bearing mice humanized with gut microbiota from patients with
  colorectal cancer. *Theranostics* 2021;11(9):4155–4170.
  https://doi.org/10.7150/thno.54476
- **A375:** Rosato RR, D'ávila-González D, Choi DS, Qian W, Chen W,
  Kozielski AJ, et al. Evaluation of anti-PD-1-based therapy against
  triple-negative breast cancer patient-derived xenograft tumors engrafted
  in humanized mouse models. *Breast Cancer Research* 2018;20(1):108.
  Figure 7c, as identified by the manuscript author.
  https://doi.org/10.1186/s13058-018-1037-4
- **CT26:** Liu C, Liu R, Wang B, Lian J, Yao Y, Sun H, et al. Blocking
  IL-17A enhances tumor response to anti-PD-1 immunotherapy in microsatellite
  stable colorectal cancer. *Journal for ImmunoTherapy of Cancer*
  2021;9(1):e001895. https://doi.org/10.1136/jitc-2020-001895

`fits.mat` is a plotting cache with the four observations and fitted RMSEs.
The dense fitted curves can be regenerated with the saved parameter estimates
and the copied MATLAB model; parameter optimisation is not part of this release.

## Interpretation

M1/M2 omit recruitment and have linear/saturated DC activation, respectively;
M3/M4 include recruitment with the corresponding activation response.
The low-signal activation slope is matched. Table entries are actual
growth-time parameter widths `W_alpha=alpha_H-alpha_LPC`, with display
factor `10^4`. Independent stability checks sampled the original windows;
the public width generator recomputes their Hopf and periodic-fold endpoints.

The general theorems and the independent reference-parameter interval
applications remain in the manuscript and Supplementary Information. Their
historical CAPD certificates and exploratory numerical runs are not inputs
to these five figures and are not bundled here.
