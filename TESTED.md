# Local reproduction test — 16 October 2026

Environment: Windows, Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1,
Matplotlib 3.10.9, MATLAB R2025a, AUTO07p through MSYS2.

The release directory was copied outside the original project tree. All
five plotting entries and the Table 1 builder ran successfully from that
isolated copy, using only its bundled `data/`. They generated five sets of
PDF/SVG/PNG files and the Markdown/CSV table. No original project imports or
absolute author-specific paths are needed by the plotting code.

The numerical generation entries were also exercised locally:

- **Figure 1:** reintegrated all four saved parameter estimates. RMSE values
  were 5.91843898, 34.69233645, 84.27914531 and 8.72264883 mm³;
  R² values were 0.99964364, 0.99009767, 0.98078196 and 0.99951754.
- **Figures 2–3:** regenerated the equilibrium tables and all sixteen
  weak-recruitment fold samples; both plots rendered from those generated data.
- **Figure 4:** the two one-parameter slices and both displayed parameter
  planes were regenerated with `generate_figure_04.py --auto`, starting from
  model equilibria without the original AUTO restart files. Independent
  correction closed the stable/saddle companions with maximum state residuals
  of approximately 3.2e-11 and 1.9e-11. The fast cycles, Floquet branch,
  refined cycle fold and equilibrium overlay were regenerated; the plot
  rendered from the regenerated continuation/shooting data.
- **Figure 5:** integrated all sixteen initial states to the original terminal
  horizons. Every initial state reached its expected equilibrium or stable
  cycle; final equilibrium distances and late-cycle peak tests passed.
  In particular B1 was integrated to tau=220000 and approached the low
  equilibrium, supporting its interpretation as a slowly damped transient.
- **Table 1:** recomputed the twelve Hopf/cycle-fold widths and the four
  recruitment-path widths from compact periodic seeds, evaluating periodic
  closure, phase and the unit transverse multiplier rather than copying
  endpoint numbers. All widths were positive and fold residual checks passed.

The release's supplied figure/continuation tables were then restored from
the original verified results, retaining only the rows used in the main
figures and table. Tests did not change the original project's data or figures.
Generated figures and raw AUTO work files are excluded from version control.

This record concerns numerical reproduction of the current five main figures
and Table 1. It is not an interval-arithmetic certification of the LLC diagrams.
