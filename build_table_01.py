"""Export manuscript Table 1 from saved results using algebraic scaling only.

Run: python build_table_01.py
Outputs: outputs/Table_01_mechanism_modulation.md and the companion CSV.
The original numerical CSV/JSON files are read-only. No solver is imported.
"""
from pathlib import Path
import csv
import json


HERE = Path(__file__).resolve().parent
ROOT = HERE
OUT = HERE / 'outputs'
MECHANISM = ROOT / 'data/equilibria'
FOLLOWUP = ROOT / 'data/sn2'
MODELS = ('M1', 'M2', 'M3', 'M4')
FIELDS = ['gamma', *(f'1e4_W_alpha_{model}' for model in MODELS), 'I']


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    config = json.loads((ROOT / 'data/parameters.json')
                        .read_text(encoding='utf-8'))['analysis']
    growth = config['dimensional']['r']
    rho = growth / config['r0_per_day']
    normalized = config['normalized']
    sigma_ref = normalized['sT'] / rho**2
    summary = read_csv(FOLLOWUP / 'decay_interaction_summary.csv')
    effects = read_csv(FOLLOWUP / 'decay_interaction_effects.csv')
    by_decay = {}
    for sample in summary:
        by_decay.setdefault(float(sample['dD_per_day']), {})[sample['model']] = sample
    effects_by_decay = {float(row['dD_per_day']): row for row in effects}
    rows = []
    for decay, samples in sorted(by_decay.items()):
        effect = effects_by_decay[decay]
        row = {'gamma': decay / growth}
        row.update({f'1e4_W_alpha_{model}':
                    1e4 * float(samples[model]['width_per_day']) / growth
                    for model in MODELS})
        row['I'] = (float(effect['recruitment_saturation_ratio']) /
                    float(effect['no_recruitment_saturation_ratio']))
        rows.append(row)

    OUT.mkdir(exist_ok=True)
    with (OUT / 'Table_01_mechanism_modulation.csv').open(
            'w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    table = ['| DC decay $\\gamma$ | M1 | M2 | M3 | M4 | Interaction $I$ |',
             '|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        cells = [f"{row['gamma']:.4f}"]
        cells.extend(f"{row[f'1e4_W_alpha_{model}']:.4g}" for model in MODELS)
        cells.append(f"{row['I']:.4g}")
        table.append('| ' + ' | '.join(cells) + ' |')

    pathway = [q for q in read_csv(MECHANISM / 'verified_windows.csv')
               if q['axis'] == 'sigma0']
    start = min(pathway, key=lambda q: float(q['sigma0']))
    end = max(pathway, key=lambda q: float(q['sigma0']))
    w_start = float(start['width_per_day']) / growth
    w_end = float(end['width_per_day']) / growth
    recruitment_limit = float(end['sigma0']) / normalized['sT']
    gamma_ref = float(start['dD_per_day']) / growth

    caption = rf'''**Table 1. Mechanism-dependent modulation of the tristable-window width.**

The M1-M4 entries report $10^4W_\alpha$, where
$W_\alpha=\alpha_H-\alpha_{{\mathrm{{LPC}}}}$ is the actual dimensionless
expansion-parameter width between Hopf and the limit point of cycles (LPC).
The factor $10^4$ is a display multiplier. M1 and M2 have no recruitment
and use linear and saturated DC activation, respectively; M3 and M4 use
the corresponding activation functions with recruitment
$\sigma=\sigma_\ast$, where $\sigma_\ast={sigma_ref:.8f}$.
The linear and saturated responses have the same low-signal slope.
The interaction contrast is
$I=(W_{{M4}}/W_{{M3}})/(W_{{M2}}/W_{{M1}})$, with all widths evaluated
at the same $\gamma$. Values above 1 indicate a larger recruitment
retention factor $W_{{M4}}/W_{{M2}}$ for saturated DC activation than
$W_{{M3}}/W_{{M1}}$ for linear activation.
Each of the twelve saved slices has numerical evidence of two stable
equilibria, a stable periodic orbit, and an accompanying saddle periodic
orbit. All parameters and tabulated quantities are dimensionless.
The displayed $\gamma$ values are rounded to four decimal places and
widths and $I$ to four significant digits; the companion CSV retains
the unrounded converted values.

Fixed parameters are $\delta={normalized['dT']/rho:.8f}$,
$\beta={normalized['lambdaD']/rho:.8f}$,
$a={normalized['KA']/rho:.8f}$,
$\kappa={normalized['KC']:.2g}$ and $n={normalized['n']}$.
M1/M2 use $\sigma=0$ and M3/M4 use $\sigma=\sigma_\ast$.
'''
    results = rf'''Across the three DC-decay slices, slower decay broadens the
Hopf-LPC coexistence windows in all four models (Table 1). At each slice,
DC saturation increases the window width, while recruitment reduces it.
The interaction contrast approaches 1 at the smaller decay parameters,
showing that the relative recruitment effect becomes similar for linear
and saturated activation within these tested slices.
'''
    recruitment = rf'''Along the recruitment pathway in Figure 4d, at
$\gamma={gamma_ref:.8f}$, increasing $\sigma/\sigma_\ast$ from 0 to
{recruitment_limit:g} reduces $W_\alpha$ from
${w_start/1e-4:.3g}\times10^{{-4}}$ to
${w_end/1e-6:.3g}\times10^{{-6}}$, retaining
{100*w_end/w_start:.2f}% of the zero-recruitment width.
'''
    provenance = 'Source data: data/sn2/decay_interaction_summary.csv and decay_interaction_effects.csv; data/equilibria/verified_windows.csv; data/parameters.json.\n'
    document = ('# Table 1\n\n' + '\n'.join(table) + '\n\n' + caption +
                '\n## Suggested Results text\n\n' + results +
                '\n## Recruitment-pathway sentence for Results\n\n' + recruitment +
                '\n' + provenance)
    (OUT / 'Table_01_mechanism_modulation.md').write_text(document, encoding='utf-8')
    print('Exported Table 1 Markdown and CSV from existing results; no simulation.')


if __name__ == '__main__':
    main()
