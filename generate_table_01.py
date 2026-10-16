"""Recompute the twelve Table 1 Hopf/fold widths from saved compact fold seeds.
The shooting residual, phase condition and transverse multiplier are recomputed.
"""
import csv,json
from model import ROOT,R0,BASE,variant,hopfs
from generate_figure_04 import fold,table

def main():
    rows=[]
    for seed in json.loads((ROOT/'data/table_seeds.json').read_text()):
        p=variant(seed['model'],dD=seed['decay']/R0)
        f=fold(seed,p);h=hopfs(p)[0][1]
        row=dict(model=seed['model'],dD_per_day=seed['decay'],lambdaT_H_per_day=h*R0,
                 lambdaT_LPC_per_day=f['alpha0']*R0,width_per_day=(h-f['alpha0'])*R0)
        if row['width_per_day']<=0 or abs(f['leading_transverse_multiplier']-1)>2e-6:
            raise RuntimeError(('Invalid cycle-fold width',row))
        rows.append(row);print(seed['model'],seed['decay'],'width',row['width_per_day'],flush=True)
    target=ROOT/'data/sn2';table(target/'decay_interaction_summary.csv',rows)
    effects=[]
    for decay in sorted(set(r['dD_per_day'] for r in rows)):
        widths={r['model']:r['width_per_day'] for r in rows if r['dD_per_day']==decay}
        effects.append(dict(dD_per_day=decay,no_recruitment_saturation_ratio=widths['M2']/widths['M1'],
                            recruitment_saturation_ratio=widths['M4']/widths['M3']))
    table(target/'decay_interaction_effects.csv',effects)
    recruitment()
    from build_table_01 import main as build
    build()

def recruitment():
    pathway=[]
    for seed in json.loads((ROOT/'data/recruitment_seeds.json').read_text()):
        p=dict(BASE,sT=seed['sigma0']);f=fold(seed,p);h=hopfs(p)[0][1]
        pathway.append(dict(tag=seed['tag'],axis='sigma0',sigma0=p['sT'],
            dD_per_day=p['dD']*R0,lambdaT_H_per_day=h*R0,
            lambdaT_LPC_per_day=f['alpha0']*R0,width_per_day=(h-f['alpha0'])*R0))
    table(ROOT/'data/equilibria/verified_windows.csv',pathway)

if __name__=='__main__':main()
