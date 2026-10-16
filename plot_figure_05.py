"""Figure 5: plot the four approved initial-state comparisons from data/."""
from pathlib import Path
import csv
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator

HERE = Path(__file__).resolve().parent
ROOT = HERE
OUT = HERE / 'outputs'
DATA = ROOT / 'data'
PARAMETERS = DATA / 'parameters.json'
SERIES = DATA / 'timeseries' / 'Figure_05_timeseries.csv'
META = DATA / 'timeseries' / 'Figure_05_parameters.json'
GREEN, BLUE, INK = '#297254', '#2463A2', '#30343B'
ORANGE = '#C87726'
DASH = (0,(4,2.5))
GREY, LIGHT = '#7D7D7D', '#DDE1E5'
MM = 1 / 25.4
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
    'mathtext.fontset': 'stix', 'font.size': 10,
    'axes.titlesize': 10, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'axes.linewidth': .75, 'axes.edgecolor': INK,
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK, 'ytick.color': INK, 'lines.linewidth': 1.35,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'tumour-immune-figure5-timeseries',
    'savefig.facecolor': 'white', 'savefig.bbox': None,
})












def axes_style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color=LIGHT, lw=.5, alpha=.65)
    ax.tick_params(direction='out', length=3, width=.65, pad=3)
    ax.spines[['top','right']].set_visible(False)



def plot():
    OUT.mkdir(exist_ok=True)
    if not SERIES.exists() or not META.exists():
        raise FileNotFoundError('First run this script with --simulate to create the four selected slices.')
    meta = json.loads(META.read_text(encoding='utf-8'))
    with SERIES.open(encoding='utf-8', newline='') as stream:
        data = list(csv.DictReader(stream))
    fig = plt.figure(figsize=(180*MM, 160*MM))
    # Four initial states per panel. The tumour coordinate resolves basin
    # selection without superposing twelve x/y/z curves; all states remain
    # available in the CSV. Draw the early transient from the true time origin.
    windows = {'A1':(0.,450.), 'A2':(0.,450.),
               'B1':(0.,1800.), 'B2':(0.,1800.)}
    for i,g in enumerate(meta['groups']):
        left = .105 if i%2==0 else .605
        bottom = .635 if i<2 else .150
        ax = fig.add_axes((left,bottom,.365,.295))
        axes_style(ax)
        start,end = windows[g['group']]
        handles = []
        for name in ('IC1','IC2','IC3','IC4'):
            group_rows = [r for r in data if r['group']==g['group'] and r['trajectory']==name]
            tt = np.array([float(r['tau']) for r in group_rows])
            keep = (tt>=start)&(tt<=end)
            line = '-' if name in ('IC1','IC3') else DASH
            outcome = g['trajectories'][name]['outcome']
            color = ORANGE if outcome=='E_H' else GREEN if outcome=='E_L' else BLUE
            target = r'$E_H$' if outcome=='E_H' else r'$E_L$' if outcome=='E_L' else r'$P_s$'
            uu = np.array([float(r['x']) for r in group_rows])
            ax.plot(tt[keep]-start,uu[keep],color=color,ls=line,lw=1.15)
            handles.append(Line2D([],[],color=color,ls=line,label=name+': '+target))
        ax.set_xlim(0.,end-start)
        ax.xaxis.set_major_locator(FixedLocator([0.,(end-start)/2,end-start]))
        ax.set_ylim(-.035,1.24)
        ax.set_yticks([0.,.5,1.])
        ax.set_ylabel(r'$x(\tau)$')
        ax.set_xlabel(r'Time $\tau$' if start==0. else r'Time $\tau-'+f'{start:g}'+'$')
        title = ('Two stable equilibria' if g['group'] in ('A1','B1')
                 else 'Stable equilibrium + cycle')
        ax.text(-.16,1.04,g['panel'],transform=ax.transAxes,
                fontsize=12,weight='bold',va='bottom')
        ax.text(0.,1.04,title+' ('+g['group']+')',transform=ax.transAxes,
                fontsize=10,va='bottom')
        ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(.01,1.),
                  ncol=2,fontsize=8.2,handlelength=1.8,columnspacing=1.2,
                  handletextpad=.4,frameon=False,borderaxespad=.15,
                  labelspacing=.3)
        if g['group']=='B1':
            inset = ax.inset_axes((.47,.56,.51,.19))
            for name in ('IC1','IC2'):
                envelope = g['trajectories'][name]['relaxation_envelope']
                for key in ('maxima','minima'):
                    arr = np.asarray(envelope[key])
                    inset.plot(arr[:,0]/1e4,arr[:,1],color=GREEN,
                               ls='-' if name=='IC1' else DASH,lw=.65)
            inset.axhline(g['equilibria'][0]['state'][0],color=GREY,lw=.6,ls=':')
            inset.set(xlim=(0,g['horizon']/1e4),ylim=(.19,.32),xticks=[0,11,22],yticks=[.2,.3])
            inset.set_title('Oscillation envelope',fontsize=7.3,pad=2)
            inset.set_xlabel(r'$\tau/10^4$',fontsize=7.3,labelpad=0)
            inset.tick_params(labelsize=7.,direction='out',length=2,pad=1)
            inset.spines[['top','right']].set_visible(False)
            inset.set_facecolor('white')
    for ext in ('pdf','svg','png'):
        fig.savefig(OUT/f'Figure_05_timeseries.{ext}',dpi=300)
    plt.close(fig)
    print('Saved Figure 5 PDF, SVG and PNG.',flush=True)


if __name__=='__main__':
    plot()
