"""Figure 4: read the bundled data and render the approved six panels.
Regenerate its data with generate_figure_04.py and auto/continue.auto.
"""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import FixedLocator, FormatStrFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE
OUT = HERE / 'outputs'
DATA = ROOT / 'data'
PERIODIC = DATA / 'periodic'
CONT = DATA / 'continuation'
TRAJ = DATA / 'trajectories'
MECH = ROOT / 'data' / 'equilibria'
PLANE_DIR = DATA / 'planes'
MM = 1 / 25.4
BLUE, ORANGE, GREEN = '#2463A2', '#C87726', '#297254'
GREY, PURPLE, INK, LIGHT = '#7D7D7D', '#74429A', '#30343B', '#DDE1E5'
DASH = (0, (4, 2.5))
LPC_LINE = (0, (.1, 2.2))
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
    'mathtext.fontset': 'stix', 'font.size': 10,
    'axes.titlesize': 10.5, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'axes.linewidth': 0.75, 'axes.edgecolor': INK,
    'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK, 'ytick.color': INK, 'lines.linewidth': 1.35,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'tumour-immune-figure4-dimensionless',
    'savefig.facecolor': 'white', 'savefig.bbox': None,
})

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def load_json(path):
    return json.loads(path.read_text(encoding='utf-8'))

def n(row, key):
    return float(row[key])

def fbool(value):
    return value == 'True' if value in ('True', 'False') else float(value) > 0.5

def style(ax, letter, title):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color=LIGHT, lw=0.5, alpha=0.65)
    ax.tick_params(direction='out', length=3, width=0.65, pad=3)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title(title, loc='left', pad=10, fontsize=10)
    ax.text(-0.16, 1.05, letter, transform=ax.transAxes,
            fontsize=12, weight='bold', va='bottom')

def paths(ax, rows, xkey, ykey, *, mode, xscale=1, xref=0,
          yscale=1, fold_common=False, markers=False):
    """Keep source order and break opposite stability/invalid segments.

    The saved subcritical Hopf--LPC branch is saddle before the unique
    LPC and stable after it, supported by existing independent checks.
    Tiny raw sign fluctuations at the zero-amplitude end are not rendered
    as additional stable windows. Every valid amplitude sample is kept.
    """
    xx = np.array([n(r, xkey) * xscale - xref for r in rows])
    yy = np.array([n(r, ykey) * yscale for r in rows])
    shared = np.array([r.get('type') in ('LP', 'HB') for r in rows])
    if mode in ('overlay', 'metrics'):
        flags = np.array([fbool(r['stable']) for r in rows])
        valid = np.ones(len(rows), dtype=bool)
    else:
        key = 'pt' if mode == 'baseline' else 'signed_pt'
        flags = np.array([int(r[key]) < 0 for r in rows])
        valid = np.array([r.get('plot_valid', r.get('valid', 'True')) == 'True'
                          and r.get('type') != 'MX' for r in rows])
        if mode == 'periodic':
            fold = next(i for i, r in enumerate(rows) if r['type'] == 'LP')
            flags = np.arange(len(rows)) > fold
    cycle = mode in ('periodic', 'metrics') or rows[0].get('kind') == 'cycle'
    for stable in (True, False):
        color = (BLUE if stable else ORANGE) if cycle else (GREEN if stable else GREY)
        mask = flags == stable
        if fold_common:
            mask |= shared
        ax.plot(xx, np.where(mask & valid, yy, np.nan), color=color,
                ls='-' if stable else DASH, lw=1.35,
                marker='o' if markers else None, ms=2.7,
                mfc='white', mew=0.75, zorder=3)

def local_guides(ax, fold_alpha, hopf_alpha):
    ax.axvspan(fold_alpha, hopf_alpha, color=PURPLE, alpha=0.065, zorder=0)
    ax.axvline(fold_alpha, color=PURPLE, lw=0.7, ls=(0, (2, 3)), zorder=1)
    ax.axvline(hopf_alpha, color=GREY, lw=0.7, ls=(0, (2, 3)), zorder=1)
    ax.set_xlim(hopf_alpha-1.85e-4, hopf_alpha+3.15e-4)
    ax.xaxis.set_major_locator(FixedLocator([.4987,.4988,.4989,.4990,.4991]))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    ax.set_xlabel(r'T-cell expansion $\alpha$')
    ax.xaxis.label.set_fontsize(9.5)

def inset_style(ax):
    ax.set_facecolor('white')
    ax.tick_params(labelsize=8, direction='out', length=2, width=0.55, pad=1.5)
    for spine in ax.spines.values():
        spine.set_color(GREY)
        spine.set_linewidth(0.6)

def panel_legend(ax, handles, *, loc='upper left', anchor=None, columns=1, size=8.2):
    """Compact legend inside its own panel, placed in unoccupied space."""
    return ax.legend(handles=handles, loc=loc, bbox_to_anchor=anchor,
                     ncol=columns, fontsize=size, handlelength=1.65,
                     handletextpad=.5, columnspacing=.8, borderaxespad=.25,
                     labelspacing=.30, frameon=True, facecolor='white',
                     edgecolor='none', framealpha=.96, borderpad=.3)





def mechanism_plane(ax, rows, plane):
    """Actual dimensionless loci plus a resolved Hopf--LPC inset."""
    xkey='sigma' if plane=='sigma_alpha' else 'beta'
    for family,color in [('SN',GREEN),('Hopf',INK),('LPC',PURPLE)]:
        branches=dict.fromkeys(q['branch'] for q in rows if q['plane']==plane and q['family']==family)
        for branch in branches:
            part=[q for q in rows if q['plane']==plane and q['branch']==branch]
            valid=np.array([q['valid']=='True' for q in part])
            xx=np.array([n(q,xkey) for q in part]);yy=np.array([n(q,'alpha') for q in part])
            ax.plot(np.where(valid,xx,np.nan),np.where(valid,yy,np.nan),color=color,
                    ls=LPC_LINE if family=='LPC' else '-',dash_capstyle='round',
                    lw=1.5 if family in ('Hopf','LPC') else 1.3,zorder=3)
    inset=ax.inset_axes((.54,.625,.44,.315))
    # Both insets resolve the saved Hopf--LPC gap in actual alpha.
    # Paired Hopf references are already stored at the LPC samples.
    for branch in dict.fromkeys(q['branch'] for q in rows if q['plane']==plane and q['family']=='LPC'):
        part=[q for q in rows if q['plane']==plane and q['branch']==branch]
        valid=np.array([q['valid']=='True' and bool(q['hopf_reference']) for q in part])
        xx=np.array([n(q,xkey) for q in part])
        yy=np.array([n(q,'alpha') for q in part])
        hh=np.array([n(q,'hopf_reference') if q['hopf_reference'] else np.nan for q in part])
        xx,yy,hh=[np.where(valid,v,np.nan) for v in (xx,yy,hh)]
        inset.plot(xx,hh,color=INK,lw=1.05,zorder=3)
        inset.plot(xx,yy,color=PURPLE,ls=LPC_LINE,lw=1.35,dash_capstyle='round',zorder=4)
        inset.fill_between(xx,yy,hh,color=PURPLE,alpha=.10,zorder=1)
    inset.set_title('Hopf-LPC detail',fontsize=7.2,pad=3)
    inset.set_ylabel(r'$\alpha$',fontsize=7.2,labelpad=0)
    if plane=='sigma_alpha':
        inset.set(xlim=(0,.02),ylim=(.49850,.49915),
                  xticks=[0,.01,.02],yticks=[.4986,.4988,.4990])
        ax.set(xlim=(-.004,.16),ylim=(.452,.56),xlabel=r'Recruitment $\sigma$',ylabel=r'T-cell expansion $\alpha$')
        ax.set_xticks([0,.04,.08,.12,.16]);ax.set_yticks([.46,.49,.52,.55])
        style(ax,'d','Recruitment-expansion plane')
    else:
        inset.set(xlim=(.1573,.15755),ylim=(.4985,.4990),xticks=[.1573,.1574,.1575],yticks=[.4986,.4988,.4990])
        ax.set(xlim=(.10,.24),ylim=(.42,.60),xlabel=r'DC activation $\beta$',ylabel=r'T-cell expansion $\alpha$')
        ax.set_xticks([.10,.14,.18,.22]);ax.set_yticks([.44,.48,.52,.56,.60])
        style(ax,'e','Expansion-DC activation plane')
    inset_style(inset)
    inset.tick_params(labelsize=7.,length=2,pad=1)
    inset.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    inset.xaxis.set_major_formatter(FormatStrFormatter('%.2f' if plane=='sigma_alpha' else '%.4f'))
    inset.set_xlabel(r'$\sigma$' if plane=='sigma_alpha' else r'$\beta$',
                     fontsize=7.2,labelpad=1)
    inset.spines[['top','right']].set_visible(False)
    # Give d's source locator enough display height to avoid hiding the
    # narrow band. The inset ticks retain the exact zoom coordinates.
    xlo,xhi=inset.get_xlim();ylo,yhi=inset.get_ylim()
    source_ylo,source_yhi=ylo,yhi
    if plane=='sigma_alpha':
        centre=.5*(ylo+yhi)
        source_ylo,source_yhi=centre-.0035,centre+.0035
    ax.add_patch(Rectangle((xlo,source_ylo),xhi-xlo,source_yhi-source_ylo,
                           fill=False,edgecolor=GREY,lw=.8,zorder=2))
    ax.add_artist(ConnectionPatch(xyA=(xhi,source_yhi),coordsA=ax.transData,
                                  xyB=(0,0) if plane=='sigma_alpha' else (1,0),
                                  coordsB=inset.transAxes,
                                  arrowstyle='-|>',mutation_scale=8,
                                  shrinkA=2,shrinkB=3,color=GREY,lw=.75,
                                  connectionstyle='arc3' if plane=='sigma_alpha'
                                  else 'angle,angleA=0,angleB=90,rad=3',
                                  clip_on=False,
                                  zorder=4.5 if plane=='sigma_alpha' else 6))
    # Keep guides and the callout out of the small tick/axis labels.
    for label in [*inset.get_xticklabels(),*inset.get_yticklabels(),
                  inset.xaxis.label,inset.yaxis.label]:
        label.set_bbox(dict(facecolor='white',edgecolor='none',pad=.1))
    handles=[Line2D([],[],color=GREEN,label='SN'),Line2D([],[],color=INK,label='Hopf'),
             Line2D([],[],color=PURPLE,ls=LPC_LINE,lw=1.5,dash_capstyle='round',label='LPC')]
    panel_legend(ax,handles,loc='upper left' if plane=='sigma_alpha' else 'lower left',
                 anchor=(.01,.985) if plane=='sigma_alpha' else (.015,.025),
                 columns=1,size=7.7)


def main():
    config = load_json(DATA / 'parameters.json')['analysis']
    r = config['dimensional']['r']
    r0 = config['r0_per_day']
    rho = r / r0
    K = config['dimensional']['K']
    baseline = read_csv(CONT / 'branches.csv')
    baseline_landmarks = read_csv(CONT / 'landmarks.csv')
    raw = read_csv(PERIODIC / 'branches.csv')
    metrics = read_csv(PERIODIC / 'branch_metrics.csv')
    overlay = read_csv(PERIODIC / 'equilibrium_overlay.csv')
    landmarks = read_csv(PERIODIC / 'landmarks.csv')
    mechanism_rows = read_csv(PLANE_DIR/'comparison.csv')
    fold = load_json(PERIODIC / 'fold_verified.json')
    hopf = load_json(PERIODIC / 'hopf_reference.json')
    alpha_h = hopf['alpha0'] / rho
    hopf_alpha = alpha_h
    fold_alpha = fold['alpha0']/rho
    stable_cycle = np.genfromtxt(PERIODIC / 'stable_0p96925_corrected.csv', delimiter=',', names=True)
    saddle_cycle = np.genfromtxt(PERIODIC / 'saddle_0p96925_corrected.csv', delimiter=',', names=True)
    stable_info = load_json(PERIODIC / 'stable_0p96925_corrected.json')
    saddle_info = load_json(PERIODIC / 'saddle_0p96925_corrected.json')
    equilibria = [q for q in read_csv(MECH / 'common_parameter_equilibria.csv')
                  if q['model'] == 'M4' and q['kind'] == 'internal']
    po = [q for q in raw if q['branch'] == 'po_fw' and q['kind'] == 'cycle']
    orbits = [(stable_cycle, BLUE, '-'), (saddle_cycle, ORANGE, DASH)]

    fig = plt.figure(figsize=(180 * MM, 185 * MM))
    a = fig.add_axes((0.105, 0.755, 0.365, 0.205))
    b = fig.add_axes((0.605, 0.755, 0.365, 0.205))
    c = fig.add_axes((0.105, 0.445, 0.365, 0.205))
    d = fig.add_axes((0.605, 0.445, 0.365, 0.205))
    e = fig.add_axes((0.105, 0.115, 0.365, 0.205))
    f = fig.add_axes((0.605, 0.115, 0.365, 0.205))

    eq_stable = Line2D([],[],color=GREEN,label='Stable eq.')
    eq_unstable = Line2D([],[],color=GREY,ls=DASH,label='Unstable eq.')
    cycle_stable = Line2D([],[],color=BLUE,label='Stable cycle')
    cycle_saddle = Line2D([],[],color=ORANGE,ls=DASH,label='Saddle cycle')
    lpc_handle = Line2D([],[],color=PURPLE,label='LPC')
    hopf_handle = Line2D([],[],color=INK,label='Hopf')

    for branch in ('eq_alpha_low_bw', 'eq_alpha_high_fw', 'eq_alpha_high_bw', 'po_alpha_h1_phys_fw'):
        rows = [q for q in baseline if q['branch'] == branch and q['parameter'] == '2'
                and q['second_parameter'] == '0' and q['chart'] == 'physical']
        for key in (('x_min', 'x_max') if rows[0]['kind'] == 'cycle' else ('x_max',)):
            paths(a, rows, 'alpha', key, mode='baseline', xscale=1/rho, fold_common=True)
    for alpha0 in (3.75, 3.9, 4.0, 4.2):
        a.axvline(alpha0 / rho, color=LIGHT, lw=0.65, ls=(0, (2, 3)), zorder=1)
    for group, alpha0 in (('A1', 3.75), ('A2', 4.0)):
        a.axvline(alpha0 / rho, color=GREY, lw=.85, ls=(0, (2, 3)), zorder=2)
        a.text(alpha0 / rho, .035, group, fontsize=8.5, ha='center', color=INK)
    for filename in ('cycle_alpha_4_metrics.csv', 'cycle_alpha_4p2_metrics.csv'):
        q = read_csv(TRAJ / filename)[0]
        for key in ('C_min_mm3', 'C_max_mm3'):
            a.plot(n(q, 'alpha0') / rho, n(q, key) / K, 'o', ms=3.6,
                   color=BLUE, mfc='white', mew=0.8, zorder=5)
    sn = next(q for q in baseline_landmarks if q['branch'] == 'eq_alpha_low_bw'
              and q['type'] == 'LP' and n(q, 'alpha') < 4)
    hb = next(q for q in baseline_landmarks if q['branch'] == 'eq_alpha_low_bw' and q['type'] == 'HB')
    a.plot(n(sn, 'alpha') / rho, n(sn, 'x_max'), 'D', ms=4.2,
           mfc='white', mec=GREEN, mew=1, zorder=5)
    a.plot(n(hb, 'alpha') / rho, n(hb, 'x_max'), 'o', ms=4.2,
           mfc='white', mec=INK, mew=1, zorder=5)
    a.annotate('SN', (n(sn, 'alpha') / rho, n(sn, 'x_max')),
               xytext=(1.80, 0.66), fontsize=8.8,
               arrowprops=dict(arrowstyle='-', lw=0.6, color=GREY))
    a.annotate('Hopf', (n(hb, 'alpha') / rho, n(hb, 'x_max')),
               xytext=(1.965, 0.57), fontsize=8.8,
               arrowprops=dict(arrowstyle='-', lw=0.6, color=GREY))
    a.text(.97,.86,r'$\gamma_f$',transform=a.transAxes,
           fontsize=8.5,ha='right')
    a.set(xlim=(1.775, 2.215), ylim=(0, 1.065),
          xlabel=r'T-cell expansion $\alpha$', ylabel=r'Tumour fraction $x$')
    a.xaxis.set_major_locator(FixedLocator([1.8, 1.9, 2.0, 2.1, 2.2]))
    a.yaxis.set_major_locator(FixedLocator([0, 0.5, 1]))
    style(a, 'a', 'SN-Hopf structure')
    panel_legend(a,[eq_stable,eq_unstable,cycle_stable],anchor=(.01,.91),columns=2,size=8.)

    for lo, hi in ((0, 0.5), (0.5, 0.95), (0.95, 1.01)):
        rows = [q for q in overlay if lo < n(q, 'x') < hi]
        paths(b, rows, 'lambdaT_per_day', 'x', mode='overlay',
              xscale=1/r)
    for key in ('x_min', 'x_max'):
        paths(b, po, 'alpha0', key, mode='periodic',
              xscale=1/rho, fold_common=True)
    fold_row = next(q for q in po if q['type'] == 'LP')
    for key in ('x_min', 'x_max'):
        b.plot(fold_alpha, n(fold_row, key), 'D', color=PURPLE, ms=3.8, zorder=5)
    hb_local = next(q for q in landmarks if q['kind'] == 'equilibrium' and q['type'] == 'HB')
    b.plot(hopf_alpha, n(hb_local, 'x_min'), 'o', color=INK, mfc='white', ms=4.2, zorder=5)
    b.text(0.03, 0.84, r'$\gamma_s$', transform=b.transAxes, fontsize=8.7)
    high = max(equilibria, key=lambda q: n(q, 'x'))
    low = min(equilibria, key=lambda q: n(q, 'x'))
    middle = next(q for q in equilibria if not fbool(q['stable']))
    b.annotate(r'$E_H$', xy=(hopf_alpha+2.3e-4, n(high, 'x')), xytext=(hopf_alpha+2.1e-4, 0.88),
               fontsize=9, arrowprops=dict(arrowstyle='-', lw=0.6, color=GREEN))
    b.text(fold_alpha, 0.48, 'LPC', fontsize=8.5, color=PURPLE, ha='center')
    b.annotate('Hopf',(hopf_alpha,n(hb_local,'x_min')),xytext=(hopf_alpha+.55e-4,.34),
               fontsize=8.5,color=INK,
               arrowprops=dict(arrowstyle='-',lw=.55,color=GREY))
    b.set(ylim=(0, 1.065), ylabel=r'Tumour fraction $x$')
    b.yaxis.set_major_locator(FixedLocator([0, 0.5, 1]))
    local_guides(b, fold_alpha, hopf_alpha)
    for group, offset in (('B1', -1.8), ('B2', 3.0)):
        b.axvline(hopf_alpha+offset*1e-4, color=GREY, lw=.85, ls=(0, (2, 3)), zorder=2)
        b.text(hopf_alpha+offset*1e-4, .035, group, fontsize=8.5,
               ha='left' if group == 'B1' else 'right', color=INK)
    style(b, 'b', 'Hopf-LPC coexistence')
    panel_legend(b,[eq_stable,eq_unstable,cycle_stable,cycle_saddle],
                 loc='upper right',anchor=(.99,.70),columns=2,size=8.)

    for array, color, line in orbits:
        c.plot(array['x'], array['y'], color=color, ls=line, lw=1.3)
    c.plot(n(low, 'x'), n(low, 'y'), 'o', color=GREEN, ms=4.7, zorder=5)
    c.set(xlim=(0.08, 0.72), ylim=(0.54, 0.95),
          xlabel=r'Tumour fraction $x$', ylabel=r'Effector state $y$')
    c.xaxis.set_major_locator(FixedLocator([0.1, 0.3, 0.5, 0.7]))
    c.yaxis.set_major_locator(FixedLocator([0.6, 0.75, 0.9]))
    style(c, 'c', 'Cycle-equilibrium coexistence')
    inset = c.inset_axes((0.50, 0.32, 0.49, 0.63))
    for array, color, line in orbits:
        inset.plot(array['x'], array['y'], color=color, ls=line, lw=0.95)
    for eq in equilibria:
        inset.plot(n(eq, 'x'), n(eq, 'y'), 'o' if fbool(eq['stable']) else 'x',
                   color=GREEN if fbool(eq['stable']) else GREY,
                   ms=4.0, mew=1.0, zorder=5)
    inset.annotate(r'$E_L$', (n(low, 'x'), n(low, 'y')),
                   xytext=(0.04, 0.93), fontsize=8.4, color=GREEN,
                   arrowprops=dict(arrowstyle='-', lw=0.5, color=GREEN))
    inset.annotate(r'$E_S$', (n(middle, 'x'), n(middle, 'y')),
                   xytext=(0.51, 0.40), fontsize=8.4, color=GREY,
                   arrowprops=dict(arrowstyle='-', lw=0.5, color=GREY))
    inset.annotate(r'$E_H$', (n(high, 'x'), n(high, 'y')),
                   xytext=(0.53, 0.075), fontsize=8.4, color=GREEN,
                   arrowprops=dict(arrowstyle='-', lw=0.5, color=GREEN))
    inset.set(xlim=(-0.035, 1.075), ylim=(-0.04, 1.06), xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    inset.set_title(r'Full $x$-$y$ projection', fontsize=8.4, pad=3)
    inset_style(inset)
    panel_legend(c,[Line2D([],[],color=BLUE,label=r'$P_s$'),
                    Line2D([],[],color=ORANGE,ls=DASH,label=r'$P_u$'),
                    Line2D([],[],color=GREEN,marker='o',ls='none',ms=3.5,label=r'$E_L$')],
                 loc='lower left',anchor=(.015,.04),size=8.5)

    normalized=config['normalized']
    fixed=dict(sigma=normalized['sT']/rho**2,beta=normalized['lambdaD']/rho,
               delta=normalized['dT']/rho,a=normalized['KA']/rho,
               gamma=.05/rho,kappa=normalized['KC'],n=4)
    mechanism_plane(d,mechanism_rows,'sigma_alpha')
    mechanism_plane(e,mechanism_rows,'alpha_beta')

    paths(f, metrics, 'alpha0', 'leading_multiplier', mode='metrics',
          xscale=1/rho, fold_common=True, markers=True)
    f.axhline(1, color=GREY, lw=0.85, ls=(0, (2, 3)), zorder=2)
    f.plot(fold_alpha, fold['leading_transverse_multiplier'], 'D', color=PURPLE, ms=4.5, zorder=5)
    f.set(ylim=(0.87, 1.009), ylabel=r'Floquet multiplier $\mu_\perp$')
    f.yaxis.set_major_locator(FixedLocator([0.9, 0.95, 1]))
    f.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    local_guides(f, fold_alpha, hopf_alpha)
    style(f, 'f', 'Floquet stability exchange')
    zoom = f.inset_axes((0.59, 0.65, 0.40, 0.34))
    paths(zoom, metrics, 'alpha0', 'leading_multiplier', mode='metrics',
          xscale=1/rho, fold_common=True, markers=True)
    zoom.axvspan(fold_alpha, hopf_alpha, color=PURPLE, alpha=0.065, zorder=0)
    zoom.axhline(1, color=GREY, lw=0.7, ls=(0, (2, 3)), zorder=1)
    zoom.plot(fold_alpha, fold['leading_transverse_multiplier'], 'D', color=PURPLE, ms=3.8, zorder=5)
    zoom.set(xlim=(hopf_alpha-1.70e-4, hopf_alpha+0.10e-4), ylim=(0.998, 1.010), xticks=[.4987,.4988], yticks=[1, 1.004, 1.008])
    zoom.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    zoom.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    zoom.text(0.50, 0.96, r'Near $\mu_\perp=1$', transform=zoom.transAxes,
              ha='center', va='top', fontsize=8.2)
    inset_style(zoom)
    # Keep parent reference guides out of the inset's exterior tick text.
    for label in zoom.get_yticklabels():
        label.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.15))

    panel_legend(f,[cycle_stable,cycle_saddle],loc='lower left',anchor=(.015,.035),size=8.2)
    OUT.mkdir(exist_ok=True)
    for extension in ('pdf', 'svg', 'png'):
        destination=OUT / f'Figure_04_dimensionless.{extension}'
        temporary=destination.with_name(f'.{destination.stem}.tmp.{extension}')
        fig.savefig(temporary, dpi=300)
        temporary.replace(destination)
    plt.close(fig)
    print('Saved retained Figure 4 in manuscript dimensionless coordinates.')

if __name__ == '__main__':
    main()
