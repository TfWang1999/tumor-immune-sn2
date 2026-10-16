"""Regenerate the four Figure 5 slices using the original DOP853 routine."""
from pathlib import Path
import csv,json
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs'
DATA=ROOT/'data'
PARAMETERS=DATA/'parameters.json'
SERIES=DATA/'timeseries/Figure_05_timeseries.csv'
META=DATA/'timeseries/Figure_05_parameters.json'

def configuration():
    c = json.loads(PARAMETERS.read_text(encoding='utf-8'))['analysis']
    rho = c['normalized']['r']
    q = c['normalized']
    p = dict(sigma=q['sT']/rho**2, delta=q['dT']/rho,
             beta=q['lambdaD']/rho, a=q['KA']/rho,
             kappa=q['KC'], n=q['n'])
    h = json.loads((DATA/'periodic/hopf_reference.json').read_text())['alpha0']/rho
    selections = [
        dict(group='A1', panel='a', figure4_panel='a', alpha=3.75/rho,
             gamma=.5/rho, horizon=900., head_end=300., tail_length=300.,
             regime='two stable equilibria', title='Two equilibria: A1'),
        dict(group='A2', panel='b', figure4_panel='a', alpha=4./rho,
             gamma=.5/rho, horizon=1500., head_end=450., tail_length=300.,
             regime='stable equilibrium and stable cycle', title='Equilibrium + cycle: A2'),
        dict(group='B1', panel='c', figure4_panel='b', alpha=h-1.8e-4,
             gamma=.05/rho, horizon=120000., head_end=600., tail_length=600.,
             regime='two stable equilibria', title='Two equilibria: B1'),
        dict(group='B2', panel='d', figure4_panel='b', alpha=h+3.e-4,
             gamma=.05/rho, horizon=120000., head_end=600., tail_length=600.,
             regime='stable equilibrium and stable cycle', title='Equilibrium + cycle: B2'),
    ]
    return rho, p, h, selections

def rhs(tau, u, p):
    """Manuscript growth-time model: tau=r*t, x=C/K, y=E/r, z=Z."""
    x, y, z = u
    h = x**p['n'] / (p['kappa']**p['n'] + x**p['n'])
    q = x*y
    return np.array([x*(1-x-y),
                     p['sigma']*h + (p['alpha']*z/(1+z)-p['delta'])*y,
                     p['beta']*q/(p['a']+q)-p['gamma']*z])

def jacobian(u, p):
    x, y, z = u
    k, n = p['kappa']**p['n'], p['n']
    q = x*y
    b = p['beta']*p['a']/(p['a']+q)**2
    return np.array([[1-2*x-y, -x, 0.],
                     [p['sigma']*n*k*x**(n-1)/(k+x**n)**2,
                      p['alpha']*z/(1+z)-p['delta'],
                      p['alpha']*y/(1+z)**2],
                     [b*y, b*x, -p['gamma']]])

def selected_equilibria(p):
    """Locate the three internal equilibria for one selected slice only."""
    from scipy.optimize import brentq
    def state(x):
        y = 1-x
        z = p['beta']/p['gamma'] * x*y/(p['a']+x*y)
        return np.array([x, y, z])
    def balance(x):
        return rhs(0., state(x), p)[1]
    grid = np.r_[np.geomspace(1.e-8, .05, 50), np.linspace(.05, 1-1.e-10, 600)]
    result = []
    for left, right in zip(grid[:-1], grid[1:]):
        if balance(left)*balance(right) < 0:
            u = state(brentq(balance, left, right, xtol=1.e-14))
            ev = np.linalg.eigvals(jacobian(u, p))
            result.append(dict(state=u.tolist(), stable=bool(ev.real.max() < 0),
                               spectral_abscissa=float(ev.real.max()),
                               eigenvalues=[[float(v.real), float(v.imag)] for v in ev]))
    if len(result) != 3:
        raise RuntimeError('Expected the three internal equilibria on this selected slice.')
    return result

def simulate():
    from scipy.integrate import solve_ivp
    from scipy.signal import find_peaks
    rho, fixed, hopf, groups = configuration()
    rows = []
    for g in groups:
        p = dict(fixed, alpha=g['alpha'], gamma=g['gamma'])
        eq = selected_equilibria(p)
        g['equilibria'] = eq
        saddle = np.array(eq[1]['state'])
        ev, vectors = np.linalg.eig(jacobian(saddle, p))
        unstable = np.flatnonzero(ev.real > 0.)
        if len(unstable) != 1 or abs(ev[unstable[0]].imag) > 1.e-10:
            raise RuntimeError('The intermediate saddle must have one real unstable direction.')
        direction = vectors[:, unstable[0]].real
        direction /= np.max(np.abs(direction))
        if direction[0] < 0:
            direction *= -1
        g['saddle_selection'] = dict(state=saddle.tolist(),
            unstable_eigenvalue=float(ev[unstable[0]].real),
            normalized_unstable_eigenvector=direction.tolist(),
            offsets={'IC1':-.03, 'IC2':-.08, 'IC3':.03, 'IC4':.08})
        initial = {name:saddle+eps*direction for name,eps in
                   g['saddle_selection']['offsets'].items()}
        if min(u.min() for u in initial.values()) <= 0:
            raise RuntimeError('Saddle-near initial states must be positive.')
        # The slow stable focus has a very small decay rate: keep the long
        # integration for classification while displaying its early transient.
        g['horizon'] = 220000. if g['group']=='B1' else g['horizon']
        g['head_end'] = 1800. if g['group'].startswith('B') else 450.
        g['trajectories'] = {}
        total = g['horizon']
        tail_start = total-g['tail_length']
        # Save dense, resolved head/tail windows and a resolved final check
        # window. Intermediate long transients are not sampled for the figure.
        head = np.linspace(0., g['head_end'], int(g['head_end']*4)+1)
        tail = np.linspace(tail_start, total, int(g['tail_length']*4)+1)
        check_start = max(g['head_end'], total-5000.)
        check = np.linspace(check_start, total, int((total-check_start)*2)+1)
        times = np.unique(np.r_[head, check, tail])
        for name, u0 in initial.items():
            print(f"Integrating {g['group']} / {name}: tau 0--{total:g}", flush=True)
            sol = solve_ivp(lambda t,u: rhs(t,u,p), (0.,total), u0,
                            method='DOP853', rtol=2.e-10, atol=2.e-12,
                            max_step=6., t_eval=times, dense_output=g['group']=='B1')
            if not sol.success:
                raise RuntimeError(sol.message)
            if sol.y.min() < -1.e-10:
                raise RuntimeError('Selected trajectory left the nonnegative state region.')
            late = sol.t >= tail_start
            yy = sol.y[0, late]
            expected = ('E_H' if name in ('IC3','IC4') else
                        ('E_L' if g['group'] in ('A1','B1') else 'stable_cycle'))
            # Classify from the computed state, not from the perturbation sign.
            distances = {tag:float(np.max(abs(sol.y[:,-1]-np.array(ref['state']))))
                         for tag,ref in [('E_L',eq[0]),('E_H',eq[-1])]}
            nearest = min(distances, key=distances.get)
            outcome = nearest if distances[nearest]<1.e-4 else 'stable_cycle'
            if outcome != expected:
                raise RuntimeError(f'{g["group"]}/{name}: got {outcome}, expected {expected}.')
            diagnostic = dict(initial_state=u0.tolist(), outcome=outcome,
                              final_state=sol.y[:,-1].tolist(),
                              final_rhs_norm=float(np.max(np.abs(rhs(total,sol.y[:,-1],p)))),
                              last_window_x_min=float(yy.min()),
                              last_window_x_max=float(yy.max()), nfev=sol.nfev)
            if outcome != 'stable_cycle':
                ref = np.array(eq[0 if outcome=='E_L' else -1]['state'])
                diagnostic['final_equilibrium_distance'] = float(np.max(abs(sol.y[:,-1]-ref)))
                if diagnostic['final_equilibrium_distance'] > 1.e-4:
                    raise RuntimeError(f'{g["group"]}: selected equilibrium trajectory has not settled.')
            else:
                ii = sol.t >= check_start
                peaks, _ = find_peaks(sol.y[0,ii])
                peak_t = sol.t[ii][peaks]
                if len(peaks) < 10:
                    raise RuntimeError(f'{g["group"]}: too few late peaks for the cycle summary.')
                diagnostic['late_period_mean'] = float(np.diff(peak_t[-10:]).mean())
                diagnostic['last_peak_heights'] = sol.y[0,ii][peaks[-6:]].tolist()
                if np.ptp(diagnostic['last_peak_heights']) > 2.e-4:
                    raise RuntimeError(f'{g["group"]}: selected oscillatory trajectory has not settled.')
            if g['group']=='B1' and outcome=='E_L':
                clock = np.arange(0.,total+1.,2.)
                xx = sol.sol(clock)[0]
                peaks,_ = find_peaks(xx)
                troughs,_ = find_peaks(-xx)
                diagnostic['relaxation_envelope'] = {
                    'maxima':[[float(clock[k]),float(xx[k])] for k in peaks],
                    'minima':[[float(clock[k]),float(xx[k])] for k in troughs]}
            g['trajectories'][name] = diagnostic
            for t,u in zip(sol.t,sol.y.T):
                if t>g['head_end']:continue
                rows.append(dict(group=g['group'], trajectory=name,
                                 tau=float(t), x=float(u[0]), y=float(u[1]), z=float(u[2])))
            print(f"  {outcome}; late x=[{yy.min():.8f}, {yy.max():.8f}]", flush=True)
        g['tail_start'] = tail_start
    SERIES.parent.mkdir(parents=True,exist_ok=True)
    with SERIES.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['group','trajectory','tau','x','y','z'])
        writer.writeheader()
        writer.writerows(rows)
    meta = dict(coordinates='tau=r*t; x=C/K; y=E/r; z=Z', rho=rho,
                fixed_parameters=fixed, alpha_hopf_slow=hopf,
                solver=dict(method='DOP853', rtol=2.e-10, atol=2.e-12, max_step=6.),
                selection_sources=['data/continuation/branches.csv',
                                   'data/continuation/landmarks.csv',
                                   'data/periodic/hopf_reference.json',
                                   'data/periodic/fold_verified.json'], groups=groups)
    META.write_text(json.dumps(meta, indent=2)+'\n', encoding='utf-8')
    print('Saved only the four selected slices; original data remain unchanged.', flush=True)

if __name__=='__main__':
    simulate()
