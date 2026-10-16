"""Regenerate Figure 4 data. Optional --auto also repeats the two AUTO slices/planes.

The default repeats independent shooting/Floquet calculations from compact
numerical seeds; it does not require AUTO or the original project tree.
"""
import argparse,csv,json,os,shlex,shutil,subprocess
import numpy as np
from scipy.optimize import root
from scipy.linalg import eigvals
from model import ROOT,BASE,R0,RHO,equilibria,hopfs,rhs,flow,transverse,correct,write_json

DATA=ROOT/'data'; TARGET=DATA/'periodic'

def table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def fold(seed,p):
    def residual(v):
        q=dict(p,lambdaT=v[4]); f=flow(v[:3],v[3],q,8e-13,.04)
        return np.r_[f['endpoint']-v[:3],rhs(v[:3],q)[0],
                     np.linalg.det(transverse(v[:3],f['M'],q)-np.eye(2))]
    ans=root(residual,np.r_[seed['state0_ref'],seed['period_tau0'],seed['alpha0']],
             method='hybr',options={'xtol':2e-9,'eps':1e-6})
    error=float(max(abs(residual(ans.x))))
    if error>2e-8:raise RuntimeError(('Fold refinement failed',ans.message,error))
    u,T,a=ans.x[:3],ans.x[3],ans.x[4];q=dict(p,lambdaT=a)
    f=flow(u,T,q,8e-13,.04);mu=eigvals(transverse(u,f['M'],q))
    return dict(alpha0=float(a),lambdaT_per_day=float(a*R0),gamma0=p['dD'],
                period_days=float(T/R0),residual=error,
                leading_transverse_multiplier=float(mu[np.argmax(abs(mu))].real),
                u0=u.tolist(),T_tau0=float(T))

def periodic():
    seeds=json.loads((DATA/'cycle_seeds.json').read_text())
    for tag in ['stable_0p96925','saddle_0p96925','cycle_alpha_4','cycle_alpha_4p2']:
        seed=seeds[tag];p=dict(BASE,lambdaT=seed['alpha0'],dD=seed['gamma0'])
        c=correct(np.array(seed['state0_ref']),seed['period_tau0'],p)
        if c['closure']>1e-9:raise RuntimeError(('Orbit did not close',tag))
        u=c['state'];t=c['time']
        if tag.startswith('cycle_'):
            K=json.loads((DATA/'parameters.json').read_text())['analysis']['dimensional']['K']
            table(DATA/'trajectories'/(tag+'_metrics.csv'),[
                dict(alpha0=p['lambdaT'],C_min_mm3=float(K*u[:,0].min()),C_max_mm3=float(K*u[:,0].max()))])
        else:
            if (abs(c['leading'])<1)!=(tag.startswith('stable')):raise RuntimeError('Wrong periodic stability')
            np.savetxt(TARGET/(tag+'_corrected.csv'),np.c_[t,t/R0,RHO*t,u[:,0],u[:,1]/RHO,u[:,2]],
                       delimiter=',',header='tau0,days,tau,x,y,z',comments='')
            write_json(TARGET/(tag+'_corrected.json'),dict(alpha0=p['lambdaT'],gamma0=p['dD'],
                period_tau0=c['T'],period_tau=c['T']*RHO,period_days=c['T']/R0,
                leading_transverse_multiplier=c['leading'],closure=c['closure'],state0_ref=c['u0']))
        print('Corrected',tag,'closure',c['closure'],flush=True)
    rows=[]
    for seed in json.loads((DATA/'branch_seeds.json').read_text()):
        p=dict(BASE,lambdaT=seed['alpha0'],dD=.05);u=np.array(seed['state0_ref'])
        f=flow(u,seed['period_tau0'],p,2e-12,.08,1801)
        mu=eigvals(transverse(u,f['M'],p));lead=float(mu[np.argmax(abs(mu))].real)
        closure=float(max(abs(f['endpoint']-u)))
        if closure>8e-8:raise RuntimeError(('Branch closure',seed['tag'],closure))
        if seed['type']!='LP' and (abs(lead)<1)!=seed['stable_AUTO']:raise RuntimeError('Branch stability')
        row=dict(alpha0=p['lambdaT'],period_days=seed['period_tau0']/R0,
                 leading_multiplier=lead,closure=closure,stable=seed['stable_AUTO'],type=seed['type'])
        for i,name in enumerate(['x','y0','z']):
            row[name+'_min']=float(f['state'][:,i].min());row[name+'_max']=float(f['state'][:,i].max())
        rows.append(row)
    table(TARGET/'branch_metrics.csv',rows)
    write_json(TARGET/'fold_verified.json',fold(seeds['LPC'],dict(BASE,dD=.05)))
    a=hopfs(dict(BASE,dD=.05))[0][1]
    write_json(TARGET/'hopf_reference.json',dict(alpha0=a,gamma0=.05,lambdaT_per_day=a*R0))
    rows=[]
    for rate in np.linspace(.193828,.19402,60):
        for e in equilibria(dict(BASE,dD=.05,lambdaT=rate/R0)):
            rows.append(dict(lambdaT_per_day=rate,x=e['x'],stable=e['stable']))
    table(TARGET/'equilibrium_overlay.csv',rows)
    print('Regenerated Floquet branch, cycle fold and equilibrium overlay.',flush=True)

def auto():
    work=ROOT/'_work/auto';work.mkdir(parents=True,exist_ok=True)
    for name in ['model.f90','continue.auto']:shutil.copyfile(ROOT/'auto'/name,work/name)
    shutil.copyfile(DATA/'parameters.json',work/'parameters.json')
    bash=os.environ.get('AUTO_MSYS_BASH')
    if bash:
        auto_dir=os.environ.get('AUTO_DIR')
        if not auto_dir:raise RuntimeError('Set AUTO_DIR to the MSYS AUTO installation.')
        posix=work.as_posix()
        if posix[1:2]==':':posix='/'+posix[0].lower()+posix[2:]
        command='export PATH=/ucrt64/bin:/usr/bin:$PATH; source '+shlex.quote(auto_dir+'/cmds/auto.env.sh')+'; cd '+shlex.quote(posix)+'; auto continue.auto'
        cmd=[bash,'-lc',command]
    else:
        executable=shutil.which('auto')
        if not executable:raise RuntimeError('Install AUTO07p and put auto in PATH; see README.')
        cmd=[executable,'continue.auto']
    log=work/'run.log'
    with log.open('w') as f:subprocess.run(cmd,cwd=work,stdout=f,stderr=subprocess.STDOUT,check=True)
    if 'Traceback (most recent call last)' in log.read_text(errors='replace'):raise RuntimeError('AUTO failed; see '+str(log))
    # Copy only exported data; keep raw AUTO continuations and executables local.
    for folder in ['continuation','periodic']:
        for path in (work/folder).glob('*.csv'):shutil.copyfile(path,DATA/folder/path.name)
    rows=list(csv.DictReader((work/'planes.csv').open()))
    for r in rows:
        if r['family']=='LPC' and r['valid']=='True':
            p=dict(BASE,sT=float(r['sigma'])*RHO**2,lambdaD=float(r['beta'])*RHO)
            hs=hopfs(p)
            r['hopf_reference']=hs[0][1]/RHO if hs else float('nan')
            r['fold_offset']=float(r['alpha'])-r['hopf_reference']
        else:r.update(hopf_reference='',fold_offset='')
    table(DATA/'planes/comparison.csv',rows)
    print('AUTO slices and the two main-text parameter planes regenerated.',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--auto',action='store_true',help='Regenerate both AUTO slices and both parameter planes first.')
    args=parser.parse_args()
    if args.auto:auto()
    periodic()
