"""Equilibrium branches and SN2 asymptotics; copied scalar algorithms.
No trajectories, parameter fitting or AUTO installation are required.
"""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from model import ROOT,BASE,R0,variant,state,alpha_at,jac,rhs,equilibria,hopfs

def save(folder,name,rows):
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def components(u,p):
    x=1-u;s=p['sT'];rho=p['r'];k=p['KC']**4
    h=x**4/(k+x**4);hp=4*k*x**3/(k+x**4)**2
    hpp=4*k*x*x*(3*k-5*x**4)/(k+x**4)**3
    D=p['dD']*p['KA']/(p['lambdaD']*rho);q=u*(1-u)
    L=p['dT']-s*h/(rho*u)
    Lu=s/rho*(hp/u+h/u**2)
    Luu=s/rho*(-hpp/u-2*hp/u**2-2*h/u**3)
    B=1+p['dD']*p['saturated']/p['lambdaD']+D/q
    Bu=-D*(1-2*u)/q**2
    Buu=2*D/q**2+2*D*(1-2*u)**2/q**3
    return L*B,Lu*B+L*Bu,Luu*B+2*Lu*Bu+L*Buu

def sn2(p):
    h1=1/(1+p['KC']**4)
    t0=2*h1/(p['r']*p['dT'])
    t=brentq(lambda t:components(p['sT']*t,p)[1],t0*.5,t0*2.5,xtol=1e-14)
    u=p['sT']*t;a,da,dda=components(u,p)
    return 1-u,a,dda

def asymptotics(p):
    h1=1/(1+p['KC']**4);hp1=4*p['KC']**4/(1+p['KC']**4)**2
    D=p['dD']*p['KA']/(p['lambdaD']*p['r'])
    C=p['dD']*p['KA']*p['dT']**2/(4*p['lambdaD']*h1)
    C0=p['dT']/2*(1+p['dD']*p['saturated']/p['lambdaD']+D*(1+hp1/h1))
    return C,C0,2*h1/(p['r']*p['dT'])

def main():
    target=ROOT/'data/equilibria'
    curves=[];summary=[];landmarks=[];folds=[];common=[]
    grid=np.unique(np.r_[np.geomspace(.003,.05,150),np.linspace(.05,.98,800),np.linspace(.98,.9997,700)])
    for name in ('M1','M2','M3','M4'):
        p=variant(name)
        derivative=lambda x:np.imag(alpha_at(x+1e-20j,p))/1e-20
        x=brentq(derivative,.45,.55,xtol=1e-14);a=float(alpha_at(x,p))
        summary.append(dict(model=name,lambdaT_SN_per_day=a*R0))
        landmarks.append(dict(model=name,kind='equilibrium',type='LP',lambdaT_per_day=a*R0,x_min=x))
        for x,a in hopfs(p):landmarks.append(dict(model=name,kind='equilibrium',type='HB',lambdaT_per_day=a*R0,x_min=x))
        for x in grid:
            a=float(alpha_at(x,p));rate=a*R0
            if not .10<=rate<=1.25:continue
            u=state(x,p);ev=np.linalg.eigvals(jac(u,dict(p,lambdaT=a)));ab=float(ev.real.max())
            curves.append(dict(model=name,lambdaT_per_day=rate,x=x,y=u[1]/BASE['r'],z=u[2],kind='internal',stable=ab<0,abscissa=ab))
        for rate in np.linspace(.10,1.25,70):
            curves.append(dict(model=name,lambdaT_per_day=rate,x=0.,y=0.,z=0.,kind='boundary_origin',stable=False,abscissa=BASE['r']))
            if name in ('M1','M2'):curves.append(dict(model=name,lambdaT_per_day=rate,x=1.,y=0.,z=0.,kind='boundary_tumour',stable=True,abscissa=-BASE['dD']))
        if name in ('M3','M4'):
            x,a,dda=sn2(p);folds.append(dict(model=name,lambdaT_upper_SN_per_day=a*R0,x=x,scalar_second_derivative=dda))
        for e in equilibria(dict(p,lambdaT=.96925)):
            common.append(dict(model=name,kind=e['kind'],x=e['x'],y=e['y'],z=e['z'],stable=e['stable']))
    for name,rows in [('global_equilibrium_branches.csv',curves),('ablation_summary.csv',summary),('ablation_landmarks.csv',landmarks),('upper_fold_checks.csv',folds),('common_parameter_equilibria.csv',common)]:save(target,name,rows)
    rows=[]
    for name in ('M3','M4'):
        p=variant(name);C,C0,t0=asymptotics(p)
        for s in (.032,.016,.008,.004,.002,.001,.0005,.00025):
            pp=dict(p,sT=s);x,a,dda=sn2(pp)
            rows.append(dict(model=name,S_per_day2=s*R0**2,sigma0=s,x_SN2=x,lambdaT_SN2_per_day=a*R0,
                leading_prediction_per_day=R0*C/s,constant_corrected_per_day=R0*(C/s+C0),C_exact=C,C0_exact=C0))
    save(ROOT/'data/sn2','sn2_weak_recruitment.csv',rows)
    print('Regenerated Figure 2 equilibrium data and 16 Figure 3 fold points.')

if __name__=='__main__':main()
