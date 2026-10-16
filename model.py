"""Isolated mechanism variants in the existing fixed-clock coordinates.
RHS/Jacobian and shooting algorithm follow publication_code/common.
No import from that package: its config module creates old output folders.
"""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.optimize import brentq, root_scalar

ROOT=Path(__file__).resolve().parent
INPUT=json.loads((ROOT/'data/parameters.json').read_text())['analysis']
BASE=dict(INPUT['normalized'],dD=.05,saturated=1.)
R0=INPUT['r0_per_day']; R=INPUT['dimensional']['r']; RHO=BASE['r']
RESULTS=ROOT/'data/equilibria'
MAP={1:'sT',2:'lambdaT',3:'dT',4:'lambdaD',5:'KA',6:'dD',7:'KC',9:'r',10:'saturated'}

def variant(name,**kw):
    return dict(BASE,sT=0. if name in ('M1','M2') else BASE['sT'],
                saturated=0. if name in ('M1','M3') else 1.,**kw)

def growth(p):
    rho=p['r']
    return dict(sigma=p['sT']/rho**2,alpha=p['lambdaT']/rho,delta=p['dT']/rho,
                beta=p['lambdaD']/rho,a=p['KA']/rho,gamma=p['dD']/rho,kappa=p['KC'],n=4)

def rhs(u,p):
    x,y,z=u; h=x**4/(p['KC']**4+x**4); q=x*y
    d=p['KA']+p['saturated']*q
    return np.array([p['r']*x*(1-x)-x*y,
        p['sT']*h+(p['lambdaT']*z/(1+z)-p['dT'])*y,
        p['lambdaD']*q/d-p['dD']*z])

def jac(u,p):
    x,y,z=u; k=p['KC']**4; q=x*y
    b=p['lambdaD']*p['KA']/(p['KA']+p['saturated']*q)**2
    return np.array([[p['r']*(1-2*x)-y,-x,0],
        [p['sT']*4*k*x**3/(k+x**4)**2,p['lambdaT']*z/(1+z)-p['dT'],p['lambdaT']*y/(1+z)**2],
        [b*y,b*x,-p['dD']]])

def state(x,p):
    y=p['r']*(1-x);q=x*y
    z=p['lambdaD']/p['dD']*q/(p['KA']+p['saturated']*q)
    return np.array([x,y,z])

def alpha_at(x,p):
    u=state(x,p); y,z=u[1:]
    h=x**4/(p['KC']**4+x**4)
    return (p['dT']-p['sT']*h/y)*(1+z)/z

def balance_poly(p):
    x=P([0,1]);y=p['r']*(1-x);q=x*y
    den=p['KA']+(p['saturated']+p['lambdaD']/p['dD'])*q
    return p['sT']*x**4*den+(p['lambdaT']*p['lambdaD']/p['dD']*q-p['dT']*den)*y*(p['KC']**4+x**4)

def equilibria(p,boundary=True):
    if p['sT']==0:
        # Immune-active roots solve q=a*delta/[b/gamma*(alpha-delta)-sat*delta].
        den=p['lambdaD']/p['dD']*(p['lambdaT']-p['dT'])-p['saturated']*p['dT']
        disc=1-4*p['KA']*p['dT']/den/p['r'] if den>0 else -1
        xs=[] if disc<-1e-12 else [(1-np.sqrt(max(0,disc)))/2,(1+np.sqrt(max(0,disc)))/2]
    else:
        roots=balance_poly(p).roots()
        xs=sorted(float(v.real) for v in roots if abs(v.imag)<2e-6 and 1e-10<v.real<1-1e-10)
        # Correct the rational balance, rather than trusting polynomial residual scaling.
        out=[]
        for x in xs:
            ans=root_scalar(lambda t:rhs(state(t,p),p)[1],x0=x,x1=x+1e-7,method='secant',xtol=3e-15)
            if ans.converged and 0<ans.root<1 and all(abs(ans.root-v)>1e-7 for v in out):out.append(ans.root)
        xs=out
    states=[('internal',state(x,p)) for x in xs]
    if boundary:
        states.insert(0,('boundary_origin',np.zeros(3)))
        if p['sT']==0:states.append(('boundary_tumour',np.array([1.,0.,0.])))
    out=[]
    for kind,u in states:
        ev=np.linalg.eigvals(jac(u,p)); ab=float(ev.real.max())
        out.append(dict(kind=kind,x=float(u[0]),y=float(u[1]/p['r']),z=float(u[2]),
                        state=u.tolist(),abscissa=ab,stable=ab<-1e-9,
                        eigenvalues=[[float(v.real),float(v.imag)] for v in ev],
                        residual=float(np.max(abs(rhs(u,p))))))
    return out

def rh(x,p):
    a=alpha_at(x,p); c=np.poly(jac(state(x,p),dict(p,lambdaT=a)))
    return c[1]*c[2]-c[3]

def hopfs(p):
    grid=np.linspace(.015,.49,450); roots=[]
    for l,r in zip(grid[:-1],grid[1:]):
        if rh(l,p)*rh(r,p)<0:
            x=brentq(lambda x:rh(x,p),l,r,xtol=1e-14)
            a=alpha_at(x,p);ev=np.linalg.eigvals(jac(state(x,p),dict(p,lambdaT=a)))
            if a>0 and np.max(abs(ev.imag))>1e-6:roots.append((x,a))
    return roots

def write_json(path,obj):
    def convert(value):
        if isinstance(value,np.generic):return value.item()
        if isinstance(value,np.ndarray):return value.tolist()
        raise TypeError(type(value).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=convert),encoding='utf-8')

from scipy.integrate import solve_ivp
from scipy.linalg import solve,eigvals,null_space

def flow(u,T,p,rtol=2e-12,max_step=.08,samples=0):
    def augmented(t,v):
        state=v[:3];J=jac(state,p);M=v[3:12].reshape(3,3)
        return np.r_[rhs(state,p),(J@M).ravel(),np.trace(J)]
    solution=solve_ivp(augmented,[0,T],np.r_[u,np.eye(3).ravel(),0.],method='DOP853',rtol=rtol,atol=rtol/100,max_step=max_step,dense_output=samples>0)
    assert solution.success and np.all(np.isfinite(solution.y))
    v=solution.y[:,-1]; out={'endpoint':v[:3],'M':v[3:12].reshape(3,3),'logdet':v[12]}
    if samples:
        out['time']=np.linspace(0,T,samples);out['state']=solution.sol(out['time'])[:3].T
    return out

def transverse(u,M,p):
    normal=np.array([p['r'],1.,0.]); f=rhs(u,p)
    Q=null_space(normal[None,:]); P=np.eye(3)-np.outer(f,normal)/np.dot(normal,f)
    return Q.T@P@M@Q

def correct(u,T,p,rtol=2e-12,max_step=.08):
    for iteration in range(15):
        result=flow(u,T,p,rtol,max_step);f0=rhs(u,p);fT=rhs(result['endpoint'],p)
        residual=np.r_[result['endpoint']-u,f0[0]]
        if max(abs(residual))<4e-11:break
        A=np.block([[result['M']-np.eye(3),fT[:,None]],[jac(u,p)[0][None,:],np.zeros((1,1))]])
        step=solve(A,-residual)
        damping=1.
        while damping>1/1024:
            ut=u+damping*step[:3];Tt=T+damping*step[3]
            if min(ut)>0 and Tt>0:
                trial=flow(ut,Tt,p,rtol,max_step)
                trial_res=np.r_[trial['endpoint']-ut,rhs(ut,p)[0]]
                if max(abs(trial_res))<max(abs(residual)):
                    u,T=ut,Tt;break
            damping/=2
        else:raise RuntimeError('Shooting line search failed')
    else:raise RuntimeError('Shooting iteration limit')
    result=flow(u,T,p,rtol,max_step,2401);M=result['M'];P=transverse(u,M,p)
    mu=eigvals(P);lead=mu[np.argmax(abs(mu))]
    result.update(u0=u,T=T,p=p,closure=float(max(abs(result['endpoint']-u))),
                  phase=float(abs(rhs(u,p)[0])),leading=float(lead.real),
                  small=float(np.exp(result['logdet'])/lead.real),
                  trivial_error=float(np.linalg.norm(M@rhs(u,p)-rhs(u,p),np.inf)/np.linalg.norm(rhs(u,p),np.inf)))
    return result
