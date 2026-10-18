"""Local Hopf diagnostics for the two M4/LLC Figure 4 slices.

Reuses model.py and the Hessian/normal-form formula of the original local
diagnostics. Only equilibrium roots and small linear systems are evaluated.
Run: python diagnose_hopf.py
"""
import json
import platform

import numpy as np
import scipy
from numpy.linalg import eig, solve

from model import ROOT, BASE, INPUT, R0, RHO, equilibria, hopfs, jac, rhs, state


def hessian(u, p):
    """Unscaled second derivatives in the existing (x, y0, z) coordinates."""
    x, y, z = u
    k, q, s = p['KC']**4, x*y, p['saturated']
    first = p['lambdaD']*p['KA']/(p['KA']+s*q)**2
    second = -2*p['lambdaD']*p['KA']*s/(p['KA']+s*q)**3
    hh = 4*k*x**2*(3*k-5*x**4)/(k+x**4)**3
    H = np.zeros((3, 3, 3))
    H[0] = [[-2*p['r'], -1, 0], [-1, 0, 0], [0, 0, 0]]
    H[1] = [[p['sT']*hh, 0, 0],
            [0, 0, p['lambdaT']/(1+z)**2],
            [0, p['lambdaT']/(1+z)**2, -2*p['lambdaT']*y/(1+z)**3]]
    H[2] = [[second*y*y, first+second*q, 0],
            [first+second*q, second*x*x, 0], [0, 0, 0]]
    return H


def eigenpair(A):
    values, vectors = eig(A)
    index = np.argmax(values.imag)
    lam, q = values[index], vectors[:, index]
    if lam.imag <= 0:
        raise RuntimeError('No positive-frequency eigenvalue pair.')
    q = q/np.linalg.norm(q)
    left_values, left_vectors = eig(A.T)
    p = left_vectors[:, np.argmin(abs(left_values-np.conj(lam)))]
    p = p/np.conj(np.vdot(p, q))
    return lam, q, p


def normal_form(A, H, D):
    """Kuznetsov convention: B=D2f, C=D3f, l1=Re(G21)/(2*omega)."""
    lam, q, p = eigenpair(A)
    omega = lam.imag
    B = lambda v, w: np.einsum('kij,i,j->k', H, v, w)
    C = lambda v, w, a: np.einsum('kijl,i,j,l->k', D, v, w, a)
    term = C(q, q, np.conj(q))-2*B(q, solve(A, B(q, np.conj(q))))
    term += B(np.conj(q), solve(2j*omega*np.eye(3)-A, B(q, q)))
    G21 = np.vdot(p, term)
    return dict(l1=float(G21.real/(2*omega)), omega=float(omega),
                G21=[float(G21.real), float(G21.imag)],
                normalization_error=float(abs(np.vdot(p, q)-1)),
                eigenvector_residual=float(np.linalg.norm(A@q-lam*q, np.inf)))


def first_lyapunov(u, p, relative_step=1e-5):
    """Return l1 in both coordinate conventions; positive time scales only.

    Third derivatives use centred differences of the analytic Hessian with
    h_j=relative_step*(1+abs(u_j)) in fixed-clock coordinates. If u=S*v,
    tau=rho*tau0, transform each derivative tensor as S^-1 D^k f[S,...,S]/rho.
    Eigenvectors are separately normalized in each coordinate convention.
    """
    A, H = jac(u, p), hessian(u, p)
    D = np.zeros((3, 3, 3, 3))
    for j in range(3):
        h = relative_step*(1+abs(u[j]))
        increment = np.zeros(3)
        increment[j] = h
        D[:, :, :, j] = (hessian(u+increment, p)-hessian(u-increment, p))/(2*h)
    rho = p['r']
    scale = np.array([1., rho, 1.])
    Ag = A*scale[None, :]/scale[:, None]/rho
    Hg = H*scale[None, :, None]*scale[None, None, :]/scale[:, None, None]/rho
    Dg = (D*scale[None, :, None, None]*scale[None, None, :, None]
          *scale[None, None, None, :]/scale[:, None, None, None]/rho)
    return dict(fixed_clock=normal_form(A, H, D), growth_clock=normal_form(Ag, Hg, Dg))


def hopf_diagnostic(p):
    """Compute a Hopf point, total alpha derivative, l1 and local checks."""
    points = hopfs(p)
    if len(points) != 1:
        raise RuntimeError(f'Expected one low-tumour Hopf point, found {points}.')
    x, alpha0 = points[0]
    p = dict(p, lambdaT=alpha0)
    u, A, H = state(x, p), jac(state(x, p), p), hessian(state(x, p), p)
    lam, q, left = eigenpair(A)
    # J*u_alpha0 + f_alpha0 = 0, including the movement of the equilibrium.
    f_alpha0 = np.array([0., u[1]*u[2]/(1+u[2]), 0.])
    u_alpha0 = solve(A, -f_alpha0)
    J_alpha0 = np.zeros((3, 3))
    J_alpha0[1, 1] = u[2]/(1+u[2])
    J_alpha0[1, 2] = u[1]/(1+u[2])**2
    total_J_alpha0 = J_alpha0+np.einsum('kij,j->ki', H, u_alpha0)
    crossing = float(np.vdot(left, total_J_alpha0@q).real)
    fd_checks = []
    for h in (2e-5, 1e-5, 5e-6):
        critical_real = []
        for sign in (-1, 1):
            shifted = dict(p, lambdaT=alpha0+sign*h)
            eq = min(equilibria(shifted, boundary=False), key=lambda e: abs(e['x']-x))
            ev = np.array([complex(*v) for v in eq['eigenvalues']])
            critical_real.append(float(ev[np.argmin(abs(ev-lam))].real))
        slope = (critical_real[1]-critical_real[0])/(2*h)
        fd_checks.append(dict(step_alpha0=h, step_alpha=h/p['r'], derivative=slope,
                              analytic_difference=abs(slope-crossing)))
    # An independent Jacobian difference checks the analytic Hessian.
    Hfd = np.zeros_like(H)
    for j in range(3):
        h = 1e-5*(1+abs(u[j]))
        increment = np.zeros(3)
        increment[j] = h
        Hfd[:, :, j] = (jac(u+increment, p)-jac(u-increment, p))/(2*h)
    hessian_error = float(np.max(abs(H-Hfd)/(1+abs(H))))
    steps = (2e-5, 1e-5, 5e-6)
    lyapunov_checks = [dict(relative_step=h, **first_lyapunov(u, p, h)) for h in steps]
    central = lyapunov_checks[1]
    spans = {clock: float(np.ptp([v[clock]['l1'] for v in lyapunov_checks]))
             for clock in ('fixed_clock', 'growth_clock')}
    coefficients = np.poly(A)
    residual = float(np.linalg.norm(rhs(u, p), np.inf))
    rh_residual = float(abs(coefficients[1]*coefficients[2]-coefficients[3]))
    eigenvalues = np.linalg.eigvals(A)
    third = float(eigenvalues[np.argmin(abs(eigenvalues.imag))].real)
    if (residual > 1e-10 or rh_residual > 1e-10 or abs(lam.real) > 1e-10
            or min(coefficients[1:]) <= 0 or third >= 0 or abs(crossing) < 1e-7
            or hessian_error > 1e-7
            or max(v['analytic_difference'] for v in fd_checks) > 1e-6):
        raise RuntimeError('Hopf spectral, derivative or Hessian check failed.')
    for clock in spans:
        value = central[clock]['l1']
        if abs(value) < 1e-7 or spans[clock] > 1e-6*max(1., abs(value)):
            raise RuntimeError('First Lyapunov coefficient is unresolved.')
        if any(np.sign(v[clock]['l1']) != np.sign(value) for v in lyapunov_checks):
            raise RuntimeError('First Lyapunov coefficient changes sign with step.')
        if central[clock]['normalization_error'] > 1e-10:
            raise RuntimeError('Eigenvector normalization failed.')
    l1 = central['growth_clock']['l1']
    if np.sign(l1) != np.sign(central['fixed_clock']['l1']):
        raise RuntimeError('Criticality changed under coordinate conversion.')
    growth_state = u/np.array([1., p['r'], 1.])
    return dict(alpha=alpha0/p['r'], gamma=p['dD']/p['r'],
                alpha0=alpha0, gamma0=p['dD'], equilibrium_growth=growth_state.tolist(),
                parameters_fixed_clock=p,
                crossing_derivative=crossing, crossing_finite_difference_checks=fd_checks,
                l1_growth=l1, l1_fixed_clock=central['fixed_clock']['l1'],
                l1_step_span=spans, lyapunov_step_checks=lyapunov_checks,
                omega_growth=lam.imag/p['r'], omega_fixed_clock=lam.imag,
                third_eigenvalue_growth=third/p['r'], third_eigenvalue_fixed_clock=third,
                rhs_residual_fixed_clock=residual, rh_residual_fixed_clock=rh_residual,
                critical_real_fixed_clock=float(lam.real), hessian_check_error=hessian_error,
                classification='supercritical' if l1 < 0 else 'subcritical',
                periodic_branch_side='alpha > alpha_H' if -l1/crossing > 0 else 'alpha < alpha_H',
                emerging_cycle_stability='stable' if l1 < 0 else 'saddle')


def main():
    results = {name: hopf_diagnostic(dict(BASE, dD=gamma0)) for name, gamma0 in
               [('fast', INPUT['normalized']['dD']), ('slow', BASE['dD'])]}
    report = dict(
        model='M4; LLC reference; only expansion varies within each DC-decay slice',
        coordinates=dict(growth='tau=r*(t-t1), state=(x,y,z), alpha=lambda_T/r',
                         fixed_clock='tau0=r0*(t-t1), state=(x,y0,z), alpha0=lambda_T/r0',
                         r0=R0, rho=RHO, y0='rho*y'),
        method=dict(
            hopf='Existing model.hopfs: scalar Routh-Hurwitz root and analytic Jacobian',
            crossing='Re<p,(f_u_alpha0 + D_u J * u_alpha0)q>; u_alpha0=-J^-1 f_alpha0',
            crossing_scaling='d Re(lambda_tau)/d alpha = d Re(lambda_tau0)/d alpha0 at fixed rho',
            lyapunov='l1=Re<p,C(q,q,qbar)-2B(q,J^-1 B(q,qbar))+B(qbar,(2*i*omega*I-J)^-1 B(q,q))>/(2*omega)',
            tensors='B=D2f analytic; C=D3f from centred differences of B in fixed-clock coordinates; tensors transformed for growth-clock result',
            normalization='Jq=i*omega*q; J.T*p=-i*omega*p; norm(q,2)=1; conjugate(p).T*q=1, separately in each coordinate convention',
            sign='l1<0: supercritical/stable small cycle; l1>0: subcritical/saddle small cycle when the third eigenvalue is negative',
            scope='Floating-point local nondegeneracy diagnostics; step spans are numerical checks, not rigorous error bounds'),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        results=results)
    path = ROOT/'data/periodic/hopf_diagnostics.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for name, result in results.items():
        print(f"{name}: alpha_H={result['alpha']:.12g}, gamma={result['gamma']:.12g}, "
              f"crossing={result['crossing_derivative']:.12g}, "
              f"l1_growth={result['l1_growth']:.12g}, "
              f"l1_fixed_clock={result['l1_fixed_clock']:.12g}, {result['classification']}")
    print('Saved', path.relative_to(ROOT))


if __name__ == '__main__':
    main()
