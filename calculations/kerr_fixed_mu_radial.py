"""Kerr scalar radial scaling on integer-mode sequences with EXACT m/L.

Fixed real omega/L, not an on-shell complex QNM frequency sequence.
Leading angular coefficient extrapolated on a disjoint high-ell set;
the fit includes an unconstrained 1/L term (does not assume it vanishes).
"""
from fractions import Fraction
import json
import numpy as np
from scipy.optimize import root
from kerr_eikonal_order_test import spheroidal_eigenvalue


def modes(mu,targets):
    mu=Fraction(mu)
    out=[]
    for target in targets:
        candidates=[]
        for ell in range(max(1,target-10),target+11):
            m=mu*Fraction(2*ell+1,2)
            if m.denominator==1 and abs(m)<=ell:
                candidates.append((abs(ell-target),ell,int(m)))
        if not candidates:
            raise ValueError('No integer-m fixed-mu sequence near requested ell')
        _,ell,m=min(candidates)
        out.append((ell,m))
    if len(set(out))!=len(out):
        raise ValueError('Repeated modes')
    return out


def angular_fit(mu,chat,targets=(240,300,360,420),degree=3):
    pairs=modes(mu,targets)
    L=np.array([ell+.5 for ell,m in pairs])
    vals=np.array([spheroidal_eigenvalue(ell,m,chat*l).real/l**2
                   for (ell,m),l in zip(pairs,L)])
    # Scale the design columns to improve numerical conditioning.
    scale=L.min()
    coef=np.polynomial.polynomial.polyfit(scale/L,vals,degree)
    return coef*np.array([scale**k for k in range(degree+1)])


def peak_solution(a,mu,references=(240,300,360,420),degree=3):
    """Real eikonal double root, angular leading term extrapolated at fixed mu."""
    mu=Fraction(mu)
    x=np.array([1/np.sqrt(27),3.])
    final_residual=None
    for rotation in np.linspace(0,a,max(2,int(abs(a)/.1)+1)):
        def equations(values):
            om,r=values
            A0=angular_fit(mu,rotation*om,references,degree)[0]
            H=r*r+rotation**2
            D=r*r-2*r+rotation**2
            bar=A0+rotation**2*om**2-2*float(mu)*rotation*om
            shift=om-float(mu)*rotation/H
            return [shift**2-D*bar/H**2,
                    4*float(mu)*rotation*r*shift/H**2-((2*r-2)/H**2-4*r*D/H**3)*bar]
        sol=root(equations,x,tol=1e-10)
        final_residual=np.linalg.norm(equations(sol.x))
        if final_residual>1e-9:
            raise RuntimeError(f'Double root failed: {sol.message}')
        x=sol.x
    if x[1]<=1+np.sqrt(1-a*a):
        raise RuntimeError('Peak inside horizon')
    return float(x[0]),float(x[1]),float(final_residual)


def measure(a=.6,mu=Fraction(2,3),omega_hat=.225,
            references=(240,300,360,420),degree=3):
    mu=Fraction(mu)
    peak,residual=None,None
    if omega_hat is None:
        omega_hat,peak,residual=peak_solution(a,mu,references,degree)
    coeff=angular_fit(mu,a*omega_hat,references,degree)
    r=np.linspace(2.5,4.,401) if peak is None else np.linspace(peak-.2,peak+.2,401)
    if r.min()<=1+np.sqrt(1-a*a):
        raise ValueError('Measurement window crosses horizon')
    H=r*r+a*a
    delta=r*r-2*r+a*a
    B=delta*r/H**2
    dB=((2*r-2)*r+delta)/H**2-4*r*r*delta/H**3
    geometric=(delta/H)*dB+B*B
    A0=coeff[0]
    q0=(omega_hat-float(mu)*a/H)**2-delta*(A0+a*a*omega_hat**2-2*float(mu)*a*omega_hat)/H**2
    pairs=modes(mu,(40,60,80,120,160))
    eps,errors=[],[]
    for ell,m in pairs:
        L=ell+.5
        angular=spheroidal_eigenvalue(ell,m,a*omega_hat*L).real/L**2
        q=(omega_hat-(m/L)*a/H)**2-delta*(angular+a*a*omega_hat**2-2*(m/L)*a*omega_hat)/H**2-geometric/L**2
        eps.append(1/L)
        errors.append(float(np.max(abs(q-q0))))
    return dict(a=a,mu=str(mu),omega_hat=omega_hat,r_peak=peak,root_residual=residual,A0=float(A0),
                A1_fit=float(coeff[1]),slope=float(np.polyfit(np.log(eps),np.log(errors),1)[0]),
                modes=pairs,errors=errors)


if __name__=='__main__':
    for mu in (Fraction(2,5),Fraction(2,3),Fraction(4,5)):
        for a in (0.,.3,.6,.9):
            print(json.dumps(measure(a,mu)))
    print('Independent reference-range checks:')
    for refs in ((180,240,300,360),(300,360,420,480)):
        print(json.dumps(measure(references=refs)))
    print('Self-consistent real eikonal peaks:')
    for a in (0.,.3,.6,.9):
        print(json.dumps(measure(a=a,omega_hat=None)))
