"""Analytic spatial and frequency derivatives of the outgoing WKB map.

No finite-difference mesh is used. This is not a QNM eigenvalue condition.
The outgoing square-root branch is appropriate outside the barrier only.
"""
from functools import lru_cache
import json
import numpy as np
import sympy as sp
from vaidya_solvability import leaver_qnm


@lru_cache(maxsize=2)
def evaluator(order=2):
    if order not in (0,2):
        raise ValueError('order must be 0 or 2')
    r,w,L,s=sp.symbols('r w L s')
    f=1-2/r
    V=f*((L**2-sp.Rational(1,4))/r**2+2*(1-s**2)/r**3)
    q=(w**2-V)/L**2
    D=lambda a:f*sp.diff(a,r)
    u=sp.sqrt(q)
    if order==2:
        u+=(5*D(q)**2-4*q*D(D(q)))/(32*L**2*q**sp.Rational(5,2))
    ux=D(u)
    uxx=D(ux)
    expressions=[u,ux,uxx]+[sp.diff(a,w) for a in (u,ux,uxx)]
    return sp.lambdify((r,w,L,s),expressions,'numpy',cse=True)


def profile(radius,omega,ell=70,spin=2,order=2):
    L=ell+.5
    u,ux,uxx,uw,uxw,uxxw=evaluator(order)(np.asarray(radius,dtype=float),omega,L,spin)
    v=ux/u
    z=uxx/u-v*v
    B=-.5*np.real(v)-L*np.imag(u)
    C=-.5*np.real(z)-L*np.imag(ux)
    Q=-(C+B*B)/L**2
    vw=uxw/u-ux*uw/u**2
    zw=uxxw/u-uxx*uw/u**2-2*v*vw
    def tangent(direction):
        dB=-.5*np.real(vw*direction)-L*np.imag(uw*direction)
        dC=-.5*np.real(zw*direction)-L*np.imag(uxw*direction)
        return -(dC+2*B*dB)/L**2
    return Q,tangent(1),tangent(1j)


def summary(ell=70,spin=2,window=(20.,50.),points=3001,order=2):
    omega=leaver_qnm(ell,0,spin)
    # Uniform tortoise sampling, same measure as the previous audit.
    from scipy.special import lambertw
    from vaidya_solvability import tortoise
    t=np.linspace(float(tortoise(window[0])),float(tortoise(window[1])),points)
    r=2*(1+lambertw(np.exp(t/2-1)).real)
    Q,da,db=profile(r,omega,ell,spin,order)
    med=lambda a:float(np.median(a))
    gram=np.array([[np.trapezoid(da*da,t),np.trapezoid(da*db,t)],
                   [np.trapezoid(da*db,t),np.trapezoid(db*db,t)]])
    norm=np.sqrt(np.trapezoid(Q*Q,t))
    eigenvalues=np.linalg.eigvalsh(gram)
    return dict(ell=ell,spin=spin,window=window,order=order,
        radial=med(abs((omega.real*da+omega.imag*db)/Q)),
        component_real=med(abs(omega.real*da/Q)),
        component_imag=med(abs(omega.imag*db/Q)),
        pointwise_normwise=med(abs(omega)*np.hypot(da,db)/abs(Q)),
        profile_L2_normwise=float(abs(omega)*np.sqrt(max(eigenvalues[-1],0))/norm),
        profile_L2_radial=float(np.sqrt(np.trapezoid((omega.real*da+omega.imag*db)**2,t))/norm),
        asymptotic=2*abs(omega)/abs(omega.imag))


if __name__=='__main__':
    for ell in (20,40,70,100):
        print(json.dumps(summary(ell)))
    for spin in (0,1,2):
        print(json.dumps(summary(spin=spin)))
    for window in ((10.,25.),(30.,70.),(200.,300.),(2.5,4.)):
        for order in (0,2):
            print(json.dumps(summary(window=window,order=order)))
