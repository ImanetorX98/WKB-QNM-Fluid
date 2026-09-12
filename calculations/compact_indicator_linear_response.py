"""Analytic first variation of the integral indicator; no new spectral fit."""
import numpy as np
from scipy.integrate import cumulative_simpson,simpson
from compact_barrier_core_test import solve_mode,bump


def response(points=16001,indicator_points=None):
    omega,solution,_,_=solve_mode(0.,tol=2e-12)
    whole=np.linspace(-1,2.4,points)
    psi=solution.sol(whole)[0]
    domega=simpson(bump(whole,2.,.4)*psi**2,x=whole)/(
        2*omega*simpson(psi**2,x=whole)+1j*(psi[0]**2+psi[-1]**2))
    # Include x=-1 to fix the ingoing frequency derivative exactly.
    left_points=points//9
    x=np.concatenate((np.linspace(-1,-.8,left_points,endpoint=False),
                      np.linspace(-.8,.8,indicator_points or points-left_points)))
    psi,chi=solution.sol(x)
    z=chi/psi
    integ=cumulative_simpson(psi**2,x=x,initial=0)
    k=(-1j*psi[0]**2-2*omega*integ)/psi**2
    h=domega*k
    qm=(omega**2-10*bump(x)).real-z.imag**2
    dqm=2*(omega*domega).real-2*z.imag*h.imag
    den=abs(omega)**2+10*bump(x)+z.imag**2
    dden=2*(omega.conjugate()*domega).real+2*z.imag*h.imag
    mask=x>=-.8
    weight=np.exp(-x[mask]**2)
    def integrate(f):
        rule=np.trapezoid if indicator_points else simpson
        return rule(weight*f[mask],x=x[mask])
    numerator=integrate(abs(qm))
    denominator=integrate(den)
    slope=(integrate(np.sign(qm)*dqm)*denominator-
           numerator*integrate(dden))/denominator**2
    return omega,domega,float(slope)


if __name__=='__main__':
    for n in (8001,16001,32001):
        omega,domega,slope=response(n)
        print(f'N={n}, omega_eta={domega}, E_eta={slope:.12g}',flush=True)
    omega,domega,slope=response(32001,indicator_points=2001)
    print(f'Same discrete indicator as original script: E_eta={slope:.12g}',flush=True)
    for step in (1e-3,3e-4,1e-4):
        plus=solve_mode(step,omega,tol=2e-12)[3]
        minus=solve_mode(-step,omega,tol=2e-12)[3]
        measured=(plus-minus)/(2*step)
        print(f'h={step:g}, E_eta_FD={measured:.12g}, '
              f'relative_error={abs(measured/slope-1):.4g}',flush=True)
