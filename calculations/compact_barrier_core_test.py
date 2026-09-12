"""Independent outgoing QNM countertest to universal local-WKB redundancy.

Smooth compact barriers; not a Schwarzschild/Einstein solution.
Both outgoing boundary conditions are exact since V=0 outside support.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


def bump(x,center=0.,width=1.):
    y=(np.asarray(x)-center)/width
    ans=np.zeros_like(y,dtype=float)
    inside=abs(y)<1
    ans[inside]=np.exp(1-1/(1-y[inside]**2))
    return ans


def solve_mode(eta,guess=3.2-.7j,tol=2e-11):
    def integrate(omega,dense=False):
        def rhs(x,y):
            v=10*bump(x)+eta*bump(x,2.,.4)
            return [y[1],(v-omega**2)*y[0]]
        return solve_ivp(rhs,(-1.,2.4),[1.+0j,-1j*omega],
            method='DOP853',rtol=tol,atol=tol/100,max_step=.035,
            dense_output=dense)
    def residual(v):
        w=complex(*v)
        s=integrate(w)
        if not s.success:
            raise RuntimeError(s.message)
        r=s.y[1,-1]-1j*w*s.y[0,-1]
        return [r.real,r.imag]
    r=root(residual,[guess.real,guess.imag],tol=1e-10)
    w=complex(*r.x)
    defect=np.linalg.norm(residual(r.x))
    if defect>1e-8:
        raise RuntimeError((r.message,defect))
    s=integrate(w,True)
    x=np.linspace(-.8,.8,2001)
    psi,chi=s.sol(x)
    z=chi/psi
    qm=(w*w-10*bump(x)).real-z.imag**2
    weight=np.exp(-x*x)
    den=abs(w)**2+abs(10*bump(x))+z.imag**2
    indicator=np.trapezoid(weight*abs(qm),x)/np.trapezoid(weight*den,x)
    return w,s,defect,float(indicator)


def main():
    w,s,_,_=solve_mode(0)
    x=np.linspace(-1,2.4,20001)
    psi=s.sol(x)[0]
    numerator=np.trapezoid(bump(x,2.,.4)*psi**2,x)
    denominator=2*w*np.trapezoid(psi**2,x)+1j*(psi[0]**2+psi[-1]**2)
    perturbative=numerator/denominator
    print('omega derivative, generalized norm:',perturbative,flush=True)
    print('| eta | omega | outgoing residual | integral indicator |')
    print('|---|---|---|---|',flush=True)
    for eta in (-.01,-.001,0.,.001,.01):
        wm,_,res,ind=solve_mode(eta,w)
        print(f'| {eta} | {wm} | {res:.3g} | {ind:.12g} |',flush=True)
    step=1e-4
    fd=(solve_mode(step,w)[0]-solve_mode(-step,w)[0])/(2*step)
    print('central derivative:',fd,'relative discrepancy:',abs(fd/perturbative-1),flush=True)
    for eta in (0.,.01):
        wa,_,_,ia=solve_mode(eta,w)
        wb,_,_,ib=solve_mode(eta,w,tol=2e-12)
        print('tolerance check:',eta,'delta omega:',abs(wa-wb),
              'delta indicator:',abs(ia-ib),flush=True)


if __name__=='__main__':
    main()
