"""Locate zeros of the regular Frobenius amplitude, not a Riccati fit."""
import mpmath as mp
from kerr_eikonal_order_test import spheroidal_eigenvalue


def build(order=50,dps=60):
    mp.mp.dps=dps
    a=mp.mpf(3)/5; rp=mp.mpf(9)/5; gap=mp.mpf(8)/5
    m=19; om=mp.mpc(3,-.25)
    eigen=complex(spheroidal_eigenvalue(28,m,complex(a*om)))
    lam=mp.mpc(eigen.real,eigen.imag)+a*a*om**2-2*a*m*om
    hp=rp**2+a*a; s=-1j*(om-m*a/hp)/(gap/hp)
    def polys(x):
        r=rp+x; H=r*r+a*a; T=gap+x
        J=(T+x)*H-2*r*x*T
        A=x*x*T*T*H*H
        B=2*s*x*T*T*H*H+x*T*J*H
        C=s*(s-1)*T*T*H*H+s*T*J*H+(om*H-m*a)**2*H*H
        C-=x*T*lam*H*H+x*T*J*r+x*x*T*T*a*a
        return A,B,C
    abc=[mp.taylor(lambda x:polys(x)[i],0,8) for i in range(3)]
    A,B,C=abc
    indicial=abs(C[0]); C[0]=0
    f=[mp.mpc(1)]
    def get(arr,n):return arr[n] if 0<=n<len(arr) else 0
    for n in range(1,order+1):
        total=sum((get(A,n-k+2)*k*(k-1)+get(B,n-k+1)*k+get(C,n-k))*f[k]
                  for k in range(n))
        f.append(-total/(A[2]*n*(n-1)+B[1]*n))
    return f,indicial


if __name__=='__main__':
    for degree,dps in ((30,50),(50,60),(70,80)):
        f,indicial=build(degree,dps)
        fun=lambda x:mp.polyval(list(reversed(f)),x)
        zero=mp.findroot(fun,(mp.mpc(.001,.001),mp.mpc(.002,-.001)),tol=mp.mpf('1e-35'))
        deriv=mp.diff(fun,zero)
        print('degree,dps=',degree,dps,'rho_zero=',mp.nstr(zero,22),
              '|rho|=',mp.nstr(abs(zero),16),'|f|=',mp.nstr(abs(fun(zero)),5),
              '|f_prime|=',mp.nstr(abs(deriv),12),'indicial=',mp.nstr(indicial,3),flush=True)
