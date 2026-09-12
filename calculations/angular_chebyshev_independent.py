"""Independent scalar angular solver: Chebyshev collocation of the ODE.

S=(1-x^2)^(|m|/2) y removes the endpoint power. Endpoints enforce the
regular singular ODE, not arbitrary Dirichlet values. No spherical basis.
For real c: branch index ell-|m| by Sturm ordering. Complex c via continuation.
"""
import numpy as np
from scipy.linalg import eig


def collocation_matrix(m,c,degree):
    j=np.arange(degree+1)
    x=np.cos(np.pi*j/degree)
    weights=np.ones(degree+1)
    weights[[0,-1]]=2
    weights*=(-1.)**j
    diff=x[:,None]-x[None,:]
    D=(weights[:,None]/weights[None,:])/(diff+np.eye(degree+1))
    D-=np.diag(D.sum(axis=1))
    m=abs(m)
    mat=-(1-x*x)[:,None]*(D@D)+2*(m+1)*x[:,None]*D
    mat=mat+np.diag(m*(m+1)-complex(c)**2*x*x)
    return x,mat


def eigenvalue(ell,m,c,degree=None,steps=16):
    if abs(m)>ell:
        raise ValueError('invalid mode')
    degree=degree or max(60,ell-abs(m)+40)
    _,mat=collocation_matrix(m,float(np.real(c)),degree)
    vals,vec=eig(mat)
    indices=np.argsort(vals.real)
    idx=indices[ell-abs(m)]
    value=vals[idx]
    if abs(value.imag)>1e-6:
        raise RuntimeError('Real-c branch not numerically real')
    v=vec[:,idx]/np.linalg.norm(vec[:,idx])
    if np.imag(c):
        for t in np.linspace(0,1,steps+1)[1:]:
            _,mat=collocation_matrix(m,np.real(c)+1j*t*np.imag(c),degree)
            vals,vec=eig(mat)
            vec/=np.linalg.norm(vec,axis=0)
            idx=np.argmax(abs(v.conj()@vec))
            value=vals[idx]; v=vec[:,idx]
    return complex(value)


def eigenvalue_mp(ell,m,c,degree=48,dps=55):
    """High-precision inverse iteration on the independent collocation matrix.

    Double collocation only initializes the branch; all matrix entries and
    the eigenpair are recomputed in arbitrary precision. No spherical solver.
    """
    import mpmath as mp
    guess=eigenvalue(ell,m,c,max(60,degree))
    _,af=collocation_matrix(m,c,degree)
    vals,vec=eig(af)
    v0=vec[:,np.argmin(abs(vals-guess))]
    with mp.workdps(dps):
        n=degree+1
        x=[mp.cos(mp.pi*j/degree) for j in range(n)]
        weights=[mp.mpf(2 if j in (0,degree) else 1)*(-1)**j for j in range(n)]
        D=mp.matrix(n)
        for i in range(n):
            for j in range(n):
                if i!=j:
                    D[i,j]=weights[i]/weights[j]/(x[i]-x[j])
            D[i,i]=-sum(D[i,j] for j in range(n) if i!=j)
        D2=D*D
        cc=mp.mpc(c)
        mat=mp.matrix(n)
        for i in range(n):
            for j in range(n):
                mat[i,j]=-(1-x[i]**2)*D2[i,j]+2*(abs(m)+1)*x[i]*D[i,j]
            mat[i,i]+=abs(m)*(abs(m)+1)-cc**2*x[i]**2
        value=mp.mpc(guess)
        v=mp.matrix([mp.mpc(z) for z in v0])
        for iteration in range(10):
            v=mp.lu_solve(mat-value*mp.eye(n),v)
            v/=mp.norm(v)
            av=mat*v
            updated=sum(mp.conj(v[j])*av[j] for j in range(n))
            residual=mp.norm(av-updated*v)
            value=updated
            if residual<mp.mpf(10)**(-dps//2):
                return complex(value),float(residual)
        raise RuntimeError(f'Inverse iteration failed, residual {residual}')


def independent_fit(ells=(28,34,40,46,52,61),chat=.6,degree=64,dps=55):
    """Unconstrained 1/L coefficient from independent integer-m collocation."""
    L=np.array(ells)+.5
    values=[]
    for ell,l in zip(ells,L):
        if (2*ell+1)%3:
            raise ValueError('Sequence is not exactly mu=2/3')
        m=(2*ell+1)//3
        val,res=eigenvalue_mp(ell,m,chat*l,degree,dps)
        values.append(val/l**2)
    scale=L.min()
    coeff=np.polynomial.polynomial.polyfit(scale/L,values,4)
    coeff*=scale**np.arange(5)
    return dict(A0=complex(coeff[0]),A1=complex(coeff[1]),A2=complex(coeff[2]),
                values=[complex(v) for v in values])


if __name__=='__main__':
    from kerr_eikonal_order_test import spheroidal_eigenvalue
    for ell,m,chat in ((10,3,.2),(20,10,.4),(40,27,.6),(61,41,.6),(40,27,.3+.2j)):
        c=chat*(ell+.5)
        ref=spheroidal_eigenvalue(ell,m,c)
        for degree in (60,90,120):
            value=eigenvalue(ell,m,c,degree)
            print(ell,m,chat,degree,value,'relative difference',abs(value/ref-1))
    print('Independent high-precision A1 fit:',independent_fit())
