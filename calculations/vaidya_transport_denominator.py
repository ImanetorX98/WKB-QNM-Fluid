"""Compare a regularized EF derivative overlap with the QNM generalized norm.

P = Reg int mu Z Z_r dr. Checks P=i*omega*Norm, in M=1 conventions.
This does not yet establish the full time-dependent solvability numerator.
"""
import json
import math
import numpy as np
from scipy.integrate import simpson
from vaidya_numerator_factored import (Asymptotic,inner_solve,outer_solve,
    horizon_series,horizon_values)
from vaidya_solvability import tortoise,leaver_qnm
from asymptotic_series import Series,inverse_lapse,antiderivative_series


def outgoing_Domega(ell,spin,w,r,terms):
    a=np.zeros(terms,dtype=complex)
    aw=np.zeros(terms,dtype=complex)
    a[0]=1
    for k in range(1,terms):
        p=k*(k-1)-ell*(ell+1)
        b=2*k*(k-2)+2*(1-spin**2)
        a[k]=(p*a[k-1]-(b*a[k-2] if k>1 else 0))/(2j*w*k)
        aw[k]=(p*aw[k-1]-(b*aw[k-2] if k>1 else 0))/(2j*w*k)-a[k]/w
    k=np.arange(terms)
    powers=float(r)**(-k)
    u,uw=a@powers,aw@powers
    up,upw=(-k*a)@(powers/r),(-k*aw)@(powers/r)
    return 1j+(1-2/r)*(upw/u-up*uw/u**2)


def horizon_p(ell,spin,w,width,terms):
    h=horizon_series(ell,spin,w,terms)
    index=np.arange(terms)
    dh=np.r_[index[1:]*h[1:],0.]
    amp=np.convolve(h,dh)[:terms]
    exp=(-2j*w)**index/np.array([float(math.factorial(int(k))) for k in index])
    coef=np.convolve(amp,exp)[:terms]
    power=index+1-4j*w
    return np.exp(-4j*w+4j*w*np.log(2))*np.sum(coef*width**power/power)


def compute(ell=2,spin=0,width=.5,match=25.,terms=20,uppers=(40.,50.,60.,70.),points=80001):
    w=leaver_qnm(ell,0,spin)
    asy=Asymptotic(ell,spin,w,terms)
    ri,h,dh,_=inner_solve(ell,spin,w,width,match,points,60)
    u,_,_,f,_,ulog,_,_=asy.at(match)
    phase=np.exp(-2j*w*tortoise(match))
    g0=phase*h[-1]/u[0]
    dg0=phase*(dh[-1]-h[-1]*(2j*w/f[0]+ulog[0]))/u[0]
    ro,g,dg=outer_solve(asy,match,95.,g0,dg0,points)
    u,up,_,f,_,_,_,_=asy.at(ro)
    C=g[np.searchsorted(ro,85.)]
    R2i=np.exp(-2j*w*tortoise(ri))*h*h
    R2o=np.exp(2j*w*tortoise(ro))*(u*g)**2
    pi=np.exp(-2j*w*tortoise(ri))*h*dh
    po=np.exp(2j*w*tortoise(ro))*u*g*(2j*w*u*g/f+up*g+u*dg)
    W=(asy.u*(2j*w*(asy.u*inverse_lapse(terms))+asy.u1)).truncate(terms)
    V=antiderivative_series(W,w,terms,terms)
    F=C*C*V.evaluate(ro)*np.exp(2j*w*tortoise(ro))
    pbase=horizon_p(ell,spin,w,width,60)+simpson(pi,x=ri)
    nbase=simpson(R2i/(1-2/ri),x=ri)
    def Dminus(om):
        hv,hp=horizon_values(ell,spin,om,width,60)
        return -1j*om+(1-2/(2+width))*hp/hv
    dw=1e-6*abs(w)
    Dm=(Dminus(w+dw)-Dminus(w-dw))/(2*dw)
    rows=[]
    for upper in uppers:
        j=np.searchsorted(ro,upper)
        P=pbase+simpson(po[:j+1],x=ro[:j+1])-F[j]
        Dp=outgoing_Domega(ell,spin,w,float(ro[j]),terms)
        norm=nbase+simpson(R2o[:j+1]/f[:j+1],x=ro[:j+1])+(Dp*R2o[j]-Dm*R2i[0])/(2*w)
        rows.append(dict(upper=upper,P=[P.real,P.imag],norm=[norm.real,norm.imag],
                         identity_error=float(abs(P/(1j*w*norm)-1))))
    return rows


if __name__=='__main__':
    for row in compute():
        print(json.dumps(row))
    from vaidya_numerator_factored import numerator
    result=numerator(uppers=(40.,))
    p=complex(*compute(uppers=(40.,))[0]['P'])
    K=result['N']/(2*p)
    print('Formal transport candidate K:',K)
    for radius in (3.,5.,10.):
        _,h,dh,_=inner_solve(2,0,result['omega'],.5,radius,1001,60)
        shape=-radius*dh[-1]/h[-1]
        print('radius, dM log G, Xi candidate:',radius,shape,K+shape)
