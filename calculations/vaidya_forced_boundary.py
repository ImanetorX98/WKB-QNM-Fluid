"""Independent forced-ODE test of the first-order boundary solvability.

L H = S - 2 K G_r; horizon-analytic H with H(2)=0 fixes its homogeneous gauge.
Compare its bilinear Wronskian to the outgoing particular primitive.
No claim of a full time-dependent waveform boundary condition.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp
from vaidya_numerator_factored import horizon_series,numerator
from vaidya_transport_denominator import compute
from vaidya_solvability import tortoise
from asymptotic_series import Series,integrand_series,antiderivative_series,inverse_lapse
from outgoing_asymptotics import series_coefficients


def correction_series(w,K,terms=60):
    h=horizon_series(2,0,w,terms)
    n=np.arange(terms)
    dh=np.r_[n[1:]*h[1:],0.]
    ddh=np.r_[n[1:-1]*n[2:]*h[2:],0.,0.]
    forcing=2*((1-K)*dh+2*ddh+np.r_[0.,ddh[:-1]])
    forcing=np.convolve(forcing,[8,12,6,1])[:terms]
    p=[4-16j*w,2-24j*w,-12j*w,-2j*w]
    H=np.zeros(terms,dtype=complex)
    for n in range(terms-1):
        rhs=(4*n*(n-1)-14)*H[n]
        if n>=1:
            rhs+=((n-1)*(n-2)-6)*H[n-1]
        for j in (1,2,3):
            if n-j+1>=0:
                rhs+=p[j]*(n-j+1)*H[n-j+1]
        H[n+1]=(forcing[n]-rhs)/(4*(n+1)*(n+1-4j*w))
    return h,H


def run(shift=0j,width=.5,rtol=1e-12,uppers=(30.,35.,40.,45.)):
    num=numerator(uppers=(40.,))
    w=num['omega']; C=num['C']
    P=complex(*compute(uppers=(40.,))[0]['P'])
    K=num['N']/(2*P)+shift
    hc,Hc=correction_series(w,K)
    def at(c):
        return np.polynomial.polynomial.polyval(width,c),np.polynomial.polynomial.polyval(width,np.arange(1,len(c))*c[1:])
    h,dh=at(hc); H,dH=at(Hc)
    def rhs(r,y):
        g,dg,z,dz=y
        f=1-2/r
        U=6/r**2+2/r**3
        drift=2/r**2-2j*w
        ddg=(U*g-drift*dg)/f
        forcing=2*(dg+r*ddg)-2*K*dg
        ddz=(U*z-drift*dz+forcing)/f
        return [dg,ddg,dz,ddz]
    sol=solve_ivp(rhs,(2+width,max(uppers)),[h,dh,H,dH],t_eval=uppers,
                  method='DOP853',rtol=rtol,atol=rtol/100)
    if not sol.success:
        raise RuntimeError(sol.message)
    r=sol.t;g,dg,z,dz=sol.y
    J=np.exp(-2j*w*tortoise(r))*(1-2/r)*(g*dz-dg*z)
    u=Series(0,series_coefficients(2,0,w,20))
    Wp=(u*(2j*w*(u*inverse_lapse(20))+u.derivative())).truncate(20)
    W=integrand_series(2,0,w,20)+(-2*K)*Wp
    V=antiderivative_series(W,w,24,20)
    outgoing=C*C*np.exp(2j*w*tortoise(r))*V.evaluate(r)
    residual=J-outgoing
    expected=-2*shift*P
    return dict(K=[K.real,K.imag],expected=[expected.real,expected.imag],
        rows=[dict(r=float(x),residual=[float(v.real),float(v.imag)],
                   absolute_error=float(abs(v-expected)),correction_ratio_abs=float(abs(zz/gg)))
              for x,v,zz,gg in zip(r,residual,z,g)])


if __name__=='__main__':
    for shift in (0j,.01+0j,.01j):
        print(json.dumps(run(shift)))
