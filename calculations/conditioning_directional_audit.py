"""Directional WKB-map sensitivity, not spectral or inverse conditioning.

Central differences compare the same prediction map: no exact-reference floor.
Near the peak these are formal nonuniform WKB diagnostics, not validated physics.
"""
import json
import numpy as np
from scipy.special import lambertw
from vaidya_solvability import tortoise, leaver_qnm
from madelung_wkb_prediction import predict


def audit(ell=70, spin=2, points=40001, step=1e-5, window=(20.,50.)):
    om=leaver_qnm(ell,0,spin)
    L=ell+.5
    t=np.linspace(float(tortoise(2.3)),float(tortoise(100.)),points)
    r=2*(1+lambertw(np.exp(t/2-1)).real)
    V=(1-2/r)*(ell*(ell+1)/r**2+2*(1-spin**2)/r**3)
    def fun(w):
        return predict((w*w-V)/L**2,1/L,t[1]-t[0],2)
    base=fun(om)
    h=step*abs(om)
    da=(fun(om+h)-fun(om-h))/(2*h)
    db=(fun(om+1j*h)-fun(om-1j*h))/(2*h)
    mask=(r>window[0])&(r<window[1])
    def median(a):
        return float(np.median(a[mask]))
    return dict(ell=ell,spin=spin,points=points,step=step,window=window,
        radial=median(abs((om.real*da+om.imag*db)/base)),
        real_component=median(abs(om.real*da/base)),
        imag_component=median(abs(om.imag*db/base)),
        pointwise_worst_normwise=median(abs(om)*np.hypot(da,db)/abs(base)),
        infinity_worst=2*abs(om)/abs(om.imag))


if __name__=='__main__':
    for ell in (20,40,70,100):
        print(json.dumps(audit(ell)))
    for points in (20001,80001):
        for step in (1e-4,1e-5):
            print(json.dumps(audit(points=points,step=step)))
    print(json.dumps(audit(window=(2.5,4.))))
