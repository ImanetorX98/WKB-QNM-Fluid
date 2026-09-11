"""Read-only diagnostic of rounding contamination in the Kerr A1 fit.

Interpolation is a diagnostic, not a physical noninteger-m eigenmode.
The rational-mu test uses integer m at every multipole, with no interpolation.
"""
import numpy as np
from kerr_eikonal_order_test import eikonal_fit, spheroidal_eigenvalue


def interpolated_fit(mu=.5,chat=.4,ells=(40,60,80,120,160,240)):
    L=np.array(ells)+.5
    values=[]
    for ell,scale in zip(ells,L):
        m0=int(np.floor(mu*scale))
        ms=np.arange(m0-1,m0+3)
        ys=[spheroidal_eigenvalue(ell,int(m),chat*scale)/scale**2 for m in ms]
        coeff=np.polynomial.polynomial.polyfit(ms-mu*scale,ys,3)
        values.append(coeff[0])
    return np.linalg.lstsq(np.array([L*0+1,1/L,1/L**2]).T,values,rcond=None)[0]


def run():
    original=eikonal_fit(.5,.4)
    fixed=interpolated_fit()
    h=.001
    derivative=(interpolated_fit(mu=.5+h)[0]-interpolated_fit(mu=.5-h)[0])/(2*h)
    rational=eikonal_fit(2/3,.4,ells=(40,61,82,121,160,241))
    return dict(original_A1=original['A1'],interpolated_A1=fixed[1],
                predicted_rounding_A1=-.25*derivative,
                integer_fixed_mu_A1=rational['A1'])


if __name__=='__main__':
    print(run())
