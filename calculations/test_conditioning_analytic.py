import unittest
import numpy as np
import mpmath as mp
from conditioning_analytic import profile, summary
from conditioning_directional_audit import audit


def independent_mp(radius,omega,ell,spin):
    """Nested high-precision differentiation; no SymPy expressions reused."""
    L=mp.mpf(ell)+mp.mpf('.5')
    r=mp.mpf(radius)
    def q(x):
        f=1-2/x
        return (omega**2-f*(ell*(ell+1)/x**2+2*(1-spin**2)/x**3))/L**2
    def u(x):
        f=1-2/x
        qx=f*mp.diff(q,x)
        qxx=f*mp.diff(lambda y:(1-2/y)*mp.diff(q,y),x)
        return mp.sqrt(q(x))+(5*qx*qx-4*q(x)*qxx)/(32*L**2*q(x)**mp.mpf('2.5'))
    a=u(r)
    ux=(1-2/r)*mp.diff(u,r)
    uxx=(1-2/r)*mp.diff(lambda x:(1-2/x)*mp.diff(u,x),r)
    B=-mp.re(ux/a)/2-L*mp.im(a)
    C=-mp.re(uxx/a-(ux/a)**2)/2-L*mp.im(ux)
    return -(C+B*B)/L**2


class AnalyticConditioning(unittest.TestCase):
    def test_independent_multiprecision(self):
        with mp.workdps(45):
            w=mp.mpc('13.56418','-0.096217')
            for r in (22,35,60):
                expected=[independent_mp(r,w,70,2),
                    mp.diff(lambda h:independent_mp(r,w+h,70,2),0),
                    mp.diff(lambda h:independent_mp(r,w+1j*h,70,2),0)]
                actual=profile(r,complex(w),70,2)
                for a,b in zip(actual,expected):
                    self.assertLess(abs(float(a)/float(b)-1),1e-9)

    def test_central_frequency_differences(self):
        w=13.56418-.096217j
        q,da,db=profile(35,w)
        h=1e-5
        for direction,derivative in ((1,da),(1j,db)):
            fd=(profile(35,w+h*direction)[0]-profile(35,w-h*direction)[0])/(2*h)
            self.assertLess(abs(fd/derivative-1),1e-6)

    def test_old_grid_agrees_with_analytic(self):
        a=summary()
        b=audit(points=20001,step=1e-4)
        self.assertLess(abs(a['pointwise_normwise']/b['pointwise_worst_normwise']-1),1e-5)
        self.assertLess(abs(a['radial']-b['radial']),1e-3)

    def test_far_field_limit(self):
        q,da,db=profile(1e6,13.5-.1j)
        self.assertLess(abs(q/(-.1**2/70.5**2)-1),1e-8)
        self.assertLess(abs(db/(-2*(-.1)/70.5**2)-1),1e-8)
        self.assertLess(abs(da/db),1e-8)

    def test_global_norm_convergence(self):
        a=summary(points=501)
        b=summary(points=3001)
        self.assertLess(abs(a['profile_L2_normwise']/b['profile_L2_normwise']-1),1e-5)
        self.assertTrue(280<b['profile_L2_normwise']<286)
        self.assertTrue(1.9<b['profile_L2_radial']<2.1)


if __name__=='__main__':
    unittest.main()
