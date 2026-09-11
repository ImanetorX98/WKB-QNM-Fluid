import unittest
import numpy as np
from conditioning_directional_audit import audit
from madelung_wkb_prediction import predict


class DirectionalConditioning(unittest.TestCase):
    def test_constant_potential_asymptote(self):
        w=12.-.1j
        L=60.
        def fun(z):
            return predict(np.full(101,(z/L)**2,dtype=complex),1/L,.1,2)[50]
        h=1e-5
        q=fun(w)
        da=(fun(w+h)-fun(w-h))/(2*h)
        db=(fun(w+1j*h)-fun(w-1j*h))/(2*h)
        self.assertAlmostEqual(q,-(w.imag/L)**2,places=14)
        self.assertAlmostEqual((w.real*da+w.imag*db)/q,2.,places=6)
        self.assertAlmostEqual(abs(w)*np.hypot(da,db)/abs(q),
                               2*abs(w)/abs(w.imag),places=5)

    def test_schwarzschild_direction_and_resolution(self):
        a=audit(ell=70,points=20001,step=1e-4)
        b=audit(ell=70,points=40001,step=1e-4)
        self.assertTrue(1.9<a['radial']<2.1)
        self.assertTrue(280<a['pointwise_worst_normwise']<286)
        self.assertLess(abs(a['pointwise_worst_normwise']/
                            b['pointwise_worst_normwise']-1),1e-4)


if __name__=='__main__':
    unittest.main()
