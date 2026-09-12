import unittest
import numpy as np
from angular_chebyshev_independent import collocation_matrix,eigenvalue,eigenvalue_mp,independent_fit
from kerr_eikonal_order_test import spheroidal_eigenvalue


class IndependentAngular(unittest.TestCase):
    def test_spherical_limit(self):
        for ell,m in ((5,2),(10,3),(20,10)):
            self.assertLess(abs(eigenvalue(ell,m,0,60)-ell*(ell+1)),1e-6)

    def test_high_precision_convergence_and_other_solver(self):
        c=.6*40.5
        a,_=eigenvalue_mp(40,27,c,52)
        b,res=eigenvalue_mp(40,27,c,64)
        ref=spheroidal_eigenvalue(40,27,c)
        self.assertLess(abs(a/b-1),1e-11)
        self.assertLess(abs(b/ref-1),1e-11)
        self.assertLess(res,1e-24)

    def test_complex_mode(self):
        c=(.3+.2j)*40.5
        a,_=eigenvalue_mp(40,27,c,52)
        ref=spheroidal_eigenvalue(40,27,c)
        self.assertLess(abs(a/ref-1),1e-10)

    def test_linear_coefficient_independently_small(self):
        fit=independent_fit()
        self.assertLess(abs(fit['A1']),1e-6)


if __name__=='__main__':
    unittest.main()
