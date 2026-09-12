import unittest
import numpy as np
from compact_barrier_core_test import bump,solve_mode


class CoreScopeTest(unittest.TestCase):
    def test_deformation_vanishes_near_peak(self):
        x=np.linspace(-1,1,101)
        np.testing.assert_array_equal(bump(x,2.,.4),np.zeros_like(x))

    def test_on_shell_response_not_fixed_by_peak(self):
        w0,_,r0,e0=solve_mode(0.)
        w1,_,r1,e1=solve_mode(.01,w0)
        self.assertLess(max(r0,r1),1e-8)
        self.assertLess(w0.imag,0)
        self.assertLess(w1.imag,0)
        self.assertGreater(abs(w1-w0),1e-3)
        self.assertGreater(abs(e1-e0),1e-5)


if __name__=='__main__':
    unittest.main()
