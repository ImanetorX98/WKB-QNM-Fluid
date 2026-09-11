import unittest
from fractions import Fraction
from kerr_fixed_mu_radial import modes, measure


class FixedMuRadial(unittest.TestCase):
    def test_integer_modes_really_hold_mu(self):
        for mu in (Fraction(2,5),Fraction(2,3),Fraction(4,5)):
            for ell,m in modes(mu,(40,60,80,120,160,240,300)):
                self.assertEqual(Fraction(2*m,2*ell+1),mu)

    def test_radial_scaling_is_second_order(self):
        for mu in (Fraction(2,5),Fraction(2,3),Fraction(4,5)):
            for a in (0.,.3,.6,.9):
                with self.subTest(a=a,mu=str(mu)):
                    out=measure(a,mu)
                    self.assertLess(abs(out['slope']-2),1e-3)
                    self.assertLess(abs(out['A1_fit']),1e-6)

    def test_reference_extrapolation_is_independent(self):
        a=measure(references=(180,240,300,360))
        b=measure(references=(300,360,420,480))
        self.assertLess(abs(a['slope']-b['slope']),1e-5)

    def test_at_selfconsistent_peak(self):
        out=measure(a=.6,omega_hat=None)
        self.assertLess(out['root_residual'],1e-9)
        self.assertLess(abs(out['slope']-2),1e-3)


if __name__=='__main__':
    unittest.main()
