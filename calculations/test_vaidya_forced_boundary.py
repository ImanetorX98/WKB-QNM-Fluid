import unittest
from vaidya_forced_boundary import run


class ForcedBoundary(unittest.TestCase):
    def test_candidate_cancels_regularized_boundary(self):
        result=run(rtol=1e-12,uppers=(30.,35.,40.))
        for row in result['rows']:
            self.assertLess(row['absolute_error'],1e-5)

    def test_detuning_produces_predicted_complex_mismatch(self):
        for shift in (.01,.01j):
            result=run(shift=shift,rtol=1e-12,uppers=(30.,35.,40.))
            self.assertGreater(abs(complex(*result['expected'])),.09)
            for row in result['rows']:
                self.assertLess(row['absolute_error'],1e-5)

    def test_tighter_tolerance_reduces_error(self):
        coarse=run(rtol=1e-10,uppers=(40.,))['rows'][0]['absolute_error']
        fine=run(rtol=1e-12,uppers=(40.,))['rows'][0]['absolute_error']
        self.assertLess(fine,coarse/10)


if __name__=='__main__':
    unittest.main()
