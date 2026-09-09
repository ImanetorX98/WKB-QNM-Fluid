import unittest
from audit_vaidya_quadrature import audit
from audit_vaidya_error_budget import budget


class QuadratureAudit(unittest.TestCase):
    def test_same_grid_series_antiderivative(self):
        coarse, fine = audit(2000), audit(4000)
        self.assertGreater(coarse['trap_error']/fine['trap_error'], 3.8)
        self.assertLess(fine['simpson_error'], 1e-8)

    def test_ode_derivative_removes_gradient_error(self):
        row = audit(4000)
        self.assertLess(row['ode_error'], 2e-7)
        self.assertGreater(row['gradient_error']/row['ode_error'], 100)

    def test_endpoint_mismatch_dominates_fine_grid_quadrature(self):
        row = budget(100000)
        self.assertLess(row['trap_error'], 1e-7)
        self.assertGreater(row['wrong_endpoints_error'], 1e-4)
        self.assertLess(row['simpson_error'], 1e-10)
        self.assertLess(row['endpoint_linear_remainder'], 1e-6)

    def test_pointwise_error_bounds_integrated_profile_error(self):
        row = budget(24000)
        self.assertLessEqual(row['discrete_profile_error'],
                             row['discrete_profile_error_bound'])
        self.assertLess(row['discrete_profile_error_bound'], 2e-6)

    def test_frobenius_initial_data_remove_profile_plateau(self):
        plane = budget(24000, horizon_order=0)
        corrected = budget(24000, horizon_order=2)
        self.assertLess(corrected['corrected_total_error'], 1e-9)
        self.assertGreater(plane['corrected_total_error']/
                           corrected['corrected_total_error'], 100)


if __name__ == '__main__':
    unittest.main()
