"""Test rapidi: python3.13 -m unittest -v test_kerr_radial_order.py"""

import unittest

import numpy as np

from kerr_radial_order_profile import (
    eikonal_solution,
    horizon,
    leading_potential,
    order_measurement,
)


class EikonalSolutionTests(unittest.TestCase):
    def test_schwarzschild_limit(self) -> None:
        # A0 e' stimato a ell_reference=400, dove vale 1-1/(4L^2)=1-1.6e-6:
        # l'errore residuo su Omega e' ~1.5e-7 ed e' quello, non il solutore.
        solution = eikonal_solution(0.0, 0.5)
        self.assertAlmostEqual(solution["omega_hat"], 1.0 / (3.0 * np.sqrt(3.0)), places=6)
        self.assertAlmostEqual(solution["r_peak"], 3.0, places=5)
        self.assertAlmostEqual(solution["A0"], 1.0, places=5)

    def test_no_mu_dependence_without_rotation(self) -> None:
        low = eikonal_solution(0.0, 0.3)
        high = eikonal_solution(0.0, 0.9)
        self.assertAlmostEqual(low["omega_hat"], high["omega_hat"], places=9)
        self.assertAlmostEqual(low["r_peak"], high["r_peak"], places=8)

    def test_double_root_condition(self) -> None:
        for spin, mu in ((0.3, 0.5), (0.6, 0.9), (0.9, 0.5)):
            solution = eikonal_solution(spin, mu)
            r0 = solution["r_peak"]
            arguments = (spin, mu, solution["omega_hat"], solution["Abar0"])
            step = 1.0e-6
            derivative = (
                leading_potential(r0 + step, *arguments)
                - leading_potential(r0 - step, *arguments)
            ) / (2.0 * step)
            self.assertLess(abs(leading_potential(r0, *arguments)), 1.0e-10)
            self.assertLess(abs(derivative), 1.0e-8)

    def test_peak_outside_horizon_and_above_equatorial_orbit(self) -> None:
        # Le orbite fotoniche sferiche a mu<1 stanno sopra quella equatoriale
        # prograda, e si avvicinano al crescere di mu.
        for spin in (0.3, 0.6, 0.9):
            equatorial = 2.0 * (1.0 + np.cos(2.0 / 3.0 * np.arccos(-spin)))
            low = eikonal_solution(spin, 0.5)["r_peak"]
            high = eikonal_solution(spin, 0.9)["r_peak"]
            self.assertGreater(low, horizon(spin))
            self.assertGreater(low, equatorial)
            self.assertGreater(high, equatorial)
            self.assertLess(high, low)


class OrderSlopeTests(unittest.TestCase):
    def test_static_control_has_slope_two(self) -> None:
        result = order_measurement(0.0, 0.5, ells=(40, 60, 80, 120))
        self.assertGreater(result["slope"], 1.8)
        self.assertLess(result["slope"], 2.2)

    def test_rotation_gives_slope_one(self) -> None:
        for spin, mu in ((0.3, 0.5), (0.6, 0.5), (0.9, 0.9)):
            result = order_measurement(spin, mu, ells=(40, 60, 80, 120))
            self.assertGreater(result["slope"], 0.85, (spin, mu))
            self.assertLess(result["slope"], 1.15, (spin, mu))


if __name__ == "__main__":
    unittest.main()
