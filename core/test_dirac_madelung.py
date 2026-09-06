"""Test rapidi: python3.13 -m unittest -v test_dirac_madelung.py"""

import unittest

import numpy as np

from dirac_madelung_profile import (
    dirac_potential,
    dirac_potential_peak,
    dirac_qnm_wkb,
    hierarchy_scaling,
    scaling_exponents,
    integrate_profile,
)


class DiracPotentialTests(unittest.TestCase):
    def test_partner_split_is_the_spin_term(self) -> None:
        # V_+ - V_- = 2 K h', cioe' esattamente il termine di spin connection.
        x = np.linspace(2.05, 20.0, 400)
        kappa_abs = 3.0
        difference = dirac_potential(x, kappa_abs, 1) - dirac_potential(x, kappa_abs, -1)
        f = 1.0 - 2.0 / x
        expected = 2.0 * kappa_abs * np.sqrt(f) * (3.0 - x) / x**3
        np.testing.assert_allclose(difference, expected, rtol=1e-12, atol=1e-14)

    def test_peak_outside_horizon_and_near_photon_sphere(self) -> None:
        # Per K grande la barriera e' dominata da K^2 h^2, che picca a x=3.
        self.assertGreater(dirac_potential_peak(1.0, 1), 2.0)
        self.assertAlmostEqual(dirac_potential_peak(200.0, 1), 3.0, delta=2.0e-2)


class DiracQNMTests(unittest.TestCase):
    def test_reproduces_cho_2003(self) -> None:
        reference = {1: complex(0.176, -0.100), 2: complex(0.379, -0.0965)}
        for kappa_abs, expected in reference.items():
            _, omega = dirac_qnm_wkb(float(kappa_abs), 0, 1, 3)
            self.assertAlmostEqual(omega.real, expected.real, delta=1.5e-3)
            self.assertAlmostEqual(omega.imag, expected.imag, delta=1.5e-3)

    def test_partners_are_isospectral(self) -> None:
        # I partner di Darboux condividono lo spettro; la WKB3 lo vede entro
        # il proprio errore di troncamento, non meglio.
        _, omega_plus = dirac_qnm_wkb(4.0, 0, 1, 3)
        _, omega_minus = dirac_qnm_wkb(4.0, 0, -1, 3)
        self.assertLess(abs(omega_plus - omega_minus), 5.0e-3)

    def test_damped_branch(self) -> None:
        _, omega = dirac_qnm_wkb(2.0, 0, 1, 3)
        self.assertGreater(omega.real, 0.0)
        self.assertLess(omega.imag, 0.0)


class MadelungClosureTests(unittest.TestCase):
    def test_closure_is_exact(self) -> None:
        data = integrate_profile(3.0, points=1200)
        self.assertLess(float(data["closure_residual"]), 1.0e-10)

    def test_finite_difference_agrees_with_algebraic_qm(self) -> None:
        # Controllo indipendente: Q_M ricostruito da A''/A numerico.
        data = integrate_profile(3.0, points=4000)
        scale = float(np.max(np.abs(np.asarray(data["q_madelung"]))))
        self.assertLess(float(data["fd_residual"]), 1.0e-3 * scale)


class HierarchyScalingTests(unittest.TestCase):
    def test_spin_term_is_exactly_first_order(self) -> None:
        # Il termine di spin e' analitico, tau*eps*h': la pendenza in eps e' 1
        # a precisione macchina.
        exponents = scaling_exponents(hierarchy_scaling(region="far"))
        self.assertAlmostEqual(exponents["spin_exponent"], 1.0, delta=1.0e-6)

    def test_madelung_is_second_order_away_from_turning_points(self) -> None:
        # Lontano dai turning point Q_M va come eps^2.  La pendenza misurata
        # resta poco sotto 2 per contaminazione residua della condizione
        # ingoing a x_min finito; il punto e' che sia nettamente > 1.
        exponents = scaling_exponents(hierarchy_scaling(region="far"))
        self.assertGreater(exponents["madelung_exponent"], 1.7)
        self.assertLess(exponents["madelung_exponent"], 2.2)

    def test_hierarchy_breaks_at_the_barrier_peak(self) -> None:
        # Al massimo di barriera i due turning point coalescono: |P| scende a
        # zero e Q_M smette di seguire eps^2.  E' il limite di validita'
        # dichiarato nel rapporto, qui misurato.
        peak_rows = hierarchy_scaling(region="peak")
        far_rows = hierarchy_scaling(region="far")
        self.assertLess(max(row["momentum_min"] for row in peak_rows), 1.0e-3)
        self.assertGreater(min(row["momentum_min"] for row in far_rows), 0.1)
        self.assertLess(
            scaling_exponents(peak_rows)["madelung_exponent"],
            scaling_exponents(far_rows)["madelung_exponent"] - 0.25,
        )


if __name__ == "__main__":
    unittest.main()
