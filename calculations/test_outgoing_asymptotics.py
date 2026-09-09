"""Test rapidi: python3.13 -m unittest -v test_outgoing_asymptotics.py"""

import sys
import unittest
from pathlib import Path

from outgoing_asymptotics import (
    log_derivative,
    log_derivative_omega,
    residual,
    series_coefficients,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from leaver_qnm import leaver_qnm  # noqa: E402


class SeriesTests(unittest.TestCase):
    def test_first_coefficient_closed_form(self) -> None:
        # a_1 = -l(l+1)/(2 i omega), dalla ricorsione con a_0=1, a_{-1}=0.
        omega = 0.5 - 0.1j
        for ell in (2, 4, 7):
            coefficients = series_coefficients(ell, 0, omega)
            expected = -ell * (ell + 1) / (2j * omega)
            self.assertAlmostEqual(abs(coefficients[1] - expected), 0.0, places=12)

    def test_residual_improves_with_terms(self) -> None:
        omega = leaver_qnm(2, 0, 0)
        values = [residual(2, 0, omega, 60.0, terms) for terms in (4, 8, 12)]
        self.assertLess(values[1], values[0])
        self.assertLess(values[2], 1.0e-9)

    def test_residual_improves_with_radius(self) -> None:
        omega = leaver_qnm(4, 0, 0)
        values = [residual(4, 0, omega, radius, 8) for radius in (30.0, 60.0, 120.0)]
        self.assertEqual(values, sorted(values, reverse=True))


class LogDerivativeTests(unittest.TestCase):
    def test_correct_branch(self) -> None:
        # Il punto dell'intero modulo: D_+ deve tendere a +i*omega.  Integrando
        # l'ODE dall'esterno verso l'interno si otteneva -i*omega, ramo sbagliato.
        for ell in (2, 4):
            omega = leaver_qnm(ell, 0, 0)
            value = log_derivative(ell, 0, omega, 120.0)
            self.assertLess(abs(value - 1j * omega) / abs(omega), 2.0e-3)
            self.assertGreater(value.imag * omega.real, 0.0)  # segno concorde

    def test_converges_to_i_omega(self) -> None:
        omega = leaver_qnm(2, 0, 0)
        errors = [
            abs(log_derivative(2, 0, omega, radius) - 1j * omega)
            for radius in (30.0, 60.0, 120.0)
        ]
        self.assertEqual(errors, sorted(errors, reverse=True))

    def test_omega_derivative_tends_to_i(self) -> None:
        # Senza code D_+ = i omega, quindi D'_+ = i: le code danno lo scarto.
        omega = leaver_qnm(2, 0, 0)
        value = log_derivative_omega(2, 0, omega, 120.0)
        self.assertLess(abs(value - 1j), 5.0e-3)


if __name__ == "__main__":
    unittest.main()
