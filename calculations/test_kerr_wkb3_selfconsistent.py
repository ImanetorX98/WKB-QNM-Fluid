"""Test rapidi: python3.13 -m unittest -v test_kerr_wkb3_selfconsistent.py"""

import sys
import unittest
from pathlib import Path

import numpy as np

from kerr_wkb3_selfconsistent import selfconsistent_frequency

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from schwarzschild_wkb import qnm_wkb  # noqa: E402


class SchwarzschildLimitTests(unittest.TestCase):
    def test_reduces_to_scalar_wkb3(self) -> None:
        # Con a=0 il solutore autoconsistente deve riprodurre il WKB3 scalare
        # gia' implementato in core/, che e' testato contro Iyer-Will.
        for ell in (10, 20, 40):
            result = selfconsistent_frequency(ell, 0.0, 0.5)
            reference = qnm_wkb(ell, 0, 0, 3).omega_M
            self.assertLess(abs(result["omega"] - reference) / abs(reference), 1.0e-9)

    def test_no_mu_dependence_without_rotation(self) -> None:
        low = selfconsistent_frequency(20, 0.0, 0.3)["omega"]
        high = selfconsistent_frequency(20, 0.0, 0.9)["omega"]
        self.assertLess(abs(low - high) / abs(low), 1.0e-9)


class KerrConvergenceTests(unittest.TestCase):
    CASES = ((0.3, 0.5), (0.6, 0.5), (0.9, 0.5))

    def test_converges_with_small_residual(self) -> None:
        for spin, mu in self.CASES:
            for ell in (20, 45):
                result = selfconsistent_frequency(ell, spin, mu)
                self.assertTrue(result["converged"], (spin, mu, ell))
                self.assertLess(result["residual"], 1.0e-9, (spin, mu, ell))

    def test_shift_from_eikonal_is_order_epsilon(self) -> None:
        # Lo spostamento rispetto alla frequenza eikonale deve ridursi con ell:
        # e' la correzione subprincipale, non un errore del solutore.
        shifts = [selfconsistent_frequency(ell, 0.6, 0.5)["shift"] for ell in (20, 45, 100)]
        self.assertEqual(shifts, sorted(shifts, reverse=True))
        self.assertLess(shifts[-1], 2.0e-3)

    def test_damped_and_prograde(self) -> None:
        for spin, mu in self.CASES:
            omega = selfconsistent_frequency(45, spin, mu)["omega"]
            self.assertGreater(omega.real, 0.0)
            self.assertLess(omega.imag, 0.0)

    def test_rotation_raises_real_part(self) -> None:
        # Il frame dragging alza Re(omega) a ell fissato.
        values = [selfconsistent_frequency(45, a, 0.5)["omega"].real for a in (0.0, 0.3, 0.6, 0.9)]
        self.assertEqual(values, sorted(values))


if __name__ == "__main__":
    unittest.main()
