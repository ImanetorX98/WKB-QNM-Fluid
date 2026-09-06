"""Test rapidi: python -m unittest -v test_wkb.py"""

import unittest

from schwarzschild_wkb import physical_units, potential_peak, qnm_wkb


class SchwarzschildWKBTests(unittest.TestCase):
    def test_peak_is_outside_horizon(self) -> None:
        self.assertGreater(potential_peak(2, 2), 2.0)

    def test_gravitational_fundamental_wkb3(self) -> None:
        # Valore Iyer--Will atteso circa 0.3732-0.0892i; tolleranza volutamente
        # piu' larga dell'errore WKB rispetto a Leaver.
        omega = qnm_wkb(2, 0, spin=2, order=3).omega_M
        self.assertAlmostEqual(omega.real, 0.3732, delta=8.0e-4)
        self.assertAlmostEqual(omega.imag, -0.0892, delta=8.0e-4)

    def test_damped_branch(self) -> None:
        omega = qnm_wkb(2, 0, spin=0, order=3).omega_M
        self.assertGreater(omega.real, 0.0)
        self.assertLess(omega.imag, 0.0)

    def test_mass_scaling(self) -> None:
        omega = qnm_wkb(2, 0, spin=2, order=3).omega_M
        f10, tau10 = physical_units(omega, 10.0)
        f20, tau20 = physical_units(omega, 20.0)
        self.assertAlmostEqual(f10 / f20, 2.0)
        self.assertAlmostEqual(tau20 / tau10, 2.0)


if __name__ == "__main__":
    unittest.main()
