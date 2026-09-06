"""Test rapidi: python3.13 -m unittest -v test_leaver.py"""

import unittest

from leaver_qnm import continued_fraction, leaver_qnm
from schwarzschild_wkb import qnm_wkb


class LeaverReferenceTests(unittest.TestCase):
    """Valori tabulati standard di Schwarzschild (Leaver 1985; tabelle di Berti)."""

    REFERENCE = {
        (2, 2, 0): complex(0.3736716844, -0.0889623157),
        (2, 2, 1): complex(0.3467109834, -0.2739148843),
        (2, 3, 0): complex(0.5994432884, -0.0927030294),
    }

    def test_matches_reference_to_seven_digits(self) -> None:
        for (spin, ell, overtone), expected in self.REFERENCE.items():
            omega = leaver_qnm(ell, overtone, spin)
            self.assertAlmostEqual(omega.real, expected.real, places=7)
            self.assertAlmostEqual(omega.imag, expected.imag, places=7)

    def test_residual_is_machine_zero_at_the_root(self) -> None:
        for spin, ell, overtone in self.REFERENCE:
            omega = leaver_qnm(ell, overtone, spin)
            residual = continued_fraction(omega, ell, spin, overtone)
            self.assertLess(abs(residual), 1.0e-10)

    def test_depth_independence(self) -> None:
        # Criterio di Leaver: aumentare la profondita' non deve muovere la radice.
        shallow = leaver_qnm(2, 0, 2, depth=800)
        deep = leaver_qnm(2, 0, 2, depth=4000)
        self.assertLess(abs(shallow - deep), 1.0e-9)


class LeaverConsistencyTests(unittest.TestCase):
    def test_wkb3_is_close_for_large_ell(self) -> None:
        # La WKB di barriera e' accurata per n << ell: serve come guess e come
        # controllo incrociato, non come sostituto.
        for ell in (8, 16):
            for overtone in (0, 1):
                exact = leaver_qnm(ell, overtone, 2)
                approximate = qnm_wkb(ell, overtone, 2, 3).omega_M
                self.assertLess(abs(exact - approximate) / abs(exact), 1.0e-4)

    def test_damped_branch_all_spins(self) -> None:
        for spin in (0, 1, 2):
            omega = leaver_qnm(max(spin, 2), 0, spin)
            self.assertGreater(omega.real, 0.0)
            self.assertLess(omega.imag, 0.0)

    def test_overtone_damping_increases(self) -> None:
        damping = [abs(leaver_qnm(8, n, 2).imag) for n in range(4)]
        self.assertEqual(damping, sorted(damping))


if __name__ == "__main__":
    unittest.main()
