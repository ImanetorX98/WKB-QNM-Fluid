#!/usr/bin/env python3
"""Test della struttura secolare esatta di H/G.  Nessun adattamento."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vaidya_secular_series import secular_structure  # noqa: E402

MODES = ((2, 0), (3, 0), (4, 0), (2, 1), (2, 2), (3, 2))


class SecularSeriesTest(unittest.TestCase):
    """I coefficienti escono dall'inversione della serie, non da un fit."""

    def test_quadratic_coefficient_is_exactly_twice_i_omega(self):
        for ell, spin in MODES:
            with self.subTest(ell=ell, spin=spin):
                row = secular_structure(ell, spin)
                self.assertLess(abs(row["quadratic"] - 2j * row["omega"]), 1.0e-12)

    def test_leading_three_orders_are_mode_independent(self):
        """r^1 e ln r stanno a 2i*om come 8 e 24, per ogni multipolo e spin.

        E' il controllo che rende il secolare geometrico e non del modo: se
        dipendesse da ell o da s non sarebbe trasporto lungo il raggio.
        """
        for ell, spin in MODES:
            with self.subTest(ell=ell, spin=spin):
                ratios = secular_structure(ell, spin)["ratios"]
                self.assertAlmostEqual(ratios["linear"].real, 8.0, places=6)
                self.assertAlmostEqual(ratios["linear"].imag, 0.0, places=6)
                self.assertAlmostEqual(ratios["log"].real, 24.0, places=6)
                self.assertAlmostEqual(ratios["log"].imag, 0.0, places=6)

    def test_logarithm_does_not_vanish(self):
        """H/G non e' un polinomio: l'adattamento polinomiale era un artefatto."""
        row = secular_structure(2, 0)
        self.assertGreater(abs(row["log"]), 1.0)

    def test_stable_against_truncation(self):
        """I tre coefficienti di testa non dipendono dal troncamento."""
        base = secular_structure(2, 0, terms=26)
        wide = secular_structure(2, 0, terms=34)
        for key in ("quadratic", "linear", "log"):
            self.assertLess(abs(base[key] / wide[key] - 1.0), 1.0e-10)


if __name__ == "__main__":
    unittest.main()
