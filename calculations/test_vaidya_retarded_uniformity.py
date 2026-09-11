#!/usr/bin/env python3
"""Test: il secolare r^2 della correzione adiabatica e' il tempo ritardato."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vaidya_retarded_uniformity import (  # noqa: E402
    outgoing_characteristic_is_exact,
    secular_coefficient_from_retardation,
    secular_coefficient_from_source,
)


class RetardedUniformityTest(unittest.TestCase):
    def test_outer_equation_is_solved_by_any_retarded_profile(self):
        """Nella regione esterna l'adiabaticita' attorno a u = v-2r_* e' esatta."""
        self.assertTrue(outgoing_characteristic_is_exact())

    def test_secular_coefficient_is_twice_i_omega(self):
        """Il coefficiente quadratico della correzione forzata vale 2 i omega."""
        omega = sp.Symbol("omega")
        self.assertEqual(
            sp.simplify(secular_coefficient_from_source() - 2 * sp.I * omega), 0)

    def test_two_derivations_agree(self):
        """Equazione forzata e sviluppo del ritardo danno lo STESSO coefficiente.

        E' il controllo che distingue una diagnosi da una coincidenza numerica:
        il secolare non e' soltanto dello stesso ordine del ritardo, ne ha anche
        il coefficiente.
        """
        omega1 = sp.Symbol("omega_1", positive=True)
        from_retard, _ = secular_coefficient_from_retardation()
        self.assertEqual(sp.simplify(from_retard - 2 * sp.I * omega1), 0)


if __name__ == "__main__":
    unittest.main()
