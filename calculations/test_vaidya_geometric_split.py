#!/usr/bin/env python3
"""Test della separazione geometria/resto per la correzione forzata di Vaidya."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vaidya_geometric_split import geometric_primitive, identity_is_exact  # noqa: E402


class GeometricSplitTest(unittest.TestCase):
    def test_identity_holds_for_arbitrary_deviation(self):
        """Esatta per rho arbitraria: nessuno sviluppo, quindi vale sulla barriera."""
        self.assertEqual(identity_is_exact(), 0)

    def test_closed_form_expansion_gives_the_universal_coefficients(self):
        """r^2/2 + 4r + 12 ln r riproduce i rapporti 8 e 24 dopo il fattore 4 i om."""
        radius = sp.Symbol("r", positive=True)
        _, expansion = geometric_primitive()
        terms = sp.expand(expansion)
        self.assertEqual(sp.simplify(terms.coeff(radius, 2)), sp.Rational(1, 2))
        self.assertEqual(sp.simplify(terms.coeff(radius, 1)), 4)
        self.assertEqual(sp.simplify(terms.coeff(sp.log(radius))), 12)


if __name__ == "__main__":
    unittest.main()
