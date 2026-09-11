#!/usr/bin/env python3
"""L'autovalore sferoidale NON ha termine O(1/L) a mu e chat fissati.

Invalida il §9 del manoscritto.  Il termine misurato in `kerr_eikonal_order_test`
e' l'arrotondamento di m all'intero: m = round(mu L) fa variare mu_eff = m/L di
+-0.25/L, e siccome A_0 dipende da mu questo genera un falso A_1 il cui SEGNO
dipende dalla parita' di ell.

Qui si sceglie mu = 2a/b con b | (2 ell + 1), cosi' m e' intero esatto.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from kerr_eikonal_order_test import spheroidal_eigenvalue  # noqa: E402


def exact_sequence(numerator, denominator, chat, count=10, start=60):
    """A/L^2 su ell tali che m = mu(2 ell+1)/2 sia intero esatto."""
    scales, values = [], []
    for ell in range(start, 10 * start):
        if (2 * ell + 1) % denominator:
            continue
        order = numerator * (2 * ell + 1) // (2 * denominator)
        assert abs(numerator * (2 * ell + 1) / (2 * denominator) - order) < 1e-12
        scale = ell + 0.5
        scales.append(scale)
        values.append(spheroidal_eigenvalue(ell, order, chat * scale) / scale**2)
        if len(scales) == count:
            break
    return np.array(scales), np.array(values)


class EigenvalueExpansionTest(unittest.TestCase):
    CASES = ((2, 3, 0.6 + 0j), (2, 3, 0.3 + 0.2j), (2, 5, 0.6 + 0j), (4, 5, 0.8 + 0j))

    def test_second_order_coefficient_is_constant(self):
        """L^2 (A/L^2 - A_0) costante  <=>  A_1 = 0."""
        for num, den, chat in self.CASES:
            with self.subTest(mu=f"{num}/{den}", chat=chat):
                scales, values = exact_sequence(num, den, chat)
                leading = ((scales[-1] ** 2 * values[-1] - scales[-2] ** 2 * values[-2])
                           / (scales[-1] ** 2 - scales[-2] ** 2))
                second = (values - leading) * scales**2
                spread = np.ptp(second.real) / abs(second.real.mean())
                self.assertLess(spread, 1.0e-3)

    def test_linear_fit_gives_zero(self):
        """Adattando A_0 + A_1/L + A_2/L^2 + A_3/L^3, esce A_1 ~ 0."""
        scales, values = exact_sequence(2, 3, 0.6 + 0j)
        design = np.vstack([scales ** (-k) for k in range(4)]).T
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        self.assertLess(abs(coefficients[1]), 1.0e-5)

    def test_rounding_reproduces_the_spurious_term_with_parity_sign(self):
        """Guardia: il falso A_1 cambia segno con la parita' di ell.

        E' la firma che identifica l'artefatto.  Se un giorno sparisse, la
        diagnosi andrebbe rifatta.
        """
        def rounded_fit(ells, mu=0.5, chat=0.6 + 0j):
            scales = np.array([e + 0.5 for e in ells])
            values = np.array([spheroidal_eigenvalue(e, int(round(mu * (e + 0.5))),
                                                     chat * (e + 0.5)) / (e + 0.5) ** 2
                               for e in ells])
            design = np.vstack([np.ones_like(scales), 1 / scales, 1 / scales**2]).T
            return np.linalg.lstsq(design, values, rcond=None)[0][1]

        even = rounded_fit((40, 60, 80, 120, 160, 240))
        odd = rounded_fit((41, 61, 81, 121, 161, 241))
        self.assertGreater(abs(even), 0.01)
        self.assertLess(abs(even + odd), 1.0e-4)     # segni opposti, moduli uguali


if __name__ == "__main__":
    unittest.main()
