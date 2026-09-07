from __future__ import annotations

import unittest

from calculations.poschl_teller_madelung_benchmark import (
    benchmark_row,
    exact_frequency,
    madelung_indicator,
)


class PoschlTellerMadelungTests(unittest.TestCase):
    def test_exact_frequency(self) -> None:
        omega = exact_frequency(2.0, 1)
        self.assertAlmostEqual(omega.real, (4.0 - 0.25) ** 0.5)
        self.assertAlmostEqual(omega.imag, -1.5)

    def test_madelung_finite_difference(self) -> None:
        indicator, residual, nodal = madelung_indicator(4.0, 0, points=3001)
        self.assertFalse(nodal)
        self.assertGreater(indicator, 0.0)
        self.assertLess(residual, 2.0e-4)

    def test_odd_mode_is_flagged_as_nodal(self) -> None:
        indicator, residual, nodal = madelung_indicator(4.0, 1, points=3001)
        self.assertTrue(nodal)
        self.assertNotEqual(indicator, indicator)
        self.assertNotEqual(residual, residual)

    def test_wkb3_improves_reference_case(self) -> None:
        row = benchmark_row(4.0, 0)
        self.assertLess(row.wkb3_relative_error, row.wkb1_relative_error)


if __name__ == "__main__":
    unittest.main()
