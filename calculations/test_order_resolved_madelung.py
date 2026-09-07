from __future__ import annotations

import unittest

from calculations.order_resolved_madelung_benchmark import (
    benchmark_mode,
    integrated_remainders,
)


class OrderResolvedMadelungTests(unittest.TestCase):
    def test_remainders_are_finite_and_positive(self) -> None:
        values = integrated_remainders(2, 2, 0.3736716844 - 0.0889623157j, points=1001)
        self.assertEqual(len(values), 3)
        for value in values:
            self.assertGreater(value, 0.0)
            self.assertLess(value, 1.0)

    def test_reference_mode(self) -> None:
        row = benchmark_mode(2, 0, 2, depth=1200)
        self.assertLess(row.wkb3_relative_error, row.wkb1_relative_error)
        self.assertGreater(row.remainder_0, 0.0)


if __name__ == "__main__":
    unittest.main()
