from __future__ import annotations

import unittest

from calculations.uniform_madelung_defect import (
    benchmark_row,
    power_law_exponent,
    uniform_defect_indicator,
)


class UniformMadelungDefectTests(unittest.TestCase):
    def test_regular_mode(self) -> None:
        defect, full, nodal = uniform_defect_indicator(4.0, 0, points=201)
        self.assertFalse(nodal)
        self.assertGreater(defect, 0.0)
        self.assertGreater(full, 0.0)

    def test_nodal_mode(self) -> None:
        defect, full, nodal = uniform_defect_indicator(4.0, 1, points=201)
        self.assertTrue(nodal)
        self.assertNotEqual(defect, defect)
        self.assertNotEqual(full, full)

    def test_leading_defect_tracks_wkb1_scaling(self) -> None:
        for overtone in (0, 2):
            rows = [
                benchmark_row(scale, overtone)
                for scale in (4.0, 6.0, 8.0, 12.0, 16.0)
            ]
            defect_power = power_law_exponent(rows, "uniform_defect")
            wkb1_power = power_law_exponent(rows, "wkb1_relative_error")
            self.assertLess(abs(defect_power - wkb1_power), 0.08)


if __name__ == "__main__":
    unittest.main()
