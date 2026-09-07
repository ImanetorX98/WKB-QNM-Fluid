from __future__ import annotations

import unittest

import numpy as np

from calculations.static_madelung_benchmark import (
    barrier_profile,
    benchmark_mode,
    inverse_tortoise,
    tortoise,
)


class StaticMadelungBenchmarkTests(unittest.TestCase):
    def test_tortoise_inverse(self) -> None:
        radii = np.array([2.01, 2.5, 3.0, 8.0, 30.0])
        reconstructed = inverse_tortoise(tortoise(radii))
        np.testing.assert_allclose(reconstructed, radii, rtol=2.0e-14, atol=2.0e-14)

    def test_barrier_profile_has_converged_curvature(self) -> None:
        omega = 0.3736716844 - 0.0889623157j
        profile = barrier_profile(2, 0, 2, omega, points=2501)
        self.assertGreater(float(profile["indicator"]), 0.0)
        self.assertLess(float(profile["fd_relative_residual"]), 3.0e-5)

    def test_wkb3_improves_fundamental_gravitational_mode(self) -> None:
        row = benchmark_mode(2, 0, 2, depth=1200)
        self.assertLess(row.wkb3_relative_error, row.wkb1_relative_error)
        self.assertLess(row.leaver_depth_shift, 1.0e-10)


if __name__ == "__main__":
    unittest.main()
