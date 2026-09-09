"""Test rapidi: python3.13 -m unittest -v test_kerr_eikonal_order.py"""

import unittest

import numpy as np

from kerr_eikonal_order_test import eikonal_fit, spheroidal_eigenvalue


class SpheroidalEigenvalueTests(unittest.TestCase):
    def test_spherical_limit(self) -> None:
        for ell, m in ((10, 3), (40, 20), (80, 72)):
            value = spheroidal_eigenvalue(ell, m, 0.0)
            self.assertAlmostEqual(value.real, ell * (ell + 1), places=8)
            self.assertAlmostEqual(value.imag, 0.0, places=12)

    def test_truncation_converged(self) -> None:
        # La troncatura della matrice non deve muovere l'autovalore.
        arguments = (120, 36, (0.9 - 0.09j) * 120.5)
        shallow = spheroidal_eigenvalue(*arguments, pad=60)
        deep = spheroidal_eigenvalue(*arguments, pad=200)
        self.assertLess(abs(shallow - deep), 1.0e-6)

    def test_mode_tracking_is_smooth(self) -> None:
        # Il tracciamento per prossimita' di autovalore falliva qui: la
        # successione in ell saltava fra modi.  Deve essere monotona.
        values = []
        for ell in (40, 60, 80, 120, 160):
            large_l = ell + 0.5
            m = int(round(0.3 * large_l))
            values.append(
                (spheroidal_eigenvalue(ell, m, (0.9 - 0.09j) * large_l) / large_l**2).real
            )
        differences = np.diff(values)
        self.assertTrue(np.all(differences > 0), values)
        self.assertLess(float(np.max(np.abs(differences))), 1.0e-3)


class EikonalOrderTests(unittest.TestCase):
    def test_no_order_one_term_without_rotation(self) -> None:
        # Il caso sferico: Langer elimina l'ordine L, e resta -1/4 a ordine L^0.
        fit = eikonal_fit(0.5, 0.0 + 0.0j)
        self.assertLess(abs(fit["A1"]), 1.0e-10)
        self.assertAlmostEqual(fit["A2"].real, -0.25, places=6)

    def test_rotation_generates_an_order_one_term(self) -> None:
        for mu, chat in ((0.5, 0.4 - 0.04j), (0.8, 0.6 - 0.06j), (0.9, 0.75 - 0.05j)):
            fit = eikonal_fit(mu, chat)
            self.assertLess(fit["residual"], 1.0e-5, (mu, chat))
            self.assertGreater(abs(fit["A1"]), 1.0e-2, (mu, chat))

    def test_order_one_term_is_not_a_langer_artefact(self) -> None:
        # L -> L+delta manda A1 -> A1 - 2 delta A0.  Una costante unica potrebbe
        # assorbirlo solo se delta=A1/(2A0) fosse comune a tutti i modi.
        shifts = [
            eikonal_fit(mu, chat)["langer_shift"]
            for mu, chat in ((0.5, 0.4 - 0.04j), (0.8, 0.6 - 0.06j), (0.9, 0.75 - 0.05j))
        ]
        array = np.array(shifts)
        spread = float(np.max(np.abs(array - array.mean())))
        self.assertGreater(spread / abs(array.mean()), 0.5)

    def test_real_part_survives_the_physical_limit(self) -> None:
        # Per un overtone fissato Im(chat)=O(eps): la parte reale di A1 deve
        # restare finita e stabile quando Im(chat) -> 0.
        real_parts = [
            eikonal_fit(0.5, 0.4 + 1j * imaginary)["A1"].real
            for imaginary in (0.0, -0.02, -0.05)
        ]
        self.assertTrue(all(value < -1.0e-2 for value in real_parts), real_parts)
        self.assertLess(max(real_parts) - min(real_parts), 5.0e-4)


if __name__ == "__main__":
    unittest.main()
