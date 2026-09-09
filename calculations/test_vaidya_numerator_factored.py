#!/usr/bin/env python3
"""Test del numeratore fattorizzato di Vaidya.

Da eseguire dalla radice del progetto:
    python3.13 -m unittest calculations.test_vaidya_numerator_factored
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from vaidya_numerator_factored import (  # noqa: E402
    Asymptotic, horizon_integral, horizon_series, numerator, outer_solve,
)
from leaver_qnm import leaver_qnm  # noqa: E402
from vaidya_solvability import tortoise  # noqa: E402

REFERENCE = 20.666545 - 40.326537j   # ell=2, s=0, n=0, regolarizzato ai due bordi


class NumeratorTest(unittest.TestCase):
    """Le proprieta' che rendono N una quantita' e non un artefatto."""

    @classmethod
    def setUpClass(cls):
        cls.base = numerator(uppers=(40.0, 50.0, 60.0, 70.0, 80.0))

    def test_independent_of_outer_cutoff(self):
        """L'indipendenza da L_+ e' il senso stesso della regolarizzazione."""
        self.assertLess(self.base["spread"], 1.0e-5)

    def test_reproduces_reference(self):
        self.assertAlmostEqual(abs(self.base["N"] / REFERENCE - 1.0), 0.0, delta=1.0e-5)

    def test_projection_does_not_vanish(self):
        """N != 0: la forzatura risonante ha componente sul modo."""
        self.assertGreater(abs(self.base["N"]), 1.0)

    def test_independent_of_horizon_width(self):
        """Il taglio interno deve cadere: sopra `width` numerico, sotto esatto."""
        values = [numerator(width=w, uppers=(40.0, 60.0, 80.0))["N"]
                  for w in (0.2, 0.5, 1.0, 1.25)]
        centre = np.mean(values)
        self.assertLess(max(abs(v - centre) for v in values) / abs(centre), 1.0e-6)

    def test_independent_of_match_point(self):
        """Il raccordo dentro/fuori non deve lasciare traccia."""
        values = [numerator(r_match=r, uppers=(40.0, 60.0, 80.0))["N"]
                  for r in (15.0, 25.0, 35.0)]
        centre = np.mean(values)
        self.assertLess(max(abs(v - centre) for v in values) / abs(centre), 1.0e-4)

    def test_works_across_spin_and_multipole(self):
        for ell, spin in ((3, 0), (2, 2), (2, 1)):
            with self.subTest(ell=ell, spin=spin):
                row = numerator(ell=ell, spin=spin, uppers=(40.0, 60.0, 80.0))
                self.assertLess(row["spread"], 1.0e-4)
                self.assertGreater(abs(row["N"]), 1.0)


class FactorisationTest(unittest.TestCase):
    """I due cambi di variabile che tolgono il range dinamico."""

    def test_frobenius_satisfies_indicial_relation(self):
        """h_1 = U(2)/(1/2 - 2 i om): la relazione indiciale al bordo interno."""
        omega = leaver_qnm(2, 0, 0)
        coefficients = horizon_series(2, 0, omega, 10)
        expected = (2 * 3 / 4.0 + 2.0 / 8.0) / (0.5 - 2j * omega)
        self.assertAlmostEqual(abs(coefficients[1] / expected - 1.0), 0.0, delta=1.0e-14)

    def test_horizon_exponent_is_integrable_only_for_the_fundamental(self):
        """Re(-4 i om) > -1 solo a n=0: la' il bordo non chiede regolarizzazione."""
        self.assertGreater(4 * leaver_qnm(2, 0, 0).imag, -1.0)
        self.assertLess(4 * leaver_qnm(2, 1, 0).imag, -1.0)

    def test_horizon_integral_matches_independent_quadrature(self):
        """Forma chiusa contro quadratura in t = ln x, dove l'integrale converge."""
        omega = leaver_qnm(2, 0, 0)
        terms, width = 60, 0.5
        coefficients = horizon_series(2, 0, omega, terms)
        index = np.arange(terms)
        first = np.concatenate([index[1:] * coefficients[1:], [0.0]])
        second = np.concatenate([index[1:-1] * index[2:] * coefficients[2:], [0.0, 0.0]])
        t = np.linspace(np.log(1e-13), np.log(width), 400001)
        x = np.exp(t)
        powers = x[:, None] ** index[None, :]
        star = (2.0 + x) + 2.0 * np.log(x / 2.0)
        integrand = (np.exp(-2j * omega * star) * (powers @ coefficients)
                     * 2.0 * ((powers @ first) + (2.0 + x) * (powers @ second)))
        from scipy.integrate import simpson
        reference = simpson(integrand * x, x=t)
        closed = horizon_integral(2, 0, omega, width, terms)
        self.assertLess(abs(closed / reference - 1.0), 1.0e-7)

    def test_series_solves_the_equation(self):
        """Il difetto R/f^2 misura quanto u risolve gia' l'ODE: deve sparire."""
        omega = leaver_qnm(2, 0, 0)
        asy = Asymptotic(2, 0, omega, 20)
        self.assertLess(abs(asy.at(30.0)[6][0]), 1.0e-11)
        self.assertLess(abs(asy.at(60.0)[6][0]), 1.0e-17)
        # e cala di ordini con il raggio: e' una serie asintotica, non un caso
        self.assertLess(abs(asy.at(60.0)[6][0]) / abs(asy.at(30.0)[6][0]), 1.0e-5)

    def test_outer_variable_is_flat(self):
        """g -> C: se non fosse piatta, la fattorizzazione non servirebbe."""
        base = numerator(uppers=(80.0,))
        omega, asy = base["omega"], Asymptotic(2, 0, base["omega"], 20)
        self.assertLess(base["flatness"], 1.0e-6)
        del omega, asy

    def test_residual_ingoing_mode_decays_at_the_ingoing_rate(self):
        """Guardia sulla diagnosi: g - C e' il modo entrante, non rumore.

        Con la BC troncata (un solo termine) l'ammissione entrante e' visibile e
        decade come |e^{-2 i om r_*}|.  Se un giorno decadesse a un altro tasso,
        la diagnosi del bordo sarebbe da rifare.
        """
        omega = leaver_qnm(2, 0, 0)
        asy = Asymptotic(2, 0, omega, 20)
        match = 25.0
        u, _, _, f, _, w, _, _ = asy.at(match)
        phase = np.exp(-2j * omega * tortoise(match))
        naive, naive_slope = 1.0, 0.0        # BC troncata: h=1, h'=0
        g0 = phase * naive / u[0]
        dg0 = phase * ((-1j * omega * naive / f[0] + naive_slope)
                       - naive * (1j * omega / f[0] + w[0])) / u[0]
        radius, g, _ = outer_solve(asy, match, 95.0, g0, dg0, 40001)
        constant = g[np.searchsorted(radius, 85.0)]

        def deviation(target):
            return abs(g[np.searchsorted(radius, target)] / constant - 1.0)

        predicted = abs(np.exp(-2j * omega * (tortoise(70.0) - tortoise(40.0))))
        self.assertAlmostEqual(deviation(70.0) / deviation(40.0) / predicted,
                               1.0, delta=0.05)


if __name__ == "__main__":
    unittest.main()
