#!/usr/bin/env python3
"""Algebra di serie asintotiche in 1/r, con derivate esatte.

MOTIVO.  Il primo tentativo di regolarizzare il numeratore della condizione di
solvibilita' costruiva l'antiderivata asintotica differenziando **numericamente**
il profilo con `np.gradient`.  Ogni iterazione amplifica il rumore: dopo cinque
passi il valore divergeva di undici ordini di grandezza, mentre lo scarto fra
punti di raccordo diversi *migliorava* fino a 1e-5.  L'indipendenza da L era
spuria, prodotta da due numeri enormi dominati dallo stesso termine divergente.

La cura e' rappresentare ogni funzione come serie troncata in 1/r e derivare
esattamente sui coefficienti.

RAPPRESENTAZIONE.  Una serie e' la coppia (offset, coeffs) che denota

    sum_k coeffs[k] * r^{-(offset+k)}

L'offset puo' essere negativo: serve per moltiplicare per r.

REGOLA APPRESA.  L'indipendenza dal parametro di regolarizzazione non basta:
il **valore** deve essere stabile separatamente al variare del troncamento.
Sono due controlli, non uno, e questo modulo espone entrambi.
"""

from __future__ import annotations

import numpy as np


class Series:
    """Serie troncata sum_k c[k] r^{-(offset+k)}."""

    __slots__ = ("offset", "coeffs")

    def __init__(self, offset: int, coeffs: np.ndarray) -> None:
        self.offset = int(offset)
        self.coeffs = np.asarray(coeffs, dtype=complex)

    def __len__(self) -> int:
        return self.coeffs.size

    def truncate(self, terms: int) -> "Series":
        if self.coeffs.size >= terms:
            return Series(self.offset, self.coeffs[:terms])
        padded = np.zeros(terms, dtype=complex)
        padded[: self.coeffs.size] = self.coeffs
        return Series(self.offset, padded)

    def __add__(self, other: "Series") -> "Series":
        offset = min(self.offset, other.offset)
        size = max(
            self.offset + len(self) - offset, other.offset + len(other) - offset
        )
        out = np.zeros(size, dtype=complex)
        out[self.offset - offset : self.offset - offset + len(self)] += self.coeffs
        out[other.offset - offset : other.offset - offset + len(other)] += other.coeffs
        return Series(offset, out)

    def __mul__(self, other: "Series | complex") -> "Series":
        if not isinstance(other, Series):
            return Series(self.offset, self.coeffs * other)
        return Series(
            self.offset + other.offset, np.convolve(self.coeffs, other.coeffs)
        )

    __rmul__ = __mul__

    def derivative(self) -> "Series":
        """d/dr: r^{-p} -> -p r^{-p-1}."""
        powers = self.offset + np.arange(len(self))
        return Series(self.offset + 1, -powers * self.coeffs)

    def evaluate(self, radius: float | np.ndarray) -> complex | np.ndarray:
        powers = self.offset + np.arange(len(self))
        return np.sum(
            self.coeffs[:, None] * np.asarray(radius)[None, :] ** (-powers[:, None]),
            axis=0,
        )

    def optimal_terms(self, radius: float) -> int:
        """Indice del termine minimo: dove troncare una serie asintotica."""
        powers = self.offset + np.arange(len(self))
        magnitudes = np.abs(self.coeffs) * radius ** (-powers.astype(float))
        finite = np.isfinite(magnitudes)
        if not finite.any():
            return len(self)
        return int(np.argmin(np.where(finite, magnitudes, np.inf))) + 1


def unit(terms: int) -> Series:
    coeffs = np.zeros(terms, dtype=complex)
    coeffs[0] = 1.0
    return Series(0, coeffs)


def inverse_lapse(terms: int) -> Series:
    """1/f = 1/(1-2/r) = sum_j 2^j r^{-j}."""
    return Series(0, 2.0 ** np.arange(terms))


def lapse(terms: int) -> Series:
    """f = 1 - 2/r."""
    coeffs = np.zeros(terms, dtype=complex)
    coeffs[0] = 1.0
    if terms > 1:
        coeffs[1] = -2.0
    return Series(0, coeffs)


def integrand_series(
    ell: int, spin: int, omega: complex, terms: int = 24
) -> Series:
    """W tale che l'integrando del numeratore sia W(r) exp(2 i omega r_*).

    Con Z = exp(2 i omega r_*) u nella regione esterna,
        P  = r (2 i omega u/f + u'),
        S  = 2 exp(2 i omega r_*) [(2 i omega/f) P + P'],
        mu Z S = 2 u [(2 i omega/f) P + P'] exp(2 i omega r_*) = W exp(...).
    """
    from outgoing_asymptotics import series_coefficients

    u = Series(0, series_coefficients(ell, spin, omega, terms))
    inv_f = inverse_lapse(terms)
    k = 2j * omega

    inner = (k * (u * inv_f) + u.derivative()).truncate(terms)
    momentum = Series(inner.offset - 1, inner.coeffs)          # P = r * inner
    bracket = (k * (momentum * inv_f) + momentum.derivative()).truncate(terms)
    return (2.0 * (u * bracket)).truncate(terms)


def antiderivative_series(
    integrand: Series, omega: complex, iterations: int, terms: int
) -> Series:
    """V con (V e^{2 i omega r_*})' = W e^{2 i omega r_*}, cioe' V' + (2i om/f)V = W.

    Ricorsione V <- (f/2i omega)(W - V'), tutte derivate esatte.
    """
    k = 2j * omega
    f = lapse(terms)
    current = Series(0, np.zeros(terms, dtype=complex))
    for _ in range(iterations):
        current = ((1.0 / k) * (f * (integrand + (-1.0) * current.derivative()))).truncate(
            terms
        )
    return current
