#!/usr/bin/env python3
"""Serie asintotica uscente a grande r, e derivata logaritmica D_+(omega).

PERCHE'.  La norma generalizzata di Leung et al. (1998), eq. (2.16), richiede
D_{+0}(omega), la derivata logaritmica della soluzione **uscente** al punto di
raccordo, e la sua derivata in omega.  Non si ottiene integrando l'ODE
dall'esterno verso l'interno: con Im(omega)<0 la soluzione uscente cresce con
r_*, quindi in quella direzione e' il ramo **recessivo** e la contaminazione dal
ramo entrante lo sovrasta.  Un tentativo diretto restituisce D_+ = -i omega
invece di +i omega: segno opposto, ramo sbagliato.

La via corretta e' la serie asintotica.  Con R = exp(i omega r_*) u(r) e
M=1, l'equazione radiale diventa

    (1 - 2/r) u'' + (2/r^2 + 2 i omega) u' - [L/r^2 + 2(1-s^2)/r^3] u = 0,

e ponendo u = sum_k a_k r^{-k}, a_0 = 1, si ottiene

    a_k = { [k(k-1) - L] a_{k-1} - [2k(k-2) + 2(1-s^2)] a_{k-2} } / (2 i omega k).

La serie e' asintotica, non convergente: va troncata al termine minimo.

Uso: python3.13 outgoing_asymptotics.py
"""

from __future__ import annotations

import numpy as np


def series_coefficients(
    ell: int, spin: int, omega: complex, terms: int = 20
) -> np.ndarray:
    """a_0..a_{terms-1} della serie uscente."""
    angular = ell * (ell + 1)
    coefficients = np.zeros(terms, dtype=complex)
    coefficients[0] = 1.0
    for k in range(1, terms):
        previous = coefficients[k - 1]
        second = coefficients[k - 2] if k >= 2 else 0.0
        coefficients[k] = (
            (k * (k - 1) - angular) * previous
            - (2 * k * (k - 2) + 2 * (1 - spin**2)) * second
        ) / (2j * omega * k)
    return coefficients


def optimal_truncation(coefficients: np.ndarray, radius: float) -> int:
    """Indice del termine minimo: dove troncare una serie asintotica."""
    magnitudes = np.abs(coefficients) / radius ** np.arange(coefficients.size)
    # ignora a_0, che e' 1 per costruzione
    return int(np.argmin(magnitudes[1:])) + 1


def amplitude_and_derivative(
    ell: int, spin: int, omega: complex, radius: float, terms: int | None = None
) -> tuple[complex, complex]:
    """u(r) e du/dr dalla serie, troncata al termine minimo."""
    coefficients = series_coefficients(ell, spin, omega)
    cut = optimal_truncation(coefficients, radius) if terms is None else terms
    powers = radius ** -np.arange(cut + 1)
    value = float(1.0) * np.dot(coefficients[: cut + 1], powers)
    derivative = np.dot(
        -np.arange(cut + 1) * coefficients[: cut + 1], powers / radius
    )
    return complex(value), complex(derivative)


def log_derivative(ell: int, spin: int, omega: complex, radius: float) -> complex:
    """D_+ = (dR/dr_*)/R = i omega + f u'/u, dalla serie."""
    value, derivative = amplitude_and_derivative(ell, spin, omega, radius)
    f = 1.0 - 2.0 / radius
    return 1j * omega + f * derivative / value


def log_derivative_omega(
    ell: int, spin: int, omega: complex, radius: float, step: float | None = None
) -> complex:
    """D'_+ = dD_+/domega, per differenza centrata sulla serie.

    La serie e' una funzione esplicita di omega: nessuna integrazione, quindi
    nessuna contaminazione di ramo.
    """
    h = step if step is not None else 1.0e-6 * abs(omega)
    return (
        log_derivative(ell, spin, omega + h, radius)
        - log_derivative(ell, spin, omega - h, radius)
    ) / (2.0 * h)


def residual(ell: int, spin: int, omega: complex, radius: float, terms: int) -> float:
    """Residuo dell'ODE per u, come controllo della serie."""
    step = 1.0e-4 * radius
    values = [
        amplitude_and_derivative(ell, spin, omega, radius + offset, terms)[0]
        for offset in (-step, 0.0, step)
    ]
    first = (values[2] - values[0]) / (2 * step)
    second = (values[2] - 2 * values[1] + values[0]) / step**2
    f = 1.0 - 2.0 / radius
    angular = ell * (ell + 1)
    lhs = (
        f * second
        + (2.0 / radius**2 + 2j * omega) * first
        - (angular / radius**2 + 2 * (1 - spin**2) / radius**3) * values[1]
    )
    return abs(lhs) / abs(values[1])


def main() -> None:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))
    from leaver_qnm import leaver_qnm

    print("=== La serie soddisfa l'equazione? residuo relativo ===")
    print("  ell    r     termini=4     termini=8     termini=12")
    for ell in (2, 4):
        omega = leaver_qnm(ell, 0, 0)
        for radius in (30.0, 60.0, 120.0):
            row = [residual(ell, 0, omega, radius, t) for t in (4, 8, 12)]
            print(f"  {ell:3d} {radius:6.0f}   {row[0]:.3e}   {row[1]:.3e}   {row[2]:.3e}")
    print()

    print("=== D_+ dalla serie: deve tendere a +i*omega, non a -i*omega ===")
    print("  ell    r       D_+ (serie)             i*omega              D'_+")
    for ell in (2, 4):
        omega = leaver_qnm(ell, 0, 0)
        for radius in (30.0, 60.0, 120.0):
            d_plus = log_derivative(ell, 0, omega, radius)
            d_prime = log_derivative_omega(ell, 0, omega, radius)
            print(
                f"  {ell:3d} {radius:6.0f}  {d_plus.real:+.6f}{d_plus.imag:+.6f}i  "
                f"{(1j*omega).real:+.6f}{(1j*omega).imag:+.6f}i  "
                f"{d_prime.real:+.5f}{d_prime.imag:+.5f}i"
            )
        print()


if __name__ == "__main__":
    main()
