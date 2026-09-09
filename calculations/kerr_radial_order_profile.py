#!/usr/bin/env python3
"""Il potenziale radiale efficace di Kerr ha un termine di ordine eps^1?

Struttura (unita' M=1, H=r^2+a^2, Delta=r^2-2r+a^2, h=sqrt(H)):

    q/L^2 = Q0(r) - eps*Delta*A1/H^2 - eps^2*[Delta*A2/H^2 + h''/h],

    Q0(r) = (Omega - mu*a/H)^2 - Delta*Abar0/H^2,
    Abar0 = A0 + chat^2 - 2*mu*chat,      chat = a*Omega,

dove A = A0 + A1/L + A2/L^2 + ... e' l'espansione eikonale dell'autovalore
sferoidale e lambda = A + a^2 omega^2 - 2 a m omega.  I due termini di
convenzione sono puramente O(L^2): entrano in Abar0, non in A1.

MISURA NON CIRCOLARE.  Definire A1 come coefficiente di 1/L e poi verificarne
lo scaling sarebbe tautologico.  Si procede invece cosi': a ogni ell si calcola
il potenziale **esatto** con l'autovalore sferoidale esatto A_{lm}(a*omega), e
lo si confronta con la sua forma eikonale di testa Q0.  Il residuo

    D(r) = q(r)/L^2 - Q0(r)

deve andare come eps^1 se esiste un termine di rotazione, e come eps^2 se non
esiste.  Il controllo a = 0 e' il caso di Schwarzschild, dove il criterio
prevede pendenza 2.

Uso: python3.13 kerr_radial_order_profile.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, fsolve

sys.path.insert(0, str(Path(__file__).resolve().parent))

from kerr_eikonal_order_test import spheroidal_eigenvalue  # noqa: E402


def horizon(spin: float) -> float:
    return 1.0 + np.sqrt(1.0 - spin**2)


def leading_potential(
    r: np.ndarray | float, spin: float, mu: float, omega_hat: float, abar0: float
) -> np.ndarray | float:
    h_squared = r**2 + spin**2
    delta = r**2 - 2.0 * r + spin**2
    return (omega_hat - mu * spin / h_squared) ** 2 - delta * abar0 / h_squared**2


def eikonal_solution(spin: float, mu: float, ell_reference: int = 400) -> dict[str, float]:
    """Risolve autoconsistentemente (Omega, r_peak, A0) nel limite eikonale.

    Condizione di QNM eikonale: Q0 ha una radice doppia al raggio dell'orbita
    fotonica, cioe' Q0(r0)=0 e dQ0/dr(r0)=0.  A0 dipende da chat=a*Omega, che
    dipende dalla soluzione: il sistema e' chiuso per iterazione.
    """

    def angular_leading(chat: float) -> float:
        large_l = ell_reference + 0.5
        m = int(round(mu * large_l))
        value = spheroidal_eigenvalue(ell_reference, m, chat * large_l)
        return float(np.real(value)) / large_l**2

    def residuals(unknowns: np.ndarray, current_spin: float) -> list[float]:
        omega_hat, r0 = unknowns
        # Il picco deve restare fuori dall'orizzonte: si penalizza l'uscita
        # invece di lasciare che il solutore scappi verso r grandi.
        r_horizon = horizon(current_spin)
        if not np.isfinite(r0) or r0 <= r_horizon + 1.0e-6 or r0 > 20.0:
            return [1.0e3, 1.0e3]
        chat = current_spin * omega_hat
        abar0 = angular_leading(chat) + chat**2 - 2.0 * mu * chat
        step = 1.0e-6
        q0 = leading_potential(r0, current_spin, mu, omega_hat, abar0)
        derivative = (
            leading_potential(r0 + step, current_spin, mu, omega_hat, abar0)
            - leading_potential(r0 - step, current_spin, mu, omega_hat, abar0)
        ) / (2.0 * step)
        return [q0, derivative]

    # Continuazione in a: si parte da Schwarzschild, dove la soluzione e' nota
    # in forma chiusa, e si sale a passi usando la soluzione precedente.
    current = np.array([1.0 / (3.0 * np.sqrt(3.0)), 3.0])
    steps = max(1, int(np.ceil(abs(spin) / 0.05)))
    for index in range(1, steps + 1):
        current_spin = spin * index / steps
        current = fsolve(residuals, current, args=(current_spin,), xtol=1.0e-12)
    omega_hat, r_peak = current
    chat = spin * omega_hat
    a0 = angular_leading(chat)
    return {
        "omega_hat": float(omega_hat),
        "r_peak": float(r_peak),
        "chat": float(chat),
        "A0": a0,
        "Abar0": a0 + chat**2 - 2.0 * mu * chat,
    }


def exact_scaled_potential(
    r: np.ndarray, ell: int, spin: float, mu: float, omega_hat: float
) -> np.ndarray:
    """q(r)/L^2 con l'autovalore sferoidale esatto a quel ell."""
    large_l = ell + 0.5
    m = int(round(mu * large_l))
    omega = large_l * omega_hat
    eigenvalue = float(np.real(spheroidal_eigenvalue(ell, m, spin * omega)))
    separation = eigenvalue + spin**2 * omega**2 - 2.0 * spin * m * omega

    h_squared = r**2 + spin**2
    delta = r**2 - 2.0 * r + spin**2
    amplitude = np.sqrt(h_squared)
    # h''/h con '' rispetto a r_*, dr_*/dr = H/Delta
    first = (delta / h_squared) * np.gradient(amplitude, r)
    second = (delta / h_squared) * np.gradient(first, r)
    geometric = second / amplitude

    q = (omega - m * spin / h_squared) ** 2 - delta * separation / h_squared**2 - geometric
    return q / large_l**2


def order_measurement(
    spin: float,
    mu: float,
    ells: tuple[int, ...] = (40, 60, 80, 120, 160),
    window: float = 0.6,
) -> dict[str, float | list]:
    solution = eikonal_solution(spin, mu)
    r_peak = solution["r_peak"]
    r = np.linspace(r_peak - window, r_peak + window, 1201)
    leading = leading_potential(r, spin, mu, solution["omega_hat"], solution["Abar0"])

    epsilons, deviations = [], []
    for ell in ells:
        exact = exact_scaled_potential(r, ell, spin, mu, solution["omega_hat"])
        interior = slice(6, -6)
        deviations.append(float(np.max(np.abs(exact[interior] - leading[interior]))))
        epsilons.append(1.0 / (ell + 0.5))
    slope = float(np.polyfit(np.log(epsilons), np.log(deviations), 1)[0])
    return {
        **solution,
        "epsilons": epsilons,
        "deviations": deviations,
        "slope": slope,
    }


def main() -> None:
    print("Pendenza in eps del residuo  D(r) = q/L^2 - Q0(r)")
    print("Pendenza 2 = nessun termine di rotazione; pendenza 1 = termine O(eps).")
    print()
    print("    a     mu   Omega_eik   r_picco    A0       pendenza")
    for spin, mu in (
        (0.0, 0.5),
        (0.0, 0.9),
        (0.3, 0.5),
        (0.6, 0.5),
        (0.9, 0.5),
        (0.6, 0.9),
        (0.9, 0.9),
    ):
        result = order_measurement(spin, mu)
        marker = "  <- controllo statico" if spin == 0.0 else ""
        print(
            f"  {spin:4.1f}  {mu:4.1f}   {result['omega_hat']:8.6f}  {result['r_peak']:8.5f}  "
            f"{result['A0']:7.5f}   {result['slope']:7.4f}{marker}"
        )
    print()
    print("Controllo di consistenza: per a=0 il limite eikonale deve dare")
    print(f"  Omega = 1/(3 sqrt 3) = {1.0/(3.0*np.sqrt(3.0)):.8f},  r_picco = 3")


if __name__ == "__main__":
    main()
