#!/usr/bin/env python3
"""Frequenza WKB3 di Kerr autoconsistente col potenziale effettivamente usato.

PERCHE'.  Il confronto fra previsione analitica di Q_M e numerica (Appendice C)
riesce su Schwarzschild e fallisce su Kerr.  La diagnosi in
`research/kerr_madelung_sensitivity_2026-09-09.md` e' che la colpa e' della
frequenza: usavamo quella **eikonale**, il cui errore per a!=0 e' ~1e-3, e Q_M
amplifica di ~1e2 l'errore su omega.  Serve una frequenza consistente col
potenziale a quel ell, non con il suo limite eikonale.

COSA SI RISOLVE.  A omega fissata l'autovalore sferoidale lambda e' un numero,
quindi

    q(r) = (omega - m a/H)^2 - Delta*lambda/H^2 - (1/h) d^2 h/dr_*^2

e' una funzione esplicita di r: le derivate rispetto a r_* si ottengono
simbolicamente iterando D_* = (Delta/H) d/dr, come gia' si fa per Schwarzschild.
La condizione di barriera di Iyer-Will, scritta in termini di q invece che di V,
e'

    q_0/sqrt(2 q_0'') = Lambda_2 - i(n+1/2)(1 + Lambda_3),

con Lambda_2, Lambda_3 costruiti dai rapporti q_4/q_2, (q_3/q_2)^2, ... che sono
invarianti sotto q <-> -V (cfr. Appendice A del manoscritto).

Il sistema e' autoconsistente perche' lambda dipende da omega: si risolve per
(omega, r_0) complessi simultaneamente, partendo dalla soluzione eikonale.

Uso: python3.13 kerr_wkb3_selfconsistent.py
"""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import fsolve

sys.path.insert(0, str(Path(__file__).resolve().parent))

from kerr_eikonal_order_test import spheroidal_eigenvalue  # noqa: E402
from kerr_madelung_profile import eikonal_frequency  # noqa: E402


@lru_cache(maxsize=None)
def _derivative_functions() -> tuple:
    """q e le sue derivate d^k q/dr_*^k, k=0..6, come funzioni di (r,a,m,lam,om)."""
    r, spin, azimuthal, separation, omega = sp.symbols("r a m lam om")
    h_squared = r**2 + spin**2
    delta = r**2 - 2 * r + spin**2
    amplitude = sp.sqrt(h_squared)

    def d_star(expression: sp.Expr) -> sp.Expr:
        return sp.together(delta / h_squared * sp.diff(expression, r))

    geometric = sp.together(d_star(d_star(amplitude)) / amplitude)
    q = (
        (omega - azimuthal * spin / h_squared) ** 2
        - delta * separation / h_squared**2
        - geometric
    )

    expressions = [sp.together(q)]
    for _ in range(6):
        expressions.append(d_star(expressions[-1]))
    symbols = (r, spin, azimuthal, separation, omega)
    return tuple(sp.lambdify(symbols, e, modules="numpy") for e in expressions)


def separation_constant(ell: int, m: int, spin: float, omega: complex) -> complex:
    """lambda = A_{lm}(a*omega) + a^2 omega^2 - 2 a m omega, con A complesso."""
    eigenvalue = complex(spheroidal_eigenvalue(ell, m, spin * omega))
    return eigenvalue + spin**2 * omega**2 - 2.0 * spin * m * omega


def iyer_will_residual(
    r_peak: complex, omega: complex, spin: float, m: int, separation: complex, overtone: int
) -> complex:
    """q_0/sqrt(2 q_0'') - Lambda_2 + i(n+1/2)(1+Lambda_3)."""
    functions = _derivative_functions()
    values = [f(r_peak, spin, m, separation, omega) for f in functions]
    q0, _, q2, q3, q4, q5, q6 = values
    alpha = overtone + 0.5
    root = np.sqrt(2.0 * q2 + 0j)

    lambda2 = (
        0.125 * (q4 / q2) * (0.25 + alpha**2)
        - (q3 / q2) ** 2 * (7.0 + 60.0 * alpha**2) / 288.0
    ) / root
    lambda3 = (
        5.0 * (q3 / q2) ** 4 * (77.0 + 188.0 * alpha**2) / 6912.0
        - (q3**2 * q4 / q2**3) * (51.0 + 100.0 * alpha**2) / 384.0
        + (q4 / q2) ** 2 * (67.0 + 68.0 * alpha**2) / 2304.0
        + (q3 * q5 / q2**2) * (19.0 + 28.0 * alpha**2) / 288.0
        - (q6 / q2) * (5.0 + 4.0 * alpha**2) / 288.0
    ) / (2.0 * q2)
    return q0 / root - lambda2 + 1j * alpha * (1.0 + lambda3)


def selfconsistent_frequency(
    ell: int,
    spin: float,
    mu: float,
    overtone: int = 0,
    tolerance: float = 1.0e-11,
) -> dict[str, complex | float]:
    """Risolve simultaneamente per (omega, r_0) complessi."""
    large_l = ell + 0.5
    m = int(round(mu * large_l))
    guess_omega, eikonal = eikonal_frequency(spin, mu, overtone, large_l)
    functions = _derivative_functions()

    def system(unknowns: np.ndarray) -> list[float]:
        omega = unknowns[0] + 1j * unknowns[1]
        r_peak = unknowns[2] + 1j * unknowns[3]
        separation = separation_constant(ell, m, spin, omega)
        # dq/dr_* = 0 al picco
        slope = functions[1](r_peak, spin, m, separation, omega)
        residual = iyer_will_residual(r_peak, omega, spin, m, separation, overtone)
        return [residual.real, residual.imag, slope.real, slope.imag]

    start = np.array(
        [guess_omega.real, guess_omega.imag, eikonal["r_peak"], 0.0], dtype=float
    )
    solution, info, status, message = fsolve(
        system, start, full_output=True, xtol=tolerance
    )
    omega = complex(solution[0], solution[1])
    r_peak = complex(solution[2], solution[3])
    return {
        "omega": omega,
        "r_peak": r_peak,
        "omega_eikonal": guess_omega,
        "shift": abs(omega - guess_omega) / abs(omega),
        "residual": float(np.max(np.abs(system(solution)))),
        "converged": status == 1,
        "message": message,
    }


def main() -> None:
    print("=== Controllo: a=0 deve riprodurre il WKB3 di Schwarzschild ===")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))
    from schwarzschild_wkb import qnm_wkb  # noqa: E402

    for ell in (10, 20, 40):
        result = selfconsistent_frequency(ell, 0.0, 0.5)
        reference = qnm_wkb(ell, 0, 0, 3).omega_M
        print(
            f"  ell={ell:3d}  autoconsistente {result['omega'].real:.9f}"
            f"{result['omega'].imag:+.9f}i   WKB3 scalare {reference.real:.9f}"
            f"{reference.imag:+.9f}i   scarto {abs(result['omega']-reference)/abs(reference):.2e}"
        )

    print()
    print("=== Kerr: quanto si sposta rispetto alla frequenza eikonale? ===")
    print("    a    mu   ell     omega autoconsistente        spostamento   residuo")
    for spin, mu in ((0.3, 0.5), (0.6, 0.5), (0.9, 0.5), (0.9, 0.9)):
        for ell in (20, 45, 100):
            result = selfconsistent_frequency(ell, spin, mu)
            flag = "" if result["converged"] else "  [NON CONVERGE]"
            print(
                f"  {spin:4.1f} {mu:4.1f} {ell:4d}   {result['omega'].real:12.7f}"
                f"{result['omega'].imag:+.7f}i   {result['shift']:.2e}   "
                f"{result['residual']:.1e}{flag}"
            )
        print()


if __name__ == "__main__":
    main()
