#!/usr/bin/env python3
"""Il fluido di Madelung di un QNM e' aperto: continuita' con sorgente.

Per psi = A exp(iS) con A, S reali e omega^2 = E + i Gamma complessa, l'ODE

    psi'' + (omega^2 - V) psi = 0

si separa in

    A'' - A S'^2 + (E - V) A = 0,                        (Hamilton-Jacobi)
    2 A' S' + A S'' + Gamma A = 0   <=>   (rho v)' = -Gamma rho,

con rho = A^2 e v = S'.  Per un modo smorzato Im(omega)<0 e Re(omega)>0 si ha
Gamma = 2 Re(omega) Im(omega) < 0, quindi la sorgente -Gamma rho e' positiva:
il flusso non si conserva.  Nessuna densita' reale positiva puo' quindi
sostenere un'interpretazione idrodinamica chiusa di un QNM.

Lo script misura il residuo di (rho v)' + Gamma rho lungo il profilo.

Uso: python3.13 open_continuity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from scalar_eikonal_scaling import integrate_profile  # noqa: E402


def continuity_residual(ell: int, spin: int = 2, overtone: int = 0) -> dict[str, float]:
    data = integrate_profile(ell, overtone, spin)
    xstar = np.asarray(data["xstar"])
    step = xstar[1] - xstar[0]
    amplitude = np.asarray(data["amplitude"])
    epsilon = float(data["epsilon"])

    # Il momento restituito e' riscalato (P = eps S'); qui serve S' fisico.
    velocity = np.asarray(data["momentum"]) / epsilon
    density = amplitude**2
    omega = complex(data["omega"])
    gamma = float(np.imag(omega**2))

    flux = density * velocity
    flux_derivative = np.gradient(flux, step)
    interior = slice(4, -4)
    residual = flux_derivative[interior] + gamma * density[interior]
    scale = np.abs(gamma) * density[interior]

    # Confronto: quanto varrebbe il residuo se si imponesse la conservazione.
    closed_residual = flux_derivative[interior]

    return {
        "ell": float(ell),
        "gamma": gamma,
        "max_relative_open": float(np.max(np.abs(residual) / scale)),
        "median_relative_closed": float(np.median(np.abs(closed_residual) / scale)),
    }


def main() -> None:
    print("Continuita' con sorgente:  (rho v)' = -Gamma rho,  Gamma = Im(omega^2)")
    print()
    print(" ell     Gamma        max|res. aperta|/(|Gamma|rho)   mediana|(rho v)'|/(|Gamma|rho)")
    for ell in (2, 4, 8, 16):
        row = continuity_residual(ell)
        print(
            f"{row['ell']:4.0f}  {row['gamma']:+.6f}   {row['max_relative_open']:26.3e}   "
            f"{row['median_relative_closed']:28.3f}"
        )
    print()
    print("Colonna 3: la legge con sorgente e' soddisfatta a precisione numerica.")
    print("Colonna 4: il flusso rho*v e' lontano dall'essere conservato, di ordine 1")
    print("rispetto alla scala |Gamma|rho.  Il fluido e' aperto, non chiuso.")


if __name__ == "__main__":
    main()
