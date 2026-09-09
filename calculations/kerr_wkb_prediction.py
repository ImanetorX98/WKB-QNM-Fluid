#!/usr/bin/env python3
"""Previsione analitica di Q_M per Kerr, e confronto con la numerica.

Il Teorema 2 del manoscritto (assenza di ordini dispari) assume q indipendente
da eps.  In Kerr non e' cosi': il potenziale riscalato ha una sua espansione

    qhat = q0 + eps*q1 + eps^2*q2,      q1 = -Delta*A1/H^2,

e rifacendo l'espansione del momento u = u0 + eps*u1 + eps^2*u2 si trova

    u0 = sqrt(q0),
    u1 = q1/(2 sqrt(q0))            <- ordine dispari, assente nel caso statico
    u2 = [q2 - u1^2 + Q[u0]]/(2 u0),  Q[u] = (u^-1/2)''/(u^-1/2).

Sostituendo nella chiusura P^2 + Q_M = Re qhat, i termini di ordine eps si
cancellano identicamente:  q1 - 2 u0 u1 = 0.  Resta la previsione

    Q_M = -eps^2 * (5 q0'^2 - 4 q0 q0'') / (16 q0^2) + O(eps^3),   ' = d/dr_*

che e' il funzionale di Bohm valutato sull'ampiezza WKB di testa A0 = q0^(-1/4).
E' una formula chiusa, senza parametri liberi, che si ricava dal solo potenziale
di testa: non richiede di integrare l'ODE.

Questo script la confronta con Q_M estratto dall'ampiezza integrata.  Il
confronto e' molto piu' stringente di un fit di pendenza: se la numerica devia,
dice **dove** e di quanto.

Uso: python3.13 kerr_wkb_prediction.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from kerr_madelung_profile import kerr_profile, tortoise_derivative  # noqa: E402
from kerr_radial_order_profile import eikonal_solution, leading_potential  # noqa: E402


def leading_madelung(
    r: np.ndarray, spin: float, mu: float, omega_hat: float, abar0: float
) -> np.ndarray:
    """-(5 q0'^2 - 4 q0 q0'')/(16 q0^2), con ' = d/dr_*.

    Il fattore eps^2 e' lasciato fuori: e' il coefficiente della previsione.
    """
    q0 = np.asarray(leading_potential(r, spin, mu, omega_hat, abar0), dtype=float)
    lapse = np.asarray(tortoise_derivative(r, spin), dtype=float)
    first = lapse * np.gradient(q0, r)
    second = lapse * np.gradient(first, r)
    return -(5.0 * first**2 - 4.0 * q0 * second) / (16.0 * q0**2)


def compare(
    ell: int, spin: float, mu: float, window: tuple[float, float] = (20.0, 50.0)
) -> dict[str, float]:
    data = kerr_profile(ell, spin, mu)
    solution = eikonal_solution(spin, mu)
    r = np.asarray(data["r"])
    epsilon = float(data["epsilon"])

    predicted = epsilon**2 * leading_madelung(
        r, spin, mu, solution["omega_hat"], solution["Abar0"]
    )
    measured = np.asarray(data["q_madelung"])

    mask = (r > window[0]) & (r < window[1])
    mask[:8] = False
    mask[-8:] = False
    relative = np.abs(measured[mask] - predicted[mask]) / np.abs(predicted[mask])
    return {
        "ell": float(ell),
        "epsilon": epsilon,
        "predicted": float(np.max(np.abs(predicted[mask]))),
        "measured": float(np.max(np.abs(measured[mask]))),
        "median_relative_error": float(np.median(relative)),
        "max_relative_error": float(np.max(relative)),
    }


def main() -> None:
    print("Previsione  Q_M = -eps^2 (5 q0'^2 - 4 q0 q0'')/(16 q0^2)")
    print("contro Q_M estratto dall'ampiezza integrata, finestra 20<r<50.")
    print()
    for spin, mu in ((0.0, 0.5), (0.3, 0.5), (0.6, 0.5), (0.9, 0.5), (0.9, 0.9)):
        print(f"=== a={spin}, mu={mu} ===")
        print("   ell     previsto        misurato      errore mediano   errore max")
        for ell in (20, 30, 45, 70, 100):
            try:
                row = compare(ell, spin, mu)
            except RuntimeError as error:
                print(f"  {ell:4d}   integrazione fallita: {error}")
                continue
            print(
                f"  {row['ell']:4.0f}   {row['predicted']:.6e}   {row['measured']:.6e}   "
                f"{row['median_relative_error']:12.2%}   {row['max_relative_error']:.2%}"
            )
        print()


if __name__ == "__main__":
    main()
