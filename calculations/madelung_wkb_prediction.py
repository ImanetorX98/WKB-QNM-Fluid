#!/usr/bin/env python3
"""Previsione analitica di Q_M dalla WKB complessa, e sua sensibilita'.

TEORIA.  Il Teorema 2 del manoscritto (assenza di ordini dispari) assume q
indipendente da eps.  Quando invece il potenziale riscalato ha una sua
espansione, come in Kerr dove qhat = q0 + eps*q1 + eps^2*q2, il momento acquista
un termine dispari u1 = q1/(2 sqrt(q0)) e la parita' si rompe.  Conviene allora
non separare qhat in ordini: si usa la serie WKB standard sul qhat **complesso**
completo,

    u = sqrt(q) + eps^2 (5 q'^2 - 4 q q'')/(32 q^(5/2)) + O(eps^4),   ' = d/dx.

AMPIEZZA.  Con omega complessa l'ampiezza reale della decomposizione di Madelung
non e' |u|^(-1/2): la parte immaginaria di u produce un fattore esponenziale.
Da psi = u^(-1/2) exp(i \\int u / eps) segue

    ln A = -1/2 ln|u| - (1/eps) \\int Im(u),

    Q_M = -eps^2 [ (ln A)'' + ((ln A)')^2 ].

Omettere il secondo termine di ln A sbaglia Q_M di un fattore ~50: e' l'errore
che questo modulo esiste per non ripetere.

RAMO.  u e -u sono entrambe soluzioni: va scelto il ramo coerente con la
condizione al contorno (uscente a destra della barriera).  Applicare lo stesso
ramo su entrambi i lati produce un errore del 10% che non svanisce con eps.

VALIDAZIONE.  Su Poschl-Teller n=0, dove Q_M = -(eps^2/4)(1+sech^2 y) e' esatto,
l'errore relativo mediano su 1<y<2.5 e'

    L        u = sqrt(q)     u = sqrt(q) + eps^2 u2
      8       3.9e-03            6.9e-04
     50       1.0e-04            5.1e-06
    300       2.8e-06            1.7e-07

cioe' la previsione e' accurata e migliora di un ordine con il termine u2.

Uso: python3.13 madelung_wkb_prediction.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


def wkb_momentum(q: np.ndarray, epsilon: float, step: float, order: int = 2) -> np.ndarray:
    """u = sqrt(q) + eps^2 (5 q'^2 - 4 q q'')/(32 q^(5/2)), q complesso."""
    if order not in (0, 2):
        raise ValueError("order deve essere 0 oppure 2")
    momentum = np.sqrt(q)
    if order == 2:
        first = np.gradient(q, step)
        second = np.gradient(first, step)
        momentum = momentum + epsilon**2 * (
            5.0 * first**2 - 4.0 * q * second
        ) / (32.0 * q ** 2.5)
    return momentum


def madelung_from_momentum(
    momentum: np.ndarray, epsilon: float, step: float
) -> np.ndarray:
    """Q_M = -eps^2 A''/A con A = |u^(-1/2) exp(i int u/eps)|.

    Le derivate di ln A sono analitiche nei termini di u: non si integra e non
    si differenzia una somma cumulativa.
    """
    log_modulus = np.log(np.abs(momentum))
    imaginary = np.imag(momentum)
    first = -0.5 * np.gradient(log_modulus, step) - imaginary / epsilon
    second = (
        -0.5 * np.gradient(np.gradient(log_modulus, step), step)
        - np.gradient(imaginary, step) / epsilon
    )
    return -(epsilon**2) * (second + first**2)


def predict(q: np.ndarray, epsilon: float, step: float, order: int = 2) -> np.ndarray:
    return madelung_from_momentum(wkb_momentum(q, epsilon, step, order), epsilon, step)


def poschl_teller_validation(
    scales: tuple[float, ...] = (8.0, 20.0, 50.0, 120.0, 300.0),
    points: int = 40001,
) -> list[dict[str, float]]:
    from poschl_teller_madelung_benchmark import exact_frequency, exact_wavefunction

    rows = []
    for scale in scales:
        epsilon = 1.0 / scale
        y = np.linspace(-3.0, 3.0, points)
        step = float(y[1] - y[0])
        psi, _ = exact_wavefunction(y, scale, 0)
        amplitude = np.abs(psi)
        exact = -(epsilon**2) * np.gradient(np.gradient(amplitude, step), step) / amplitude

        omega = exact_frequency(scale, 0) / scale
        q = omega**2 - 1.0 / np.cosh(y) ** 2
        # ramo uscente: solo y>0
        window = (y > 1.0) & (y < 2.5)
        row = {"L": scale}
        for order in (0, 2):
            predicted = predict(q, epsilon, step, order)
            row[f"order{order}"] = float(
                np.median(np.abs(predicted - exact)[window] / np.abs(exact[window]))
            )
        rows.append(row)
    return rows


def main() -> None:
    print("Validazione su Poschl-Teller n=0 (Q_M esatto noto in forma chiusa)")
    print("Errore relativo mediano su 1<y<2.5, ramo uscente")
    print()
    print("     L     u = sqrt(q)    u = sqrt(q) + eps^2 u2")
    for row in poschl_teller_validation():
        print(f"  {row['L']:6.1f}    {row['order0']:.3e}      {row['order2']:.3e}")
    print()
    print("La previsione non richiede di integrare l'ODE: usa solo il potenziale")
    print("complesso e la frequenza.  E' pero' fortemente sensibile a quest'ultima:")
    print("vedi research/kerr_madelung_sensitivity_2026-09-09.md.")


if __name__ == "__main__":
    main()
