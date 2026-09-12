#!/usr/bin/env python3
"""N2 — l'audit del §9.6 rifatto su sequenze a m intero ESATTO.

MOTIVO.  La tabella del §9.6 usa mu = 1/2, dove m = round(mu L) arrotonda
**sempre** di 0.25 e mu_eff = m/L non e' mu.  E' lo stesso difetto che ha
prodotto il falso A_1 del §9.3.  La conclusione qualitativa del §9.6 — nulla di
patologico nell'ampiezza — non dipende da quelle cifre, ma le cifre si'.

LA CURA NON RICHIEDE CODICE NUOVO.  Con mu = 2/3 e ell congruo a 1 modulo 3 si ha
m = mu(2 ell + 1)/2 = (2 ell + 1)/3, **intero esatto**, e `int(round(...))` non
arrotonda nulla.  Anche il riferimento angolare di default, ell = 400, soddisfa
la condizione.  Serve quindi solo scegliere i parametri, e verificarlo invece di
assumerlo — cosa che questo script fa con un assert esplicito.

CHE COSA MISURA.  Lo scarto fra Q_M previsto analiticamente (§5 del manoscritto)
e Q_M estratto dalla soluzione, con due frequenze:
  * quella **eikonale**, cioe' il limite di testa;
  * quella **autoconsistente**, che risolve la condizione di barriera con
    l'autovalore sferoidale esatto a quel ell.
Le derivate del riferimento sono analitiche, non per differenze finite: e' la
correzione gia' introdotta in `kerr_madelung_analytic.py`.

Uso: python3.13 calculations/claude_kerr_amplitude_exact_m.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from kerr_madelung_analytic import analytic_profile  # noqa: E402
from kerr_madelung_profile import eikonal_frequency  # noqa: E402
from kerr_wkb3_selfconsistent import selfconsistent_frequency  # noqa: E402
from madelung_wkb_prediction import predict  # noqa: E402

MU_NUM, MU_DEN = 2, 3          # mu = 2/3
WINDOW = (20.0, 50.0)


def exact_order(ell: int) -> int:
    """m = mu (2 ell + 1)/2, verificato intero.  Nessun arrotondamento."""
    value = MU_NUM * (2 * ell + 1) / (2 * MU_DEN)
    assert abs(value - round(value)) < 1e-12, (
        f"ell={ell} non da' m intero con mu={MU_NUM}/{MU_DEN}: usare ell = 1 mod {MU_DEN}")
    return int(round(value))


def prediction_error(ell: int, spin: float, omega: complex) -> float:
    """Scarto mediano fra Q_M previsto e Q_M esatto, nella finestra."""
    data = analytic_profile(ell, spin, MU_NUM / MU_DEN, omega=omega)
    radius, grid = np.asarray(data["r"]), np.asarray(data["rstar"])
    reference = np.asarray(data["q_madelung"]).real
    scale = ell + 0.5
    epsilon = 1.0 / scale
    step = float(grid[1] - grid[0])

    m = exact_order(ell)
    from kerr_eikonal_order_test import spheroidal_eigenvalue
    eigen = complex(spheroidal_eigenvalue(ell, m, spin * omega))
    separation = eigen + spin**2 * omega**2 - 2.0 * spin * m * omega

    h_squared = radius**2 + spin**2
    delta = radius**2 - 2.0 * radius + spin**2
    q = ((omega - m * spin / h_squared) ** 2
         - delta * separation / h_squared**2) / scale**2
    predicted = predict(q, epsilon, step, 2)

    mask = (radius > WINDOW[0]) & (radius < WINDOW[1])
    return float(np.median(np.abs(predicted - reference)[mask]
                           / np.abs(reference[mask])))


def main() -> None:
    print("N2 — §9.6 su sequenze a m intero esatto\n")
    print(f"mu = {MU_NUM}/{MU_DEN};  ell scelti congrui a 1 modulo {MU_DEN},")
    print("cosi' che m sia intero e `int(round(...))` non arrotondi nulla.\n")

    print("  verifica preliminare: m e' esatto?")
    for ell in (40, 61, 79, 100):
        m = exact_order(ell)
        print(f"    ell={ell:4d}  L={ell+0.5:6.1f}  mu*L={MU_NUM*(2*ell+1)/(2*MU_DEN):8.1f}"
              f"  m={m:4d}  scarto={abs(MU_NUM*(2*ell+1)/(2*MU_DEN)-m):.1e}")

    print("\n   a    ell      om eikonale         om autoconsistente     "
          "spostamento   err.eikonale  err.autoconsistente")
    for spin in (0.0, 0.3, 0.6, 0.9):
        for ell in (40, 61, 100):
            scale = ell + 0.5
            eikonal, _ = eikonal_frequency(spin, MU_NUM / MU_DEN, 0, scale)
            try:
                row = selfconsistent_frequency(ell, spin, MU_NUM / MU_DEN)
            except Exception as error:  # noqa: BLE001
                print(f"  {spin:4.1f} {ell:4d}   FALLITO: {type(error).__name__}")
                continue
            if not row["converged"]:
                print(f"  {spin:4.1f} {ell:4d}   non convergente: {row['message'][:40]}")
                continue
            consistent = row["omega"]
            try:
                crude = prediction_error(ell, spin, eikonal)
                fine = prediction_error(ell, spin, consistent)
            except Exception as error:  # noqa: BLE001
                print(f"  {spin:4.1f} {ell:4d}   profilo fallito: {type(error).__name__}")
                continue
            print(f"  {spin:4.1f} {ell:4d}   {eikonal.real:8.4f}{eikonal.imag:+8.4f}i   "
                  f"{consistent.real:8.4f}{consistent.imag:+8.4f}i   "
                  f"{row['shift']:.2e}     {crude:.2e}      {fine:.2e}")


if __name__ == "__main__":
    main()
