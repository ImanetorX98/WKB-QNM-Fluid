#!/usr/bin/env python3
"""Condizionamento di Q_M rispetto alla frequenza, su Schwarzschild.

--------------------------------------------------------------------------
VERDETTO: L'AMPLIFICAZIONE E' 2, NON 10^2.  IL §5 DEL MANOSCRITTO E' FALSO.
--------------------------------------------------------------------------
Misurata a ell = 20, 40, 70, 100, per s = 0 e s = 2, su tre finestre diverse:
l'amplificazione vale **2.0** ovunque, e tende a 2 esattamente al crescere della
perturbazione.  Non ~10^2, e non dipendente da nulla.

La ragione e' una regola della catena.  Nella regione esterna Q_M e' dominato da
-eps^2 (Im u / eps)^2 con Im u -> Im(omega)/L, cioe' Q_M ∝ (Im omega)^2.  Un
errore relativo delta su omega ne produce quindi **2 delta** su Q_M.  Fattore 2,
analitico, nessuna ostruzione.

Il §5 del manoscritto — promosso a risultato di testa il 9 settembre — afferma
un'amplificazione di due ordini di grandezza e ne deduce che nessun funzionale
costruito su Q_M possa predire la frequenza.  **La premessa non regge.**  Resta
in piedi il §6 (E_M = c|Lambda_3| allo 0.15%), che e' misurato
indipendentemente, e la sua ragione strutturale vera, che non e' il
condizionamento ma l'argomento di Riccati: un'espansione della Riccati non puo'
contenere piu' della WKB che essa e'.

--------------------------------------------------------------------------
TRE DIFETTI DELLA RICOSTRUZIONE, TUTTI TROVATI E CORRETTI
--------------------------------------------------------------------------
Il numero 10^2 non era riproducibile perche' il codice originale non e' nel
repository.  Ricostruendolo, tre difetti si sono nascosti l'uno dietro l'altro —
lo stesso schema del §11.

1. **Condizione al bordo dell'orizzonte.**  Integrare psi da r = 2 + delta con la
   sola BC entrante di testa semina il modo spurio.  Lo scarto a frequenza esatta
   scendeva di 160 volte stringendo delta da 1e-4 a 1e-10: contaminazione, non
   modello.  Qui la BC e' **Frobenius completa** (`horizon_values`).

2. **Griglia uniforme in r invece che in x_*.**  `predict` fa differenze finite
   assumendo passo costante, ma dx_* = dr/f varia di un fattore 5 fra r=2.5 e
   r=90.  Qui si integra **in x_*** con r come variabile di stato.

3. **Troncamento della serie di Frobenius.**  I coefficienti crescono con
   ell(ell+1): a ell=70, sessanta termini a x=0.5 lasciano un residuo che si
   propaga come modo spurio e fa **oscillare** Q_M di 5.9e-5 contro un valore di
   1.9e-6 — cambiando persino segno.  Con 120 termini l'ampiezza scende a 1.7e-8.
   Regola adottata: due termini per multipolo.

Solo dopo il terzo il pavimento **cala** con ell (1.7e-4, 4.4e-5, 1.5e-5, 7.1e-6
a ell = 20, 40, 70, 100) invece di crescere, ed e' a quel punto che
l'amplificazione si legge pulita.

--------------------------------------------------------------------------
COSTRUZIONE
--------------------------------------------------------------------------
Riferimento: Q_M dalla soluzione esatta alla frequenza di Leaver, con derivate
**analitiche**.  Con psi = e^{-i om r_*} h, h di ordine uno all'orizzonte,

    ln A     = Im(om) r_* + ln|h|,
    (ln A)'  = Im(om) + f Re(h'/h),
    (ln A)'' = f f_r Re(h'/h) + f^2 Re(h''/h - (h'/h)^2),

con h'' dall'ODE.  Nessuna differenza finita su una funzione che oscilla
cinquecento volte nella finestra.

Previsione: Q_M dalla serie WKB sul potenziale complesso, con la frequenza
perturbata di una frazione relativa delta.

Uso: python3.13 calculations/madelung_conditioning_schwarzschild.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from leaver_qnm import leaver_qnm  # noqa: E402
from madelung_wkb_prediction import predict  # noqa: E402
from vaidya_numerator_factored import horizon_values  # noqa: E402
from vaidya_solvability import tortoise  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402


def exact_quantum_potential(ell, spin, omega, width, r_max, points, terms=None):
    """Q_M dalla soluzione esatta, con derivate analitiche e BC di Frobenius.

    Restituisce (r, x_*, Q_M) su griglia UNIFORME IN x_*.
    """
    scale = ell + 0.5
    epsilon = 1.0 / scale
    if terms is None:
        # I coefficienti di Frobenius crescono con ell(ell+1): a ell=70 sessanta
        # termini a x=0.5 lasciano un residuo che si propaga come modo spurio e
        # fa oscillare Q_M di 5.9e-5 contro un valore di 1.9e-6.  Con 120
        # l'ampiezza scende a 1.7e-8.  Regola: due termini per multipolo.
        terms = max(60, 2 * ell)

    def curvature(r, h, dh):
        lapse = 1.0 - 2.0 / r
        potential = ell * (ell + 1) / r**2 + 2.0 * (1 - spin**2) / r**3
        return (potential * h - (2.0 / r**2 - 2j * omega) * dh) / lapse

    # Integrazione in x_*, con r come variabile di stato: la griglia deve essere
    # uniforme in x_* perche' le differenze finite di `predict` siano valide.
    # Con una griglia uniforme in r il passo dx_* = dr/f varia di un fattore 5
    # fra r=2.5 e r=90, e l'errore entra allo stesso ordine di Q_M.
    start = 2.0 + width
    h0, dh0 = horizon_values(ell, spin, omega, width, terms)
    grid = np.linspace(float(tortoise(start)), float(tortoise(r_max)), points)

    def rhs(_t, y):
        h, dh, r = y
        radius = float(r.real)
        lapse = 1.0 - 2.0 / radius
        return np.array([lapse * dh, lapse * curvature(radius, h, dh), lapse],
                        dtype=complex)

    solution = solve_ivp(
        rhs, (grid[0], grid[-1]),
        np.array([h0, dh0, start], dtype=complex),
        t_eval=grid, method="DOP853", rtol=1.0e-13, atol=1.0e-16,
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    h, dh = solution.y[0], solution.y[1]
    radius = solution.y[2].real
    ddh = curvature(radius, h, dh)
    lapse, dlapse = 1.0 - 2.0 / radius, 2.0 / radius**2

    logarithmic = dh / h
    first = omega.imag + lapse * logarithmic.real
    second = (lapse * dlapse * logarithmic.real
              + lapse**2 * (ddh / h - logarithmic**2).real)
    return radius, grid, -(epsilon**2) * (second + first**2)


def conditioning(
    ell: int = 70,
    spin: int = 2,
    window: tuple[float, float] = (20.0, 50.0),
    deltas: tuple[float, ...] = (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2),
    width: float = 0.3,
    r_max: float = 90.0,
    points: int = 200001,
) -> dict:
    omega = leaver_qnm(ell, 0, spin)
    scale = ell + 0.5
    epsilon = 1.0 / scale

    r, xstar, reference = exact_quantum_potential(
        ell, spin, omega, width, r_max, points)
    step = float(xstar[1] - xstar[0])

    lapse = 1.0 - 2.0 / r
    potential = lapse * (ell * (ell + 1) / r**2 + 2.0 * (1 - spin**2) / r**3)
    mask = (r > window[0]) & (r < window[1])

    rows = []
    for delta in deltas:
        q = ((omega * (1.0 + delta)) / scale) ** 2 - potential / scale**2
        predicted = predict(q, epsilon, step, 2)
        rows.append((delta, float(np.median(
            np.abs(predicted - reference)[mask] / np.abs(reference[mask])))))
    return {"omega": omega, "epsilon": epsilon, "rows": rows}


def main() -> None:
    print("Condizionamento di Q_M rispetto alla frequenza — Schwarzschild, s=2")
    print("Riferimento: Q_M esatto, BC di Frobenius e derivate analitiche.")
    print()
    for ell in (20, 40, 70, 100):
        out = conditioning(ell=ell)
        floor = out["rows"][0][1]
        print(f"  ell={ell:3d}   M*omega = {out['omega'].real:.6f}{out['omega'].imag:+.6f}i"
              f"   pavimento = {floor:.3e}")
        print("     delta      errore      amplificazione")
        for delta, error in out["rows"][1:]:
            print(f"    {delta:.0e}   {error:.3e}      {(error - floor) / delta:9.1f}")
        print()


if __name__ == "__main__":
    main()
