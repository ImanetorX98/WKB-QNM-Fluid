#!/usr/bin/env python3
"""Condizionamento di Q_M rispetto alla frequenza, su Schwarzschild.

MOTIVO.  Il §5 del manoscritto afferma che il potenziale quantistico amplifica di
~10^2 l'errore sulla frequenza, sulla base di una misura a ell=70 riportata in
`research/kerr_madelung_sensitivity_2026-09-09.md`.  **Il codice che la produce
non era nel repository**, e la tabella non era quindi riproducibile: questo
modulo colma la lacuna.

L'ha resa evidente il tentativo di disegnarne la figura.  Ripetendo la stessa
misura su Poschl-Teller, dove Q_M e' noto in forma chiusa, l'amplificazione
risulta **1.41 e costante in L** — cioe' assente.  Poschl-Teller e' un caso
speciale: Q_M = -(eps^2/4)(1+sech^2 y) non dipende affatto dalla frequenza, e
percio' non puo' esibire il fenomeno.  La misura va fatta dove il fenomeno vive.

COSTRUZIONE.  Con eps = 1/L, L = ell+1/2, coordinata tortoise x_*:

    q  = (omega/L)^2 - V/L^2,
    u  = sqrt(q) + eps^2 u_2,
    ln A = -1/2 ln|u| - (1/eps) \\int Im(u) dx_*,
    Q_M^pred = -eps^2 [ (ln A)'' + ((ln A)')^2 ].

Il riferimento e' Q_M estratto dalla soluzione numerica dell'equazione master
alla frequenza **esatta** di Leaver: A = |psi|, Q_M = -eps^2 A''/A.  Si perturba
la frequenza nella sola previsione e si misura lo scarto mediano nella finestra
oscillatoria esterna alla barriera.

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
from vaidya_solvability import radial_mode, tortoise  # noqa: E402


def conditioning(
    ell: int = 70,
    spin: int = 2,
    window: tuple[float, float] = (20.0, 50.0),
    deltas: tuple[float, ...] = (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2),
    x_max: float = 90.0,
    points: int = 120001,
) -> dict:
    """Errore su Q_M quando la frequenza usata nella previsione e' sbagliata."""
    omega = leaver_qnm(ell, 0, spin)
    scale = ell + 0.5
    epsilon = 1.0 / scale

    data = radial_mode(ell, spin, omega, 2.0001, x_max, points)
    r, xstar = data["x"], data["xstar"]
    step = float(xstar[1] - xstar[0])

    amplitude = np.abs(data["R"])
    reference = -(epsilon**2) * np.gradient(np.gradient(amplitude, step), step) / amplitude

    # Langer: ell(ell+1) -> L^2.  Non e' cosmetico — la differenza vale
    # -1/(4L^2) volte f/r^2, cioe' **lo stesso ordine di Q_M**, e usando
    # ell(ell+1) lo scarto relativo cresce come L^2 invece di calare.
    lapse = 1.0 - 2.0 / r
    potential = lapse * (scale**2 / r**2 + 2.0 * (1 - spin**2) / r**3)
    mask = (r > window[0]) & (r < window[1])

    rows = []
    for delta in deltas:
        q = ((omega * (1.0 + delta)) / scale) ** 2 - potential / scale**2
        predicted = predict(q, epsilon, step, 2)
        rows.append((delta, float(np.median(
            np.abs(predicted - reference)[mask] / np.abs(reference[mask])))))
    return {"omega": omega, "epsilon": epsilon, "rows": rows}


def main() -> None:
    print("Condizionamento di Q_M rispetto alla frequenza — Schwarzschild")
    print("Riferimento: Q_M dalla soluzione numerica alla frequenza esatta di Leaver.")
    print()
    for ell in (20, 40, 70, 100):
        out = conditioning(ell=ell)
        floor = out["rows"][0][1]
        print(f"  ell={ell:3d}   M*omega = {out['omega'].real:.6f}{out['omega'].imag:+.6f}i"
              f"   pavimento = {floor:.2e}")
        print("     delta      errore     amplificazione")
        for delta, error in out["rows"][1:]:
            amp = (error - floor) / delta
            print(f"    {delta:.0e}   {error:.3e}     {amp:8.1f}")
        print()


if __name__ == "__main__":
    main()
