#!/usr/bin/env python3
"""Degrado della chiusura di Madelung al crescere dell'overtone.

Il Teorema 1 richiede psi != 0 lungo il contorno: il momento u = -i eps psi'/psi
ha un polo a ogni zero di psi, e Q_M = -eps^2 A''/A vi diverge.  Per omega reale
gli zeri sarebbero sull'asse; per omega complessa sono spostati nel piano, ma
piu' alto e' l'overtone piu' si avvicinano al contorno di integrazione.

L'ipotesi da testare e' quindi che la gerarchia non degradi in modo liscio con n
ma per avvicinamento di zeri, con una firma precisa: |psi| sviluppa minimi
profondi e |Q_M| vi produce picchi, mentre lontano da quelli la chiusura resta
buona.  Se e' cosi', il limite non e' l'accuratezza WKB ma la geometria degli
zeri, cioe' la stessa struttura che l'exact WKB descrive con le linee di Stokes.

Uso: python3.13 overtone_degradation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from scalar_eikonal_scaling import integrate_profile, scaling_table, slopes  # noqa: E402


def envelope_dips(amplitude: np.ndarray) -> tuple[float, int]:
    """Profondita' massima dei minimi di |psi| rispetto all'inviluppo liscio.

    L'ampiezza di un QNM cresce esponenzialmente; per vedere gli zeri serve
    togliere quella tendenza.  Si sottrae una retta ai minimi quadrati da
    log|psi| e si guardano i residui negativi.
    """
    log_amplitude = np.log(amplitude)
    index = np.arange(log_amplitude.size, dtype=float)
    trend = np.polyval(np.polyfit(index, log_amplitude, 1), index)
    residual = log_amplitude - trend
    interior = residual[5:-5]
    local_minima = int(
        np.sum((interior[1:-1] < interior[:-2]) & (interior[1:-1] < interior[2:]))
    )
    return float(-np.min(interior)), local_minima


def overtone_row(ell: int, overtone: int, spin: int = 2) -> dict[str, float]:
    data = integrate_profile(ell, overtone, spin)
    x = np.asarray(data["x"])
    amplitude = np.asarray(data["amplitude"])
    q_madelung = np.abs(np.asarray(data["q_madelung"]))
    omega = complex(data["omega"])
    mask = (x > 20.0) & (x < 50.0)

    depth, minima = envelope_dips(amplitude)
    window = q_madelung[mask]
    return {
        "n": float(overtone),
        "gamma": float(np.imag(omega**2)),
        "ratio_im_re": float(abs(omega.imag / omega.real)),
        "dip_depth": depth,
        "local_minima": float(minima),
        "q_median": float(np.median(window)),
        "q_max": float(np.max(window)),
        "spikiness": float(np.max(window) / np.median(window)),
        "fd_residual": float(data["fd_residual"]),
    }


def main() -> None:
    for ell in (8, 16):
        print(f"=== Regge-Wheeler s=2, ell={ell}, finestra 20<x<50 ===")
        print("  n     Gamma    |Im/Re|   dip max   minimi   med|Q_M|    max/med   res.FD")
        for overtone in range(5):
            row = overtone_row(ell, overtone)
            print(
                f"{row['n']:3.0f}  {row['gamma']:+9.4f}  {row['ratio_im_re']:7.4f}  "
                f"{row['dip_depth']:8.4f}  {row['local_minima']:6.0f}  "
                f"{row['q_median']:10.3e}  {row['spikiness']:8.2f}  {row['fd_residual']:.1e}"
            )
        print()

    print("=== pendenza in eps di Q_M al variare dell'overtone ===")
    print("(ell = 4,8,16,32; se la gerarchia regge, la pendenza resta ~2)")
    print("  n    pendenza eps^2 v2    pendenza Q_M")
    for overtone in range(4):
        rows = scaling_table(ell_values=(4, 8, 16, 32), overtone=overtone)
        exponents = slopes(rows)
        print(
            f"{overtone:3d}  {exponents['subleading_exponent']:17.4f}  "
            f"{exponents['madelung_exponent']:14.4f}"
        )
    print()
    print("Nota: la WKB di barriera perde accuratezza per n >~ ell, quindi le")
    print("righe con n alto vanno lette come struttura, non come QNM esatti.")


if __name__ == "__main__":
    main()
