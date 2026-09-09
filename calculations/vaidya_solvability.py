#!/usr/bin/env python3
"""Condizione di solvibilita' per la correzione non adiabatica di Vaidya.

DOMANDA.  A ordine Mdot il residuo dell'ansatz adiabatico e' 2 d_r d_M Z
(vedi `research/vaidya_adiabatic_2026-09-09.md`).  La correzione Z_1 risolve

    L_M Z_1 = -2 d_r d_M Z,

con L_M l'operatore congelato.  Ma alla frequenza QNM L_M e' **singolare**: Z
sta nel suo nucleo.  La forzatura e' quindi risonante, e la componente della
sorgente parallela al modo non produce Z_1 ma uno spostamento di frequenza.
Solo la componente ortogonale genera una correzione vera, ed e' li' che una
dipendenza dalla storia potrebbe entrare.

Questo script calcola la proiezione.

IMPOSTAZIONE.  Da 2 psi_vr + d_r(f psi_r) - U psi = 0, con
psi = Z(r;M) exp(-i int omega dv), l'ordine Mdot^0 da'

    L_M Z = f Z'' + (f' - 2 i omega) Z' - U Z = 0,     ' = d/dr,

che con il fattore integrante mu = exp(-2 i omega r_*) e' in forma di
Sturm-Liouville: L_M e' autoaggiunto rispetto al prodotto **bilineare** (senza
coniugazione) <u,v> = int mu u v dr.  E' il prodotto corretto per problemi di
risonanza non autoaggiunti.

SCALA.  Schwarzschild ha un solo parametro, quindi d_M e' una derivata di scala:
con x=r/M, Omega=M omega, Z(r;M)=G(x) e G = exp(i Omega x_*) R,

    d_M Z = -(x/M) G',      d_r d_M Z = -(1/M^2) (x G')'.

La proiezione da valutare e' dunque, a meno di fattori di M,

    I = int mu Z (x G')' dx = int exp(-2 i Omega x_*) G (x G')' dx.

REGOLARIZZAZIONE.  L'integrale QNM non converge: R diverge a entrambi i bordi.
Lo script non finge il contrario: calcola I su finestre crescenti e riporta se
converge, diverge o si riduce a termini di bordo.

Uso: python3.13 vaidya_solvability.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from leaver_qnm import leaver_qnm  # noqa: E402
from schwarzschild_wkb import schwarzschild_potential  # noqa: E402


def tortoise(x: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(x) + 2.0 * np.log(np.asarray(x) / 2.0 - 1.0)


def radial_mode(
    ell: int,
    spin: int,
    omega: complex,
    x_min: float,
    x_max: float,
    points: int,
) -> dict[str, np.ndarray]:
    """R(x) e derivate, da integrazione con condizione entrante all'orizzonte."""
    xstar = np.linspace(float(tortoise(x_min)), float(tortoise(x_max)), points)
    psi0 = np.exp(-1j * omega * xstar[0])

    def rhs(_t: float, state: np.ndarray) -> np.ndarray:
        psi, chi, x = state
        xr = float(x.real)
        f = 1.0 - 2.0 / xr
        return np.array(
            [chi, (schwarzschild_potential(xr, ell, spin) - omega**2) * psi, f],
            dtype=complex,
        )

    solution = solve_ivp(
        rhs,
        (xstar[0], xstar[-1]),
        np.array([psi0, -1j * omega * psi0, x_min], dtype=complex),
        t_eval=xstar,
        method="DOP853",
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return {"x": solution.y[2].real, "xstar": xstar, "R": solution.y[0], "dR": solution.y[1]}


def projection(
    ell: int,
    spin: int = 0,
    window: tuple[float, float] = (4.0, 20.0),
    x_min: float = 2.0001,
    x_max: float = 60.0,
    points: int = 24000,
) -> dict[str, complex | float]:
    """I = int exp(-2 i Omega x_*) G (x G')' dx sulla finestra data."""
    omega = leaver_qnm(ell, 0, spin)
    data = radial_mode(ell, spin, omega, x_min, x_max, points)
    x, xstar, R, dR_star = data["x"], data["xstar"], data["R"], data["dR"]

    f = 1.0 - 2.0 / x
    phase = np.exp(1j * omega * xstar)
    amplitude = phase * R                      # G
    # G' rispetto a x:  dR/dx = (1/f) dR/dx_*
    d_amplitude = phase * (1j * omega * R / f + dR_star / f)

    product = x * d_amplitude
    d_product = np.gradient(product, x)        # (x G')'
    integrand = np.exp(-2j * omega * xstar) * amplitude * d_product

    mask = (x > window[0]) & (x < window[1])
    mask[:4] = False
    mask[-4:] = False
    value = np.trapezoid(integrand[mask], x[mask])
    scale = np.trapezoid(np.abs(integrand[mask]), x[mask])
    return {
        "ell": float(ell),
        "omega": omega,
        "I": complex(value),
        "abs_I": abs(value),
        "scale": float(scale),
        "ratio": abs(value) / scale if scale else float("nan"),
    }


def main() -> None:
    print("Proiezione della sorgente non adiabatica sul modo QNM")
    print("I = int exp(-2 i Omega x_*) G (x G')' dx,  prodotto bilineare con peso mu")
    print()
    print("=== L'integrale converge allargando la finestra? ===")
    print("  ell   finestra        |I|            |I|/int|integrando|")
    for ell in (2, 4):
        for hi in (10.0, 20.0, 30.0, 40.0):
            row = projection(ell, 0, window=(4.0, hi))
            print(
                f"  {ell:3d}   4 < x < {hi:4.0f}   {row['abs_I']:.6e}    {row['ratio']:.4f}"
            )
        print()

    print("=== La cancellazione e' esatta o parziale? ===")
    print("Se I fosse nullo per struttura, |I|/int|integrando| sarebbe ~0.")
    print("Un rapporto O(1) significa che la sorgente NON e' ortogonale al modo.")
    print()
    for ell in (2, 3, 4, 6):
        row = projection(ell, 0, window=(4.0, 25.0))
        print(
            f"  ell={row['ell']:.0f}  M*omega = {row['omega'].real:.6f}"
            f"{row['omega'].imag:+.6f}i   |I| = {row['abs_I']:.4e}   "
            f"rapporto = {row['ratio']:.4f}"
        )


if __name__ == "__main__":
    main()
