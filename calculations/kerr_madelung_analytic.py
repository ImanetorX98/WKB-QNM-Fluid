#!/usr/bin/env python3
"""Q_M su Kerr con derivate analitiche: audit del §9.6 del manoscritto.

MOTIVO.  Il §9.6 attribuisce il miglioramento della previsione di Q_M, passando
dalla frequenza eikonale a quella autoconsistente, all'amplificazione ~10^2 del
§5.  **Quell'amplificazione e' falsa: vale 2**
(`madelung_conditioning_schwarzschild.py`).  Inoltre `kerr_madelung_profile.py`
calcola A''/A con `np.gradient` applicato due volte a |psi|, che e' il difetto
numerico piu' grave trovato nell'audit di Schwarzschild.

Questo modulo rifa' la stessa misura con derivate **analitiche** e verifica
separatamente la sensibilita' alla condizione al bordo.

--------------------------------------------------------------------------
COSTRUZIONE
--------------------------------------------------------------------------
La ODE integrata da' psi e chi = dpsi/dr_*, e psi'' = -q psi.  Quindi, con
L = log|psi| e ' = d/dr_*,

    L'  = Re(chi/psi),
    L'' = Re(-q - (chi/psi)^2),
    Q_M = -eps^2 [ L'' + (L')^2 ] .

Nessuna differenza finita.  Il confronto fra questa e la versione a
`np.gradient` isola quanto del risultato del §9.6 fosse numerico.

Uso: python3.13 calculations/kerr_madelung_analytic.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from kerr_madelung_profile import eikonal_frequency, kerr_profile  # noqa: E402
from kerr_eikonal_order_test import spheroidal_eigenvalue  # noqa: E402


def analytic_profile(ell, spin, mu, omega=None, horizon_offset=1.0e-6,
                     points=60001, r_max_factor=18.0):
    """Q_M da derivate analitiche, stessa geometria di `kerr_profile`."""
    large_l = ell + 0.5
    epsilon = 1.0 / large_l
    m = int(round(mu * large_l))
    if omega is None:
        omega = eikonal_frequency(spin, mu, 0, large_l)[0]
    eigen = complex(spheroidal_eigenvalue(ell, m, spin * omega))
    separation = eigen + spin**2 * omega**2 - 2.0 * spin * m * omega

    r_plus = 1.0 + np.sqrt(max(1.0 - spin**2, 0.0))
    reference = kerr_profile(ell, spin, mu)
    r0, r_max = float(reference["r_peak"]), float(reference["r_max"])

    def potential(r):
        h_squared = r**2 + spin**2
        delta = r**2 - 2.0 * r + spin**2
        h = np.sqrt(h_squared)
        d2h = (delta / h_squared * (1.0 / h - r**2 / h**3)
               - 2.0 * r * delta / h_squared**2 * (r / h)
               + (2.0 * r - 2.0) / h_squared * (r / h))
        return ((omega - m * spin / h_squared) ** 2
                - delta * separation / h_squared**2 - d2h / h)

    def rhs(_t, y):
        psi, chi, r = y
        radius = float(r.real)
        h_squared = radius**2 + spin**2
        delta = radius**2 - 2.0 * radius + spin**2
        return np.array([chi, -potential(radius) * psi, delta / h_squared],
                        dtype=complex)

    r_start = r_plus + horizon_offset * (r0 - r_plus)
    sample = np.linspace(r_start, r_max, 40000)
    span = float(np.trapezoid(
        (sample**2 + spin**2) / (sample**2 - 2.0 * sample + spin**2), sample))
    horizon_frequency = spin / (r_plus**2 + spin**2)
    wavenumber = omega - m * horizon_frequency

    grid = np.linspace(0.0, span, points)
    solution = solve_ivp(rhs, (grid[0], grid[-1]),
                         np.array([1.0 + 0j, -1j * wavenumber, r_start], dtype=complex),
                         t_eval=grid, method="DOP853", rtol=2.0e-12, atol=2.0e-14)
    if not solution.success:
        raise RuntimeError(solution.message)

    psi, chi = solution.y[0], solution.y[1]
    radius = solution.y[2].real
    logarithmic = chi / psi
    first = logarithmic.real
    second = (-potential(radius) - logarithmic**2).real
    return {"r": radius, "rstar": grid, "omega": omega, "epsilon": epsilon,
            "q_madelung": -(epsilon**2) * (second + first**2)}


def compare(ell, spin, mu, window=(20.0, 50.0)):
    """Derivate analitiche contro `np.gradient`, sullo stesso modo."""
    analytic = analytic_profile(ell, spin, mu)
    finite = kerr_profile(ell, spin, mu)
    out = {}
    for name, data in (("analitico", analytic), ("np.gradient", finite)):
        r = np.asarray(data["r"])
        mask = (r > window[0]) & (r < window[1])
        out[name] = float(np.median(np.asarray(data["q_madelung"])[mask].real))
    out["scarto relativo"] = abs(out["analitico"] / out["np.gradient"] - 1.0)
    return out


def boundary_sensitivity(ell, spin, mu, offsets=(1e-4, 1e-6, 1e-8),
                         window=(20.0, 50.0)):
    """La BC troncata all'orizzonte lascia traccia?  Difetto 1 dell'audit."""
    values = []
    for offset in offsets:
        data = analytic_profile(ell, spin, mu, horizon_offset=offset)
        r = np.asarray(data["r"])
        mask = (r > window[0]) & (r < window[1])
        values.append(float(np.median(np.asarray(data["q_madelung"])[mask].real)))
    centre = float(np.mean(values))
    return {"valori": values, "dispersione": max(abs(v - centre) for v in values) / abs(centre)}


def main() -> None:
    print("Audit del §9.6: Q_M su Kerr con derivate analitiche")
    print()
    print("1. Derivate analitiche contro np.gradient (finestra 20<r<50)")
    print("     a   mu  ell     analitico      np.gradient     scarto")
    for spin, mu, ell in ((0.0, 0.0, 30), (0.6, 0.5, 30), (0.6, 0.5, 70), (0.9, 0.5, 70)):
        row = compare(ell, spin, mu)
        print(f"   {spin:4.1f} {mu:4.1f} {ell:4d}   {row['analitico']:+.6e}  "
              f"{row['np.gradient']:+.6e}   {row['scarto relativo']:.2e}")
    print()
    print("2. Sensibilita' alla condizione al bordo dell'orizzonte")
    print("     a   mu  ell    dispersione su offset 1e-4..1e-8")
    for spin, mu, ell in ((0.0, 0.0, 30), (0.6, 0.5, 30), (0.9, 0.5, 70)):
        row = boundary_sensitivity(ell, spin, mu)
        print(f"   {spin:4.1f} {mu:4.1f} {ell:4d}    {row['dispersione']:.2e}")


if __name__ == "__main__":
    main()
