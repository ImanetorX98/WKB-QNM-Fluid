#!/usr/bin/env python3
"""Scaling eikonale del termine di Madelung per Regge-Wheeler su Schwarzschild.

Con L=ell+1/2, eps=1/L, omega=L*Omega, il potenziale RW si riscrive

    V_s/L^2 = h^2 + eps^2 v2,
    h^2 = f/x^2,   v2 = f[2(1-s^2)/x^3 - 1/(4x^2)],   f = 1-2/x,

quindi la correzione al termine eikonale entra a ordine eps^2: nel settore
bosonico NON esiste un termine di ordine eps^1.  La decomposizione di Madelung
psi = A exp(iS/eps), P = S', chiude come

    Re(Omega^2) = P^2 + h^2 + eps^2 v2 + Q_M,     Q_M = -eps^2 A''/A,

e Q_M deve mostrare pendenza 2 in eps, la stessa di eps^2 v2.  E' il contrasto
col caso di Dirac (`core/dirac_madelung_profile.py`), dove un termine di ordine
eps^1 esiste ed e' la connessione di spin.

Uso: python3.13 scalar_eikonal_scaling.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "core"))

from schwarzschild_wkb import qnm_wkb, schwarzschild_potential  # noqa: E402


def tortoise(x: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(x) + 2.0 * np.log(np.asarray(x) / 2.0 - 1.0)


def eikonal_terms(x: np.ndarray, spin: int) -> tuple[np.ndarray, np.ndarray]:
    """(h^2, v2): termine eikonale e coefficiente della correzione eps^2."""
    f = 1.0 - 2.0 / x
    h_squared = f / x**2
    v2 = f * (2.0 * (1.0 - spin**2) / x**3 - 0.25 / x**2)
    return h_squared, v2


def integrate_profile(
    ell: int,
    overtone: int = 0,
    spin: int = 2,
    order: int = 3,
    x_min: float = 2.0001,
    x_max: float = 60.0,
    points: int = 8000,
) -> dict[str, np.ndarray | complex | float]:
    """Integra psi in x_* uniforme e decompone nei termini della chiusura."""
    result = qnm_wkb(ell, overtone, spin, order)
    omega = result.omega_M
    omega_sq = omega**2
    xstar_min = float(tortoise(x_min))
    xstar_max = float(tortoise(x_max))
    xstar_eval = np.linspace(xstar_min, xstar_max, points)

    psi0 = np.exp(-1j * omega * xstar_min)
    chi0 = -1j * omega * psi0

    def rhs(_xstar: float, y: np.ndarray) -> np.ndarray:
        psi, chi, x = y
        f = 1.0 - 2.0 / x.real
        potential = schwarzschild_potential(x.real, ell, spin)
        return np.array([chi, (potential - omega_sq) * psi, f], dtype=complex)

    solution = solve_ivp(
        rhs,
        (xstar_min, xstar_max),
        np.array([psi0, chi0, x_min], dtype=complex),
        t_eval=xstar_eval,
        rtol=2.0e-10,
        atol=2.0e-12,
        method="DOP853",
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    psi = solution.y[0]
    chi = solution.y[1]
    x_eval = solution.y[2].real

    large_l = ell + 0.5
    epsilon = 1.0 / large_l
    amplitude = np.abs(psi)
    h_squared, v2 = eikonal_terms(x_eval, spin)
    momentum = epsilon * np.imag(chi / psi)
    subleading = epsilon**2 * v2
    target = np.full_like(x_eval, float(np.real(omega_sq)) * epsilon**2)
    q_madelung = target - h_squared - subleading - momentum**2

    step = xstar_eval[1] - xstar_eval[0]
    amplitude_second = np.gradient(np.gradient(amplitude, step), step)
    q_madelung_fd = -(epsilon**2) * amplitude_second / amplitude
    interior = slice(2, -2)

    return {
        "x": x_eval,
        "xstar": xstar_eval,
        "amplitude": amplitude,
        "momentum": momentum,
        "h_squared": h_squared,
        "subleading": subleading,
        "q_madelung": q_madelung,
        "q_madelung_fd": q_madelung_fd,
        "target": target,
        "omega": omega,
        "epsilon": float(epsilon),
        "ell": int(ell),
        "fd_residual": float(
            np.max(np.abs(q_madelung[interior] - q_madelung_fd[interior]))
        ),
    }


def scaling_table(
    ell_values: tuple[int, ...] = (2, 4, 8, 16, 32),
    spin: int = 2,
    overtone: int = 0,
    far_window: tuple[float, float] = (20.0, 50.0),
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for ell in ell_values:
        data = integrate_profile(ell, overtone, spin)
        x = np.asarray(data["x"])
        mask = (x > far_window[0]) & (x < far_window[1])
        rows.append(
            {
                "ell": float(ell),
                "epsilon": float(data["epsilon"]),
                "momentum_min": float(np.min(np.abs(np.asarray(data["momentum"])[mask]))),
                "subleading_max": float(np.max(np.abs(np.asarray(data["subleading"])[mask]))),
                "q_madelung_max": float(np.max(np.abs(np.asarray(data["q_madelung"])[mask]))),
                "fd_residual": float(data["fd_residual"]),
            }
        )
        rows[-1]["spin_over_madelung"] = (
            rows[-1]["subleading_max"] / rows[-1]["q_madelung_max"]
        )
    return rows


def slopes(rows: list[dict[str, float]]) -> dict[str, float]:
    log_eps = np.log([row["epsilon"] for row in rows])
    return {
        "subleading_exponent": float(
            np.polyfit(log_eps, np.log([r["subleading_max"] for r in rows]), 1)[0]
        ),
        "madelung_exponent": float(
            np.polyfit(log_eps, np.log([r["q_madelung_max"] for r in rows]), 1)[0]
        ),
    }


def main() -> None:
    label = {0: "scalare", 1: "elettromagnetico", 2: "gravitazionale assiale"}
    for spin in (0, 1, 2):
        rows = scaling_table(spin=spin)
        exponents = slopes(rows)
        coefficient = 2 * (1 - spin**2)
        print(f"--- Regge-Wheeler, s={spin} ({label[spin]}), n=0, finestra 20<x<50 ---")
        print(f"    v2 = f[{coefficient:+d}/x^3 - 1/(4x^2)]")
        print(" ell     eps      min|P|     max|eps^2 v2|     max|Q_M|    v2/Q_M    res.FD")
        for row in rows:
            print(
                f"{row['ell']:4.0f}  {row['epsilon']:7.4f}  {row['momentum_min']:9.6f}  "
                f"{row['subleading_max']:14.6e}  {row['q_madelung_max']:12.6e}  "
                f"{row['spin_over_madelung']:7.4f}  {row['fd_residual']:.1e}"
            )
        print(
            f"pendenza in eps:  eps^2 v2 = {exponents['subleading_exponent']:.4f}, "
            f"Q_M = {exponents['madelung_exponent']:.4f}   (entrambe attese 2)"
        )
        print()
    print("Lo spin non entra in h^2: tutta la sua traccia vive a ordine eps^2,")
    print("cioe' esattamente dove compare il potenziale di Madelung.")
    print("Nessun termine di ordine eps^1 esiste nel settore bosonico.")


if __name__ == "__main__":
    main()
