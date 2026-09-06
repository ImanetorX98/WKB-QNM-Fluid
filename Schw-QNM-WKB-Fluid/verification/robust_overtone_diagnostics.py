#!/usr/bin/env python3
"""Diagnostica ben condizionata del degrado di regolarita' con l'overtone.

Il test di sensibilita' (vedi manuscript, sez. 7bis) mostra che la profondita'
dei quasi-zeri di psi cambia di un ordine di grandezza per un errore dell'1% su
omega: e' il residuo di una quasi-cancellazione fra rami, mal condizionato per
costruzione.  Numero e posizione degli zeri sono invece stabili.

Si costruiscono quindi tre osservabili che non guardano il massimo puntuale:

  zero_count  numero di minimi locali di |psi| (misura di conteggio)
  fd_l1       ||Q_M - Q_M_fd||_1 / ||Q_M||_1 sul profilo (misura integrale)
  bad_frac    frazione di punti in cui |Q_M - Q_M_fd| > 0.1 |Q_M|
              (misura dell'insieme guasto, non della sua profondita')

Tutte e tre integrano o contano invece di campionare il punto peggiore.
Lo script verifica prima che siano stabili sotto perturbazione di omega, poi
che discriminino ancora l'overtone.

Uso: python3.13 robust_overtone_diagnostics.py
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


def integrate_at_omega(
    ell: int,
    omega: complex,
    spin: int = 2,
    x_min: float = 2.0001,
    x_max: float = 60.0,
    points: int = 12000,
) -> dict[str, np.ndarray | float]:
    """Come scalar_eikonal_scaling.integrate_profile ma con omega imposta.

    Serve a poter perturbare omega senza passare dalla formula WKB.
    """
    xstar_min, xstar_max = float(tortoise(x_min)), float(tortoise(x_max))
    xstar = np.linspace(xstar_min, xstar_max, points)
    omega_sq = omega**2
    psi0 = np.exp(-1j * omega * xstar_min)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        psi, chi, x = y
        f = 1.0 - 2.0 / x.real
        return np.array(
            [chi, (schwarzschild_potential(x.real, ell, spin) - omega_sq) * psi, f],
            dtype=complex,
        )

    solution = solve_ivp(
        rhs,
        (xstar_min, xstar_max),
        np.array([psi0, -1j * omega * psi0, x_min], dtype=complex),
        t_eval=xstar,
        rtol=2.0e-10,
        atol=2.0e-12,
        method="DOP853",
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    psi, chi = solution.y[0], solution.y[1]
    x_eval = solution.y[2].real
    epsilon = 1.0 / (ell + 0.5)
    amplitude = np.abs(psi)

    f = 1.0 - 2.0 / x_eval
    h_squared = f / x_eval**2
    v2 = f * (2.0 * (1.0 - spin**2) / x_eval**3 - 0.25 / x_eval**2)
    momentum = epsilon * np.imag(chi / psi)
    target = np.full_like(x_eval, float(np.real(omega_sq)) * epsilon**2)
    q_madelung = target - h_squared - epsilon**2 * v2 - momentum**2

    step = xstar[1] - xstar[0]
    q_fd = -(epsilon**2) * np.gradient(np.gradient(amplitude, step), step) / amplitude
    return {
        "x": x_eval,
        "xstar": xstar,
        "amplitude": amplitude,
        "q_madelung": q_madelung,
        "q_fd": q_fd,
        "step": step,
    }


def robust_observables(data: dict[str, np.ndarray | float]) -> dict[str, float]:
    amplitude = np.asarray(data["amplitude"])
    q = np.asarray(data["q_madelung"])
    q_fd = np.asarray(data["q_fd"])
    step = float(data["step"])
    interior = slice(5, -5)

    zero_count = int(
        np.sum((amplitude[1:-1] < amplitude[:-2]) & (amplitude[1:-1] < amplitude[2:]))
    )
    discrepancy = np.abs(q - q_fd)[interior]
    magnitude = np.abs(q)[interior]
    fd_l1 = float(np.trapezoid(discrepancy, dx=step) / np.trapezoid(magnitude, dx=step))
    bad_frac = float(np.mean(discrepancy > 0.1 * magnitude))
    return {"zero_count": float(zero_count), "fd_l1": fd_l1, "bad_frac": bad_frac}


def stability_check(ell: int, overtone: int, spin: int = 2) -> None:
    omega = qnm_wkb(ell, overtone, spin, 3).omega_M
    print(f"  ell={ell}, n={overtone}")
    print("    perturbazione   zeri   fd_l1        bad_frac")
    reference: dict[str, float] | None = None
    for label, factor in (("nessuna", 1.0), ("+0.1%", 1.001), ("+1%", 1.01), ("+5%", 1.05)):
        values = robust_observables(integrate_at_omega(ell, omega * factor, spin))
        if reference is None:
            reference = values
        spread = max(
            abs(values[k] - reference[k]) / max(reference[k], 1e-12)
            for k in ("fd_l1", "bad_frac")
        )
        print(
            f"    {label:13} {values['zero_count']:4.0f}   {values['fd_l1']:.4e}   "
            f"{values['bad_frac']:.4f}   (scarto max {spread:.1%})"
        )


def overtone_trend(ell: int, spin: int = 2, n_max: int = 4) -> None:
    print(f"  ell={ell}")
    print("    n    zeri   fd_l1        bad_frac")
    for overtone in range(n_max + 1):
        omega = qnm_wkb(ell, overtone, spin, 3).omega_M
        values = robust_observables(integrate_at_omega(ell, omega, spin))
        print(
            f"    {overtone:2d}   {values['zero_count']:4.0f}   "
            f"{values['fd_l1']:.4e}   {values['bad_frac']:.4f}"
        )


def main() -> None:
    print("=== 1. Le osservabili sono stabili sotto perturbazione di omega? ===")
    for ell, overtone in ((16, 1), (8, 2)):
        stability_check(ell, overtone)
        print()

    print("=== 2. Discriminano ancora l'overtone? ===")
    for ell in (8, 16):
        overtone_trend(ell)
        print()


if __name__ == "__main__":
    main()
