#!/usr/bin/env python3
"""Integra il modo radiale e rende esplicito il potenziale di Madelung.

Normalizzazione Schrödinger-like:
    -1/2 psi'' + (V/2) psi = (Omega^2/2) psi,
dove i primi sono rispetto a x_*=r_*/M.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from schwarzschild_wkb import qnm_wkb, schwarzschild_potential


def tortoise(x: np.ndarray | float) -> np.ndarray | float:
    """x_*=r_*/M per x=r/M>2; la costante additiva e' convenzionale."""
    return np.asarray(x) + 2.0 * np.log(np.asarray(x) / 2.0 - 1.0)


def integrate_profile(
    ell: int,
    overtone: int,
    spin: int,
    order: int,
    x_min: float,
    x_max: float,
    points: int,
) -> dict[str, np.ndarray | complex]:
    if not 2.0 < x_min < x_max:
        raise ValueError("servono 2 < x_min < x_max")
    if points < 50:
        raise ValueError("points deve essere almeno 50")

    result = qnm_wkb(ell, overtone, spin, order)
    omega = result.omega_M
    omega_sq = omega**2
    x_eval = np.linspace(x_min, x_max, points)
    xstar_min = float(tortoise(x_min))

    # Condizione locale ingoing all'orizzonte per exp(-i omega t):
    # psi ~ exp(-i Omega x_*), chi=dpsi/dx_*=-i Omega psi.
    psi0 = np.exp(-1j * omega * xstar_min)
    chi0 = -1j * omega * psi0

    def rhs(x: float, y: np.ndarray) -> np.ndarray:
        psi, chi = y
        f = 1.0 - 2.0 / x
        potential = schwarzschild_potential(x, ell, spin)
        return np.array([chi / f, (potential - omega_sq) * psi / f], dtype=complex)

    solution = solve_ivp(
        rhs,
        (x_min, x_max),
        np.array([psi0, chi0], dtype=complex),
        t_eval=x_eval,
        rtol=2.0e-10,
        atol=2.0e-12,
        method="DOP853",
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    psi = solution.y[0]
    chi = solution.y[1]
    xstar = tortoise(x_eval)
    potential = schwarzschild_potential(x_eval, ell, spin)
    amplitude = np.abs(psi)

    # y=psi'/psi=A'/A+i S'. Dall'ODE:
    # A''/A = Re(V-Omega^2) + (S')^2.
    logarithmic_derivative = chi / psi
    phase_gradient = np.imag(logarithmic_derivative)
    q_madelung = -0.5 * (np.real(potential - omega_sq) + phase_gradient**2)
    external = 0.5 * potential
    kinetic = 0.5 * phase_gradient**2
    hj_sum = kinetic + external + q_madelung
    target = np.full_like(x_eval, 0.5 * np.real(omega_sq))

    return {
        "x": x_eval,
        "xstar": xstar,
        "psi": psi,
        "amplitude": amplitude,
        "potential": potential,
        "phase_gradient": phase_gradient,
        "q_madelung": q_madelung,
        "external": external,
        "kinetic": kinetic,
        "hj_sum": hj_sum,
        "target": target,
        "omega": omega,
    }


def save_csv(data: dict[str, np.ndarray | complex], path: Path) -> None:
    psi = np.asarray(data["psi"])
    table = np.column_stack(
        [
            data["x"],
            data["xstar"],
            psi.real,
            psi.imag,
            data["amplitude"],
            data["potential"],
            data["phase_gradient"],
            data["external"],
            data["kinetic"],
            data["q_madelung"],
            data["hj_sum"],
            data["target"],
        ]
    )
    header = "x,xstar,Re_psi,Im_psi,A,V,Sprime,V_over_2,K_phase,Q_M,HJ_sum,Re_Omega2_over_2"
    np.savetxt(path, table, delimiter=",", header=header, comments="")


def save_plot(data: dict[str, np.ndarray | complex], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(path.parent / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xstar = np.asarray(data["xstar"])
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    axes[0].plot(xstar, data["potential"], label=r"$V_s$")
    axes[0].set_ylabel(r"$M^2 V_s$")
    axes[0].legend()

    amplitude = np.asarray(data["amplitude"])
    axes[1].semilogy(xstar, amplitude / amplitude.max(), label=r"$A/\max A$")
    axes[1].set_ylabel("ampiezza normalizzata")
    axes[1].legend()

    axes[2].plot(xstar, data["external"], label=r"$V_s/2$")
    axes[2].plot(xstar, data["kinetic"], label=r"$(S')^2/2$")
    axes[2].plot(xstar, data["q_madelung"], label=r"$Q_M=-A''/(2A)$")
    axes[2].plot(xstar, data["hj_sum"], "k--", linewidth=1.2, label="somma HJ")
    axes[2].set_xlabel(r"$x_*=r_*/M$")
    axes[2].set_ylabel("energia adimensionale")
    axes[2].legend(ncol=2)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spin", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--ell", type=int, default=2)
    parser.add_argument("--n", type=int, default=0, dest="overtone")
    parser.add_argument("--order", type=int, choices=(1, 3), default=3)
    parser.add_argument("--x-min", type=float, default=2.02)
    parser.add_argument("--x-max", type=float, default=35.0)
    parser.add_argument("--points", type=int, default=2500)
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = integrate_profile(
        args.ell,
        args.overtone,
        args.spin,
        args.order,
        args.x_min,
        args.x_max,
        args.points,
    )
    stem = f"madelung_s{args.spin}_l{args.ell}_n{args.overtone}_wkb{args.order}"
    csv_path = args.output_dir / f"{stem}.csv"
    save_csv(data, csv_path)
    if not args.no_plot:
        save_plot(data, args.output_dir / f"{stem}.png")

    omega = complex(data["omega"])
    residual = np.max(np.abs(np.asarray(data["hj_sum"]) - np.asarray(data["target"])))
    print(f"M omega = {omega.real:.12f} {omega.imag:+.12f} i")
    print(f"CSV: {csv_path}")
    if not args.no_plot:
        print(f"PNG: {args.output_dir / f'{stem}.png'}")
    print(f"max |K + V/2 + Q_M - Re(Omega^2)/2| = {residual:.3e}")


if __name__ == "__main__":
    main()
