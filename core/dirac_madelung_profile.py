#!/usr/bin/env python3
"""Profilo di Madelung numerico per Dirac massless su Schwarzschild.

Estende `madelung_profile.py` ai potenziali partner di Darboux

    V_tau = W^2 + sigma D_x W = K^2 h^2 + tau K h',
    W = kappa h,  h = sqrt(f)/x,  f = 1-2/x,  x = r/M,
    D_x = f d/dx  (derivata tortoise),  tau = sigma*sign(kappa).

Nella scalatura eikonale della nota (K=|kappa|, eps=1/K, omega=K*Omega) la
componente scalare obbedisce a

    eps^2 Z'' + [Omega^2 - h^2 - tau eps h'] Z = 0,

e la decomposizione Z = A exp(i S/eps), P = S', da' la chiusura esatta

    Re(Omega^2) = P^2 + h^2 + tau eps h' + Q_M,   Q_M = -eps^2 A''/A.

Il termine di spin connection e' O(eps), il funzionale di Madelung O(eps^2):
il loro rapporto deve quindi crescere linearmente in K.  E' questa la firma
che distingue Dirac dal caso scalare, dove il termine O(eps) e' assente.
"""

from __future__ import annotations

import argparse
import os
from functools import lru_cache
from pathlib import Path
from typing import Callable

import cmath
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar


def tortoise(x: np.ndarray | float) -> np.ndarray | float:
    """x_*=r_*/M per x=r/M>2; la costante additiva e' convenzionale."""
    return np.asarray(x) + 2.0 * np.log(np.asarray(x) / 2.0 - 1.0)


def eikonal_h(x: np.ndarray | float) -> np.ndarray | float:
    """h = sqrt(f)/x, la massa spaziale adimensionale del sistema di Dirac."""
    x_arr = np.asarray(x, dtype=float)
    return np.sqrt(1.0 - 2.0 / x_arr) / x_arr


def eikonal_h_prime(x: np.ndarray | float) -> np.ndarray | float:
    """h' = D_x h = sqrt(f)(3-x)/x^3, con D_x = f d/dx."""
    x_arr = np.asarray(x, dtype=float)
    return np.sqrt(1.0 - 2.0 / x_arr) * (3.0 - x_arr) / x_arr**3


def dirac_potential(x: np.ndarray | float, kappa_abs: float, tau: int) -> np.ndarray | float:
    """M^2 V_tau(x) per Dirac massless; tau = sigma*sign(kappa) = +/-1."""
    if tau not in (1, -1):
        raise ValueError("tau deve essere +1 oppure -1")
    if kappa_abs <= 0:
        raise ValueError("K=|kappa| deve essere positivo")
    h = eikonal_h(x)
    value = kappa_abs**2 * h**2 + tau * kappa_abs * eikonal_h_prime(x)
    return float(value) if np.ndim(value) == 0 else value


@lru_cache(maxsize=None)
def _derivative_functions(tau: int) -> tuple[Callable[[float, float], float], ...]:
    """V_tau e derivate d^k V/dx_*^k, k=1,...,6, come funzioni di (x, K)."""
    x, k = sp.symbols("x K", positive=True)
    f = 1 - 2 / x
    root_f = sp.sqrt(f)
    expr = k**2 * f / x**2 + tau * k * root_f * (3 - x) / x**3
    expressions = [sp.simplify(expr)]
    for _ in range(6):
        expressions.append(sp.simplify(f * sp.diff(expressions[-1], x)))
    return tuple(sp.lambdify((x, k), item, modules="numpy") for item in expressions)


def dirac_potential_peak(kappa_abs: float, tau: int) -> float:
    """Posizione x=r/M del massimo della barriera fuori dall'orizzonte."""
    result = minimize_scalar(
        lambda x: -dirac_potential(x, kappa_abs, tau),
        bounds=(2.0 + 1.0e-8, 50.0),
        method="bounded",
        options={"xatol": 1.0e-14},
    )
    if not result.success:
        raise RuntimeError(f"ricerca del massimo fallita: {result.message}")
    return float(result.x)


def dirac_qnm_wkb(
    kappa_abs: float, overtone: int = 0, tau: int = 1, order: int = 3
) -> tuple[float, complex]:
    """Omega=M*omega al primo o terzo ordine WKB per il partner V_tau.

    Restituisce (x_peak, M*omega).  I due partner sono isospettrali, quindi il
    valore per tau=-1 deve coincidere con quello per tau=+1 entro l'errore WKB.
    """
    if order not in (1, 3):
        raise ValueError("questo modulo implementa order=1 oppure order=3")
    x0 = dirac_potential_peak(kappa_abs, tau)
    values = np.array([float(fn(x0, kappa_abs)) for fn in _derivative_functions(tau)])
    v0, _, v2, v3, v4, v5, v6 = values
    if v2 >= 0:
        raise RuntimeError(f"il punto trovato non e' un massimo: V''={v2}")

    alpha = overtone + 0.5
    root_curvature = np.sqrt(-2.0 * v2)
    lambda2 = 0.0
    lambda3 = 0.0
    if order == 3:
        lambda2 = (
            0.125 * (v4 / v2) * (0.25 + alpha**2)
            - (v3 / v2) ** 2 * (7.0 + 60.0 * alpha**2) / 288.0
        ) / root_curvature
        lambda3 = 1.0 / (-2.0 * v2) * (
            5.0 * (v3 / v2) ** 4 * (77.0 + 188.0 * alpha**2) / 6912.0
            - (v3**2 * v4 / v2**3) * (51.0 + 100.0 * alpha**2) / 384.0
            + (v4 / v2) ** 2 * (67.0 + 68.0 * alpha**2) / 2304.0
            + (v3 * v5 / v2**2) * (19.0 + 28.0 * alpha**2) / 288.0
            - (v6 / v2) * (5.0 + 4.0 * alpha**2) / 288.0
        )

    omega_sq = v0 + root_curvature * lambda2 - 1j * alpha * root_curvature * (1.0 + lambda3)
    omega = cmath.sqrt(omega_sq)
    if omega.real < 0:
        omega = -omega
    if omega.imag > 0:
        omega = omega.conjugate()
    return x0, omega


def integrate_profile(
    kappa_abs: float,
    overtone: int = 0,
    tau: int = 1,
    order: int = 3,
    x_min: float = 2.02,
    x_max: float = 35.0,
    points: int = 2500,
) -> dict[str, np.ndarray | complex | float]:
    """Integra Z_tau su una griglia uniforme in x_* e decompone in Madelung.

    L'integrazione usa x_* come variabile indipendente (portandosi dietro x)
    cosi' la griglia e' equispaziata: A''/A puo' allora essere ricalcolato per
    differenze finite come controllo indipendente della chiusura algebrica.
    """
    if not 2.0 < x_min < x_max:
        raise ValueError("servono 2 < x_min < x_max")
    if points < 50:
        raise ValueError("points deve essere almeno 50")

    _, omega = dirac_qnm_wkb(kappa_abs, overtone, tau, order)
    omega_sq = omega**2
    xstar_min = float(tortoise(x_min))
    xstar_max = float(tortoise(x_max))
    xstar_eval = np.linspace(xstar_min, xstar_max, points)

    # Ingoing all'orizzonte per exp(-i omega t): Z ~ exp(-i Omega x_*).
    z0 = np.exp(-1j * omega * xstar_min)
    chi0 = -1j * omega * z0

    def rhs(_xstar: float, y: np.ndarray) -> np.ndarray:
        z, chi, x = y
        f = 1.0 - 2.0 / x.real
        potential = dirac_potential(x.real, kappa_abs, tau)
        return np.array([chi, (potential - omega_sq) * z, f], dtype=complex)

    solution = solve_ivp(
        rhs,
        (xstar_min, xstar_max),
        np.array([z0, chi0, x_min], dtype=complex),
        t_eval=xstar_eval,
        rtol=2.0e-10,
        atol=2.0e-12,
        method="DOP853",
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    z = solution.y[0]
    chi = solution.y[1]
    x_eval = solution.y[2].real

    epsilon = 1.0 / kappa_abs
    amplitude = np.abs(z)
    h = np.asarray(eikonal_h(x_eval))
    h_prime = np.asarray(eikonal_h_prime(x_eval))

    # y=Z'/Z=A'/A+i S'/eps, con ' rispetto a x_*.  Il momento riscalato e'
    # P=eps*Im(Z'/Z), e la parte reale dell'ODE riscalata da' A''/A.
    logarithmic_derivative = chi / z
    momentum = epsilon * np.imag(logarithmic_derivative)
    spin_term = tau * epsilon * h_prime
    eikonal_term = h**2
    omega_scaled_sq = omega_sq * epsilon**2
    target = np.full_like(x_eval, float(np.real(omega_scaled_sq)))

    # Chiusura algebrica: Q_M = Re(Omega^2) - h^2 - tau eps h' - P^2.
    q_madelung = target - eikonal_term - spin_term - momentum**2
    # Stima indipendente: Q_M = -eps^2 A''/A per differenze finite.
    step = xstar_eval[1] - xstar_eval[0]
    amplitude_second = np.gradient(np.gradient(amplitude, step), step)
    q_madelung_fd = -(epsilon**2) * amplitude_second / amplitude

    interior = slice(2, -2)
    closure_residual = float(
        np.max(np.abs(q_madelung + eikonal_term + spin_term + momentum**2 - target))
    )
    fd_residual = float(
        np.max(np.abs(q_madelung[interior] - q_madelung_fd[interior]))
    )

    return {
        "x": x_eval,
        "xstar": xstar_eval,
        "z": z,
        "amplitude": amplitude,
        "potential": np.asarray(dirac_potential(x_eval, kappa_abs, tau)),
        "h": h,
        "h_prime": h_prime,
        "momentum": momentum,
        "eikonal_term": eikonal_term,
        "spin_term": spin_term,
        "q_madelung": q_madelung,
        "q_madelung_fd": q_madelung_fd,
        "target": target,
        "omega": omega,
        "omega_scaled": omega * epsilon,
        "epsilon": epsilon,
        "kappa_abs": float(kappa_abs),
        "tau": int(tau),
        "closure_residual": closure_residual,
        "fd_residual": fd_residual,
    }


def hierarchy_scaling(
    kappa_values: tuple[float, ...] = (2.0, 4.0, 8.0, 16.0),
    tau: int = 1,
    overtone: int = 0,
    region: str = "far",
    far_window: tuple[float, float] = (20.0, 50.0),
    peak_window: float = 4.0,
    x_min: float = 2.0001,
    x_max: float = 60.0,
    points: int = 8000,
) -> list[dict[str, float]]:
    """Ampiezza dei termini O(eps) e O(eps^2) al variare di K.

    `region="far"` misura in una finestra radiale lontana dai turning point,
    dove |P| resta lontano da zero e la serie locale e' valida: li' il termine
    di spin va come eps e Q_M come eps^2.

    `region="peak"` misura attorno al massimo di barriera, dove per un QNM i
    due turning point coalescono: li' Q_M *non* segue eps^2 e la gerarchia
    locale si rompe.  Le due misure insieme sono la diagnosi, non una sola.

    Nota numerica: `x_min` deve stare vicino all'orizzonte.  La condizione
    ingoing Z~exp(-i Omega x_*) e' esatta solo per V->0, e il ramo riflesso
    che si introduce a x_min troppo grande contamina Q_M con un termine che
    non scala in eps -- effetto visibile come pendenza < 2 a K alto.
    """
    if region not in ("far", "peak"):
        raise ValueError("region deve essere 'far' oppure 'peak'")
    rows: list[dict[str, float]] = []
    for kappa_abs in kappa_values:
        x_peak = dirac_potential_peak(kappa_abs, tau)
        data = integrate_profile(kappa_abs, overtone, tau, 3, x_min, x_max, points)
        x = np.asarray(data["x"])
        if region == "far":
            mask = (x > far_window[0]) & (x < far_window[1])
        else:
            xstar = np.asarray(data["xstar"])
            mask = np.abs(xstar - float(tortoise(x_peak))) <= peak_window * 0.5
        if not mask.any():
            raise RuntimeError("finestra di misura vuota")
        momentum = np.asarray(data["momentum"])[mask]
        spin = float(np.max(np.abs(np.asarray(data["spin_term"])[mask])))
        madelung = float(np.max(np.abs(np.asarray(data["q_madelung"])[mask])))
        rows.append(
            {
                "K": float(kappa_abs),
                "epsilon": float(data["epsilon"]),
                "x_peak": x_peak,
                "momentum_min": float(np.min(np.abs(momentum))),
                "spin_max": spin,
                "q_madelung_max": madelung,
                "ratio": spin / madelung if madelung else float("nan"),
                "fd_residual": float(data["fd_residual"]),
            }
        )
    return rows


def scaling_exponents(rows: list[dict[str, float]]) -> dict[str, float]:
    """Pendenze log-log in eps=1/K dei due termini della gerarchia.

    Attese lontano dai turning point: 1 per lo spin, 2 per il Madelung.
    """
    log_eps = np.log([row["epsilon"] for row in rows])
    return {
        "spin_exponent": float(np.polyfit(log_eps, np.log([r["spin_max"] for r in rows]), 1)[0]),
        "madelung_exponent": float(
            np.polyfit(log_eps, np.log([r["q_madelung_max"] for r in rows]), 1)[0]
        ),
    }


def save_csv(data: dict[str, np.ndarray | complex | float], path: Path) -> None:
    z = np.asarray(data["z"])
    table = np.column_stack(
        [
            data["x"],
            data["xstar"],
            z.real,
            z.imag,
            data["amplitude"],
            data["potential"],
            data["momentum"],
            data["eikonal_term"],
            data["spin_term"],
            data["q_madelung"],
            data["q_madelung_fd"],
            data["target"],
        ]
    )
    header = (
        "x,xstar,Re_Z,Im_Z,A,V_tau,P,h2,spin_term,Q_M,Q_M_fd,Re_Omega2"
    )
    np.savetxt(path, table, delimiter=",", header=header, comments="")


def save_plot(data: dict[str, np.ndarray | complex | float], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(path.parent / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xstar = np.asarray(data["xstar"])
    kappa_abs = data["kappa_abs"]
    tau = data["tau"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)

    axes[0].plot(xstar, data["potential"], label=rf"$M^2V_\tau$, $K={kappa_abs:g}$, $\tau={tau:+d}$")
    axes[0].set_ylabel(r"$M^2 V_\tau$")
    axes[0].legend()

    amplitude = np.asarray(data["amplitude"])
    axes[1].semilogy(xstar, amplitude / amplitude.max(), label=r"$A/\max A$")
    axes[1].set_ylabel("ampiezza normalizzata")
    axes[1].legend()

    axes[2].plot(xstar, np.asarray(data["momentum"]) ** 2, label=r"$P^2$")
    axes[2].plot(xstar, data["eikonal_term"], label=r"$h^2$")
    axes[2].plot(xstar, data["spin_term"], label=r"$\tau\varepsilon h'$  ($O(\varepsilon)$)")
    axes[2].plot(xstar, data["q_madelung"], label=r"$Q_M$  ($O(\varepsilon^2)$)")
    axes[2].plot(xstar, data["target"], "k--", linewidth=1.2, label=r"$\mathrm{Re}\,\Omega^2$")
    axes[2].set_xlabel(r"$x_*=r_*/M$")
    axes[2].set_ylabel("termini riscalati")
    axes[2].legend(ncol=2)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kappa", type=float, default=2.0, help="K=|kappa|, intero per i modi fisici")
    parser.add_argument("--tau", type=int, choices=(1, -1), default=1, help="tau=sigma*sign(kappa)")
    parser.add_argument("--n", type=int, default=0, dest="overtone")
    parser.add_argument("--order", type=int, choices=(1, 3), default=3)
    parser.add_argument("--x-min", type=float, default=2.0001)
    parser.add_argument("--x-max", type=float, default=60.0)
    parser.add_argument("--points", type=int, default=2500)
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--scaling", action="store_true", help="stampa la tabella di scaling in K")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = integrate_profile(
        args.kappa, args.overtone, args.tau, args.order, args.x_min, args.x_max, args.points
    )
    sign = "p" if args.tau > 0 else "m"
    stem = f"dirac_madelung_K{args.kappa:g}_tau{sign}_n{args.overtone}_wkb{args.order}"
    csv_path = args.output_dir / f"{stem}.csv"
    save_csv(data, csv_path)
    if not args.no_plot:
        save_plot(data, args.output_dir / f"{stem}.png")

    omega = complex(data["omega"])
    print(f"K = {args.kappa:g}, tau = {args.tau:+d}, eps = {data['epsilon']:.6f}")
    print(f"M omega  = {omega.real:.12f} {omega.imag:+.12f} i")
    print(f"Omega=eps*M omega = {complex(data['omega_scaled']).real:.12f} "
          f"{complex(data['omega_scaled']).imag:+.12f} i")
    print(f"CSV: {csv_path}")
    if not args.no_plot:
        print(f"PNG: {args.output_dir / f'{stem}.png'}")
    print(f"residuo di chiusura        = {data['closure_residual']:.3e}")
    print(f"max |Q_M - (-eps^2 A''/A)| = {data['fd_residual']:.3e}")

    if args.scaling:
        for region, label in (("far", "lontano dai turning point"), ("peak", "al massimo di barriera")):
            rows = hierarchy_scaling(tau=args.tau, overtone=args.overtone, region=region)
            exponents = scaling_exponents(rows)
            print()
            print(f"--- {label} (region={region}) ---")
            print("  K      eps      min|P|     max|tau eps h'|      max|Q_M|     rapporto   res.FD")
            for row in rows:
                print(
                    f"{row['K']:5.1f}  {row['epsilon']:7.4f}  {row['momentum_min']:9.6f}  "
                    f"{row['spin_max']:16.9e}  {row['q_madelung_max']:12.6e}  "
                    f"{row['ratio']:9.4f}  {row['fd_residual']:.1e}"
                )
            print(
                f"pendenza in eps:  spin = {exponents['spin_exponent']:.4f} (attesa 1), "
                f"Madelung = {exponents['madelung_exponent']:.4f} (attesa 2)"
            )


if __name__ == "__main__":
    main()
