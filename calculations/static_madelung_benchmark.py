#!/usr/bin/env python3
"""Benchmark statico del diagnostico di Madelung contro l'errore WKB.

La soluzione radiale usata per costruire il diagnostico non e' la soluzione
WKB: la frequenza viene da Leaver e il dato iniziale ingoing viene dalla sua
serie di Frobenius.  Il profilo e' poi integrato attraverso la barriera.

Con L=ell+1/2 ed epsilon=1/L si usa

    P = epsilon Im(psi'/psi),
    Q_M = -epsilon^2 A''/A,

e si misura il candidato

    E_M = int w |Q_M| dx_* /
          int w (|epsilon omega|^2 + epsilon^2 |V| + |P|^2) dx_*.

La finestra w e' gaussiana, centrata al massimo del potenziale, con scala
sigma=sqrt(V0/(-V0'')); entrambe le derivate sono rispetto a x_*=r_*/M.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from leaver_qnm import coefficients, leaver_qnm  # noqa: E402
from schwarzschild_wkb import (  # noqa: E402
    _derivative_functions,
    potential_peak,
    qnm_wkb,
    schwarzschild_potential,
)


@dataclass(frozen=True)
class BenchmarkRow:
    spin: int
    ell: int
    overtone: int
    epsilon: float
    omega_re: float
    omega_im: float
    wkb1_relative_error: float
    wkb3_relative_error: float
    madelung_indicator: float
    madelung_peak: float
    fd_relative_residual: float
    leaver_depth_shift: float


def tortoise(x: float | np.ndarray) -> float | np.ndarray:
    value = np.asarray(x)
    result = value + 2.0 * np.log(value / 2.0 - 1.0)
    return float(result) if result.ndim == 0 else result


def inverse_tortoise(xstar: float | np.ndarray) -> float | np.ndarray:
    value = np.asarray(xstar)
    result = 2.0 * (1.0 + lambertw(np.exp(value / 2.0 - 1.0)).real)
    return float(result) if result.ndim == 0 else result


def frobenius_log_derivative(
    x: float, omega: complex, ell: int, spin: int, terms: int = 120
) -> complex:
    """Restituisce psi'/psi in x_* dalla serie di Leaver presso l'orizzonte."""

    series = np.empty(terms, dtype=complex)
    series[0] = 1.0
    alpha0, beta0, _ = coefficients(0, omega, ell, spin)
    series[1] = -beta0 / alpha0
    for n in range(1, terms - 1):
        alpha, beta, gamma = coefficients(n, omega, ell, spin)
        series[n + 1] = -(beta * series[n] + gamma * series[n - 1]) / alpha

    z = 1.0 - 2.0 / x
    powers = z ** np.arange(terms)
    radial_sum = np.dot(series, powers)
    radial_sum_z = np.dot(np.arange(1, terms) * series[1:], powers[:-1])
    dz_dr = 2.0 / x**2
    log_prefactor_r = 1j * omega * (
        1.0 - 2.0 / (x - 2.0) + 4.0 / x
    )
    log_derivative_r = log_prefactor_r + radial_sum_z * dz_dr / radial_sum
    return (1.0 - 2.0 / x) * log_derivative_r


def barrier_profile(
    ell: int,
    overtone: int,
    spin: int,
    omega: complex,
    points: int = 4001,
    half_widths: float = 2.25,
    window_scale: float = 1.0,
) -> dict[str, np.ndarray | float]:
    x_peak = potential_peak(ell, spin)
    derivatives = _derivative_functions(ell, spin)
    v0 = float(derivatives[0](x_peak))
    v2 = float(derivatives[2](x_peak))
    sigma = float(np.sqrt(v0 / (-v2)))
    center = float(tortoise(x_peak))
    xstar = np.linspace(
        center - half_widths * sigma,
        center + half_widths * sigma,
        points,
    )
    x_left = float(inverse_tortoise(xstar[0]))
    initial_log_derivative = frobenius_log_derivative(
        x_left, omega, ell, spin
    )

    def rhs(_xstar: float, state: np.ndarray) -> np.ndarray:
        psi, chi, radius = state
        radius_real = float(radius.real)
        f = 1.0 - 2.0 / radius_real
        potential = schwarzschild_potential(radius_real, ell, spin)
        return np.array(
            [chi, (potential - omega**2) * psi, f], dtype=complex
        )

    solution = solve_ivp(
        rhs,
        (float(xstar[0]), float(xstar[-1])),
        np.array([1.0 + 0.0j, initial_log_derivative, x_left], dtype=complex),
        t_eval=xstar,
        method="DOP853",
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    psi = solution.y[0]
    chi = solution.y[1]
    radius = solution.y[2].real
    potential = np.asarray(schwarzschild_potential(radius, ell, spin))
    epsilon = 1.0 / (ell + 0.5)
    logarithmic_derivative = chi / psi
    momentum = epsilon * np.imag(logarithmic_derivative)

    # Conseguenza della Riccati esatta, non della WKB di barriera.
    amplitude_curvature = (
        potential - float(np.real(omega**2))
        + np.imag(logarithmic_derivative) ** 2
    )
    q_madelung = -(epsilon**2) * amplitude_curvature

    # Controllo indipendente della curvatura tramite differenze finite su |psi|.
    amplitude = np.abs(psi)
    step = float(xstar[1] - xstar[0])
    amplitude_second_fd = np.gradient(np.gradient(amplitude, step), step)
    q_madelung_fd = -(epsilon**2) * amplitude_second_fd / amplitude
    interior = slice(4, -4)
    fd_scale = max(float(np.max(np.abs(q_madelung[interior]))), 1.0e-15)
    fd_relative_residual = float(
        np.max(np.abs(q_madelung_fd[interior] - q_madelung[interior]))
        / fd_scale
    )

    if window_scale <= 0.0:
        raise ValueError("window_scale deve essere positivo")
    weight = np.exp(-0.5 * ((xstar - center) / (window_scale * sigma)) ** 2)
    omega_scale = abs(epsilon * omega) ** 2
    denominator_density = (
        omega_scale + epsilon**2 * np.abs(potential) + np.abs(momentum) ** 2
    )
    numerator = np.trapezoid(weight * np.abs(q_madelung), xstar)
    denominator = np.trapezoid(weight * denominator_density, xstar)

    return {
        "xstar": xstar,
        "radius": radius,
        "q_madelung": q_madelung,
        "momentum": momentum,
        "indicator": float(numerator / denominator),
        "peak": float(np.max(np.abs(q_madelung))),
        "fd_relative_residual": fd_relative_residual,
    }


def benchmark_mode(
    ell: int,
    overtone: int,
    spin: int,
    depth: int = 3000,
    window_scale: float = 1.0,
) -> BenchmarkRow:
    wkb1 = qnm_wkb(ell, overtone, spin, order=1).omega_M
    wkb3 = qnm_wkb(ell, overtone, spin, order=3).omega_M
    exact = leaver_qnm(
        ell, overtone, spin, guess=wkb3, depth=depth, tolerance=2.0e-13
    )
    shallower = leaver_qnm(
        ell,
        overtone,
        spin,
        guess=exact,
        depth=max(500, depth // 2),
        tolerance=2.0e-13,
    )
    profile = barrier_profile(
        ell, overtone, spin, exact, window_scale=window_scale
    )
    return BenchmarkRow(
        spin=spin,
        ell=ell,
        overtone=overtone,
        epsilon=float(1.0 / (ell + 0.5)),
        omega_re=float(exact.real),
        omega_im=float(exact.imag),
        wkb1_relative_error=float(abs(wkb1 - exact) / abs(exact)),
        wkb3_relative_error=float(abs(wkb3 - exact) / abs(exact)),
        madelung_indicator=float(profile["indicator"]),
        madelung_peak=float(profile["peak"]),
        fd_relative_residual=float(profile["fd_relative_residual"]),
        leaver_depth_shift=float(abs(shallower - exact) / abs(exact)),
    )


def correlation(
    rows: list[BenchmarkRow], predictor_name: str, error_name: str
) -> dict[str, float]:
    diagnostic = np.log10([getattr(row, predictor_name) for row in rows])
    errors = np.log10([getattr(row, error_name) for row in rows])
    pearson = pearsonr(diagnostic, errors)
    spearman = spearmanr(diagnostic, errors)
    return {
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def controlled_comparison(
    rows: list[BenchmarkRow], error_name: str
) -> dict[str, float]:
    """Confronta E_M con epsilon^2, condizionando separatamente su ogni n."""

    diagnostic = np.log10([row.madelung_indicator for row in rows])
    baseline = np.log10([row.epsilon**2 for row in rows])
    errors = np.log10([getattr(row, error_name) for row in rows])
    overtones = np.array([row.overtone for row in rows])
    diagnostic_residual = np.empty_like(diagnostic)
    error_residual = np.empty_like(errors)
    for overtone in np.unique(overtones):
        group = overtones == overtone
        diagnostic_residual[group] = diagnostic[group] - np.polyval(
            np.polyfit(baseline[group], diagnostic[group], 1), baseline[group]
        )
        error_residual[group] = errors[group] - np.polyval(
            np.polyfit(baseline[group], errors[group], 1), baseline[group]
        )
    partial = pearsonr(diagnostic_residual, error_residual)

    def loocv_rmse(predictor: np.ndarray) -> float:
        predictions = np.empty_like(errors)
        unique_overtones = np.unique(overtones)
        overtone_columns = np.column_stack(
            [overtones == value for value in unique_overtones[1:]]
        ) if len(unique_overtones) > 1 else np.empty((len(rows), 0))
        design = np.column_stack(
            [np.ones(len(rows)), predictor, overtone_columns]
        ).astype(float)
        for index in range(len(errors)):
            train = np.arange(len(errors)) != index
            coefficients_fit = np.linalg.lstsq(
                design[train], errors[train], rcond=None
            )[0]
            predictions[index] = design[index] @ coefficients_fit
        return float(np.sqrt(np.mean((errors - predictions) ** 2)))

    return {
        "partial_r": float(partial.statistic),
        "partial_p": float(partial.pvalue),
        "madelung_loocv_rmse_dex": loocv_rmse(diagnostic),
        "epsilon2_loocv_rmse_dex": loocv_rmse(baseline),
    }


def default_modes(
    max_ell: int, overtones: tuple[int, ...]
) -> list[tuple[int, int, int]]:
    return [
        (spin, ell, overtone)
        for overtone in overtones
        for spin in (0, 1, 2)
        for ell in range(max(1, spin), max_ell + 1)
    ]


def save_csv(rows: list[BenchmarkRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-ell", type=int, default=8)
    parser.add_argument(
        "--n",
        type=int,
        action="append",
        dest="overtones",
        help="overtone; ripetere l'opzione per una scansione combinata",
    )
    parser.add_argument("--depth", type=int, default=3000)
    parser.add_argument(
        "--window-scale",
        type=float,
        default=1.0,
        help="larghezza gaussiana in unita' di sqrt(V0/-V0'')",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "data" / "static_madelung_benchmark.csv",
    )
    args = parser.parse_args()

    overtones = tuple(args.overtones) if args.overtones else (0,)
    modes = default_modes(args.max_ell, overtones)
    rows = [
        benchmark_mode(ell, n, spin, args.depth, args.window_scale)
        for spin, ell, n in modes
    ]
    save_csv(rows, args.output)

    print(" s ell n       Re(Mw)       Im(Mw)       E_M       err.WKB1    err.WKB3")
    for row in rows:
        print(
            f" {row.spin:d} {row.ell:3d} {row.overtone:1d}  "
            f"{row.omega_re:11.7f} {row.omega_im:+11.7f}  "
            f"{row.madelung_indicator:9.3e}  "
            f"{row.wkb1_relative_error:9.3e}  {row.wkb3_relative_error:9.3e}"
        )
    print()
    for error_name in ("wkb1_relative_error", "wkb3_relative_error"):
        stats = correlation(rows, "madelung_indicator", error_name)
        baseline_stats = correlation(rows, "epsilon", error_name)
        controlled = controlled_comparison(rows, error_name)
        print(
            f"log(E_M) vs log({error_name}): "
            f"Pearson r={stats['pearson_r']:+.4f} (p={stats['pearson_p']:.2e}), "
            f"Spearman rho={stats['spearman_rho']:+.4f} "
            f"(p={stats['spearman_p']:.2e})"
        )
        print(
            f"  controllo epsilon: r(E_M,error | epsilon^2)="
            f"{controlled['partial_r']:+.4f} (p={controlled['partial_p']:.2e}); "
            f"r(epsilon,error)={baseline_stats['pearson_r']:+.4f}"
        )
        print(
            f"  LOOCV RMSE: E_M={controlled['madelung_loocv_rmse_dex']:.4f} dex, "
            f"epsilon^2={controlled['epsilon2_loocv_rmse_dex']:.4f} dex"
        )
    print(
        f"max residuale FD relativo = {max(row.fd_relative_residual for row in rows):.3e}"
    )
    print(
        f"max shift Leaver per profondita = {max(row.leaver_depth_shift for row in rows):.3e}"
    )
    print(f"CSV: {args.output}")


if __name__ == "__main__":
    main()
