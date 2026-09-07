#!/usr/bin/env python3
"""Test fuori Schwarzschild su una barriera di Poschl--Teller esattamente risolta.

Per V(y)=L^2 sech^2(y), con alpha=1, le frequenze QNM esatte sono

    omega = sqrt(L^2-1/4) - i (n+1/2).

La funzione d'onda QNM e' un ipergeometrico terminante. Questo permette di
calcolare il diagnostico di Madelung senza usare Leaver e di verificare se la
correlazione osservata in Schwarzschild sopravvive a una famiglia indipendente.
"""

from __future__ import annotations

import argparse
import cmath
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class PoschlTellerRow:
    barrier_scale: float
    overtone: int
    epsilon: float
    omega_re: float
    omega_im: float
    wkb1_relative_error: float
    wkb3_relative_error: float
    madelung_indicator: float
    fd_relative_residual: float
    nodal: bool


def exact_frequency(barrier_scale: float, overtone: int) -> complex:
    if barrier_scale <= 0.5:
        raise ValueError("serve L>1/2 per una parte reale oscillatoria")
    return complex(
        np.sqrt(barrier_scale**2 - 0.25), -(overtone + 0.5)
    )


def wkb_frequency(barrier_scale: float, overtone: int, order: int) -> complex:
    """WKB Iyer--Will per V=L^2 sech^2(y), alpha=1."""

    if order not in (1, 3):
        raise ValueError("order deve essere 1 oppure 3")
    height = barrier_scale**2
    alpha_n = overtone + 0.5
    v0 = height
    v2 = -2.0 * height
    v3 = 0.0
    v4 = 16.0 * height
    v5 = 0.0
    v6 = -272.0 * height
    root_curvature = np.sqrt(-2.0 * v2)
    lambda2 = 0.0
    lambda3 = 0.0
    if order == 3:
        lambda2 = (
            0.125 * (v4 / v2) * (0.25 + alpha_n**2)
            - (v3 / v2) ** 2 * (7.0 + 60.0 * alpha_n**2) / 288.0
        ) / root_curvature
        lambda3 = 1.0 / (-2.0 * v2) * (
            5.0 * (v3 / v2) ** 4 * (77.0 + 188.0 * alpha_n**2) / 6912.0
            - (v3**2 * v4 / v2**3) * (51.0 + 100.0 * alpha_n**2) / 384.0
            + (v4 / v2) ** 2 * (67.0 + 68.0 * alpha_n**2) / 2304.0
            + (v3 * v5 / v2**2) * (19.0 + 28.0 * alpha_n**2) / 288.0
            - (v6 / v2) * (5.0 + 4.0 * alpha_n**2) / 288.0
        )
    omega_squared = (
        v0
        + root_curvature * lambda2
        - 1j * alpha_n * root_curvature * (1.0 + lambda3)
    )
    result = cmath.sqrt(omega_squared)
    return result if result.real >= 0.0 else -result


def terminating_hypergeometric(
    overtone: int, b: complex, c: complex, z: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """2F1(-n,b;c;z) e derivata in z come polinomi finiti."""

    result = np.ones_like(z, dtype=complex)
    derivative = np.zeros_like(z, dtype=complex)
    coefficient = 1.0 + 0.0j
    power = np.ones_like(z, dtype=float)
    for index in range(1, overtone + 1):
        coefficient *= (
            (-overtone + index - 1) * (b + index - 1)
            / ((c + index - 1) * index)
        )
        power = power * z
        result = result + coefficient * power
        derivative = derivative + index * coefficient * z ** (index - 1)
    return result, derivative


def exact_wavefunction(
    y: np.ndarray, barrier_scale: float, overtone: int
) -> tuple[np.ndarray, np.ndarray]:
    """Restituisce psi e la derivata logaritmica esatta d psi/dy / psi."""

    omega = exact_frequency(barrier_scale, overtone)
    spectral_lambda = np.sqrt(barrier_scale**2 - 0.25)
    z = 0.5 * (1.0 + np.tanh(y))
    b = -overtone - 2j * spectral_lambda
    c = 0.5 - overtone - 1j * spectral_lambda
    polynomial, polynomial_z = terminating_hypergeometric(overtone, b, c, z)
    prefactor = np.exp(
        (-0.5j * omega) * (np.log(z) + np.log1p(-z))
    )
    psi = prefactor * polynomial
    dz_dy = 2.0 * z * (1.0 - z)
    log_derivative = (
        1j * omega * np.tanh(y) + dz_dy * polynomial_z / polynomial
    )
    return psi, log_derivative


def madelung_indicator(
    barrier_scale: float,
    overtone: int,
    points: int = 4001,
    window_scale: float = 1.0,
) -> tuple[float, float, bool]:
    # sigma=sqrt(V0/-V0'')=1/sqrt(2) in the y coordinate.
    sigma = 1.0 / np.sqrt(2.0)
    y = np.linspace(-2.25 * sigma, 2.25 * sigma, points)
    step = float(y[1] - y[0])
    psi, logarithmic_derivative = exact_wavefunction(y, barrier_scale, overtone)
    amplitude = np.abs(psi)
    nodal = bool(np.min(amplitude) / np.max(amplitude) < 1.0e-10)
    if nodal:
        return float("nan"), float("nan"), True
    epsilon = 1.0 / barrier_scale
    momentum = epsilon * np.imag(logarithmic_derivative)
    normalized_potential = 1.0 / np.cosh(y) ** 2
    omega_normalized = exact_frequency(barrier_scale, overtone) / barrier_scale
    q_madelung = (
        float(np.real(omega_normalized**2))
        - normalized_potential
        - momentum**2
    )

    amplitude_second = np.full_like(amplitude, np.nan)
    amplitude_second[2:-2] = (
        -amplitude[4:]
        + 16.0 * amplitude[3:-1]
        - 30.0 * amplitude[2:-2]
        + 16.0 * amplitude[1:-3]
        - amplitude[:-4]
    ) / (12.0 * step**2)
    q_madelung_fd = -(epsilon**2) * amplitude_second / amplitude
    interior = slice(4, -4)
    scale = max(float(np.max(np.abs(q_madelung[interior]))), 1.0e-15)
    fd_residual = float(
        np.max(np.abs(q_madelung_fd[interior] - q_madelung[interior])) / scale
    )

    weight = np.exp(-0.5 * (y / (window_scale * sigma)) ** 2)
    denominator_density = (
        abs(omega_normalized) ** 2 + normalized_potential + momentum**2
    )
    numerator = np.trapezoid(weight * np.abs(q_madelung), y)
    denominator = np.trapezoid(weight * denominator_density, y)
    return float(numerator / denominator), fd_residual, False


def benchmark_row(barrier_scale: float, overtone: int) -> PoschlTellerRow:
    exact = exact_frequency(barrier_scale, overtone)
    wkb1 = wkb_frequency(barrier_scale, overtone, 1)
    wkb3 = wkb_frequency(barrier_scale, overtone, 3)
    indicator, fd_residual, nodal = madelung_indicator(barrier_scale, overtone)
    return PoschlTellerRow(
        barrier_scale=float(barrier_scale),
        overtone=overtone,
        epsilon=float(1.0 / barrier_scale),
        omega_re=float(exact.real),
        omega_im=float(exact.imag),
        wkb1_relative_error=float(abs(wkb1 - exact) / abs(exact)),
        wkb3_relative_error=float(abs(wkb3 - exact) / abs(exact)),
        madelung_indicator=indicator,
        fd_relative_residual=fd_residual,
        nodal=nodal,
    )


def correlation(rows: list[PoschlTellerRow], error: str) -> tuple[float, float]:
    regular = [row for row in rows if not row.nodal]
    if len(regular) < 3:
        return float("nan"), float("nan")
    x = np.log10([row.madelung_indicator for row in regular])
    y = np.log10([getattr(row, error) for row in regular])
    return float(pearsonr(x, y).statistic), float(spearmanr(x, y).statistic)


def save_csv(rows: list[PoschlTellerRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, action="append", dest="overtones")
    parser.add_argument(
        "--scales",
        default="1.5,2,3,4,6,8,12,16",
        help="valori separati da virgole di sqrt(V0)/alpha",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "data" / "poschl_teller_benchmark.csv",
    )
    args = parser.parse_args()
    overtones = tuple(args.overtones) if args.overtones else (0,)
    scales = tuple(float(item) for item in args.scales.split(","))
    rows = [benchmark_row(scale, overtone) for overtone in overtones for scale in scales]
    save_csv(rows, args.output)

    print("   L n       E_M       err.WKB1    err.WKB3    residuale FD")
    for row in rows:
        indicator = "   nodale" if row.nodal else f"{row.madelung_indicator:10.3e}"
        residual = "       n/a" if row.nodal else f"{row.fd_relative_residual:10.3e}"
        print(
            f"{row.barrier_scale:5.1f} {row.overtone:1d}  "
            f"{indicator} {row.wkb1_relative_error:10.3e} "
            f"{row.wkb3_relative_error:10.3e} {residual}"
        )
    print()
    for overtone in overtones:
        subset = [row for row in rows if row.overtone == overtone]
        for error in ("wkb1_relative_error", "wkb3_relative_error"):
            pearson, spearman = correlation(subset, error)
            print(
                f"n={overtone}, E_M vs {error}: "
                f"Pearson={pearson:+.4f}, Spearman={spearman:+.4f}"
            )
    print(f"CSV: {args.output}")


if __name__ == "__main__":
    main()
