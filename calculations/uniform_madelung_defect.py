#!/usr/bin/env python3
"""Difetto di Madelung rispetto alla forma uniforme parabolico-cilindrica.

Il test usa la barriera esatta di Poschl--Teller. Per l'equazione normalizzata

    eps^2 psi'' + [omega_hat^2 - sech^2(y)] psi = 0

la forma quadratica al massimo e' q_PC=q0+y^2. La sua soluzione locale e' una
combinazione di D_nu(zeta) e D_nu(-zeta), con

    zeta = exp(-i pi/4) sqrt(2/eps) y,
    nu   = i q0/(2 eps) - 1/2.

Si confrontano le curvature di ampiezza esatta e uniforme senza espandere in
1/q, evitando la singolarita' artificiale dei residui WKB locali.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from calculations.poschl_teller_madelung_benchmark import (
    exact_frequency,
    exact_wavefunction,
    wkb_frequency,
)


@dataclass(frozen=True)
class UniformDefectRow:
    barrier_scale: float
    overtone: int
    epsilon: float
    uniform_defect: float
    full_madelung: float
    defect_over_full: float
    wkb1_relative_error: float
    wkb3_relative_error: float
    nodal: bool


def parabolic_cylinder_log_derivative(
    y: np.ndarray,
    barrier_scale: float,
    overtone: int,
    initial_log_derivative: complex = 0.0j,
) -> np.ndarray:
    """Derivata logaritmica della soluzione quadratica con dati al massimo.

    L'integrazione e' equivalente a scegliere la combinazione di funzioni
    parabolico-cilindriche che soddisfa il dato di Cauchy a y=0. Per i modi
    pari della barriera simmetrica il dato e' psi(0)=1, psi'(0)=0.
    """

    epsilon = 1.0 / barrier_scale
    omega_normalized = exact_frequency(barrier_scale, overtone) / barrier_scale
    q0 = omega_normalized**2 - 1.0
    center = int(np.argmin(np.abs(y)))
    if abs(y[center]) > 1.0e-12:
        raise ValueError("la griglia uniforme deve contenere y=0")

    def rhs(position: float, state: np.ndarray) -> np.ndarray:
        psi, derivative = state
        q_uniform = q0 + position**2
        return np.array(
            [derivative, -q_uniform * psi / epsilon**2], dtype=complex
        )

    initial = np.array([1.0 + 0.0j, initial_log_derivative], dtype=complex)
    right = solve_ivp(
        rhs,
        (0.0, float(y[-1])),
        initial,
        t_eval=y[center:],
        method="DOP853",
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    left_grid = y[: center + 1][::-1]
    left = solve_ivp(
        rhs,
        (0.0, float(y[0])),
        initial,
        t_eval=left_grid,
        method="DOP853",
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    if not right.success or not left.success:
        raise RuntimeError(right.message if not right.success else left.message)
    psi = np.concatenate((left.y[0][::-1][:-1], right.y[0]))
    derivative = np.concatenate((left.y[1][::-1][:-1], right.y[1]))
    return derivative / psi


def uniform_defect_indicator(
    barrier_scale: float,
    overtone: int,
    points: int = 801,
    window_scale: float = 1.0,
) -> tuple[float, float, bool]:
    sigma = 1.0 / np.sqrt(2.0)
    y = np.linspace(-2.25 * sigma, 2.25 * sigma, points)
    exact_psi, exact_log_derivative = exact_wavefunction(
        y, barrier_scale, overtone
    )
    amplitude = np.abs(exact_psi)
    nodal = bool(np.min(amplitude) / np.max(amplitude) < 1.0e-10)
    if nodal:
        return float("nan"), float("nan"), True

    epsilon = 1.0 / barrier_scale
    omega_normalized = exact_frequency(barrier_scale, overtone) / barrier_scale
    q_exact = omega_normalized**2 - 1.0 / np.cosh(y) ** 2
    q0 = omega_normalized**2 - 1.0
    q_uniform = q0 + y**2

    exact_momentum = epsilon * np.imag(exact_log_derivative)
    exact_q_madelung = np.real(q_exact) - exact_momentum**2
    center = int(np.argmin(np.abs(y)))
    uniform_log_derivative = parabolic_cylinder_log_derivative(
        y,
        barrier_scale,
        overtone,
        initial_log_derivative=complex(exact_log_derivative[center]),
    )
    uniform_momentum = epsilon * np.imag(uniform_log_derivative)
    uniform_q_madelung = np.real(q_uniform) - uniform_momentum**2
    defect = exact_q_madelung - uniform_q_madelung

    # La forma parabolico-cilindrica e' uniforme nella regione di coalescenza
    # y=O(sqrt(epsilon)), non su una finestra macroscopica fissata.
    uniform_width = window_scale * sigma * np.sqrt(epsilon)
    weight = np.exp(-0.5 * (y / uniform_width) ** 2)
    normalization_density = (
        abs(omega_normalized) ** 2
        + 1.0 / np.cosh(y) ** 2
        + exact_momentum**2
    )
    normalization = np.trapezoid(weight * normalization_density, y)
    defect_indicator = np.trapezoid(weight * np.abs(defect), y) / normalization
    full_indicator = (
        np.trapezoid(weight * np.abs(exact_q_madelung), y) / normalization
    )
    return float(defect_indicator), float(full_indicator), False


def benchmark_row(barrier_scale: float, overtone: int) -> UniformDefectRow:
    exact = exact_frequency(barrier_scale, overtone)
    wkb1 = wkb_frequency(barrier_scale, overtone, 1)
    wkb3 = wkb_frequency(barrier_scale, overtone, 3)
    defect, full, nodal = uniform_defect_indicator(barrier_scale, overtone)
    return UniformDefectRow(
        barrier_scale=float(barrier_scale),
        overtone=overtone,
        epsilon=float(1.0 / barrier_scale),
        uniform_defect=defect,
        full_madelung=full,
        defect_over_full=float(defect / full) if not nodal else float("nan"),
        wkb1_relative_error=float(abs(wkb1 - exact) / abs(exact)),
        wkb3_relative_error=float(abs(wkb3 - exact) / abs(exact)),
        nodal=nodal,
    )


def correlation(
    rows: list[UniformDefectRow], predictor: str, error: str
) -> tuple[float, float]:
    regular = [row for row in rows if not row.nodal]
    if len(regular) < 3:
        return float("nan"), float("nan")
    x = np.log10([getattr(row, predictor) for row in regular])
    y = np.log10([getattr(row, error) for row in regular])
    return float(pearsonr(x, y).statistic), float(spearmanr(x, y).statistic)


def power_law_exponent(
    rows: list[UniformDefectRow], field: str, minimum_scale: float = 4.0
) -> float:
    """Fit log(field)=constant+p log(L) in the requested asymptotic tail."""

    regular = [
        row
        for row in rows
        if not row.nodal
        and row.barrier_scale >= minimum_scale
        and getattr(row, field) > 0.0
    ]
    if len(regular) < 2:
        return float("nan")
    scales = np.log([row.barrier_scale for row in regular])
    values = np.log([getattr(row, field) for row in regular])
    return float(np.polyfit(scales, values, 1)[0])


def save_csv(rows: list[UniformDefectRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, action="append", dest="overtones")
    parser.add_argument(
        "--scales", default="1.5,2,3,4,6,8,12,16"
    )
    parser.add_argument("--points", type=int, default=801)
    parser.add_argument("--fit-min-scale", type=float, default=4.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "data" / "uniform_madelung_defect.csv",
    )
    args = parser.parse_args()
    overtones = tuple(args.overtones) if args.overtones else (0,)
    scales = tuple(float(item) for item in args.scales.split(","))
    rows = [
        benchmark_row(scale, overtone)
        if args.points == 801
        else _benchmark_with_points(scale, overtone, args.points)
        for overtone in overtones
        for scale in scales
    ]
    save_csv(rows, args.output)

    print("   L n    DeltaQ_unif     Q_full    DeltaQ/Q    err.WKB1   err.WKB3")
    for row in rows:
        if row.nodal:
            print(
                f"{row.barrier_scale:5.1f} {row.overtone:1d}       nodale       nodale       n/a "
                f"{row.wkb1_relative_error:10.3e} {row.wkb3_relative_error:10.3e}"
            )
        else:
            print(
                f"{row.barrier_scale:5.1f} {row.overtone:1d}  "
                f"{row.uniform_defect:12.3e} {row.full_madelung:10.3e} "
                f"{row.defect_over_full:10.3e} {row.wkb1_relative_error:10.3e} "
                f"{row.wkb3_relative_error:10.3e}"
            )
    print()
    for overtone in overtones:
        subset = [row for row in rows if row.overtone == overtone]
        for error in ("wkb1_relative_error", "wkb3_relative_error"):
            pearson, spearman = correlation(subset, "uniform_defect", error)
            print(
                f"n={overtone}, DeltaQ_unif vs {error}: "
                f"Pearson={pearson:+.4f}, Spearman={spearman:+.4f}"
            )
        print(
            f"n={overtone}, esponenti L>={args.fit_min_scale:g}: "
            f"DeltaQ_unif={power_law_exponent(subset, 'uniform_defect', args.fit_min_scale):+.4f}, "
            f"Q_full={power_law_exponent(subset, 'full_madelung', args.fit_min_scale):+.4f}, "
            f"WKB1={power_law_exponent(subset, 'wkb1_relative_error', args.fit_min_scale):+.4f}, "
            f"WKB3={power_law_exponent(subset, 'wkb3_relative_error', args.fit_min_scale):+.4f}"
        )
    print(f"CSV: {args.output}")


def _benchmark_with_points(
    barrier_scale: float, overtone: int, points: int
) -> UniformDefectRow:
    exact = exact_frequency(barrier_scale, overtone)
    wkb1 = wkb_frequency(barrier_scale, overtone, 1)
    wkb3 = wkb_frequency(barrier_scale, overtone, 3)
    defect, full, nodal = uniform_defect_indicator(
        barrier_scale, overtone, points=points
    )
    return UniformDefectRow(
        barrier_scale=float(barrier_scale),
        overtone=overtone,
        epsilon=float(1.0 / barrier_scale),
        uniform_defect=defect,
        full_madelung=full,
        defect_over_full=float(defect / full) if not nodal else float("nan"),
        wkb1_relative_error=float(abs(wkb1 - exact) / abs(exact)),
        wkb3_relative_error=float(abs(wkb3 - exact) / abs(exact)),
        nodal=nodal,
    )


if __name__ == "__main__":
    main()
