#!/usr/bin/env python3
"""Residui di Madelung risolti per ordine sui QNM di Schwarzschild.

Per la chiusura esatta

    q = u^2 - eps^2 (u^(-1/2))''/u^(-1/2)

si formano u^(0)=u0, u^(1)=u0+eps^2 u2 e
u^(2)=u0+eps^2 u2+eps^4 u4. Il residuo R_N della chiusura misura localmente
cio' che il troncamento non ha ancora assorbito. A differenza del diagnostico
costruito da |psi|, usa soltanto il potenziale e la frequenza di riferimento.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from leaver_qnm import leaver_qnm  # noqa: E402
from schwarzschild_wkb import (  # noqa: E402
    _derivative_functions,
    potential_peak,
    qnm_wkb,
)


@dataclass(frozen=True)
class RemainderRow:
    spin: int
    ell: int
    overtone: int
    epsilon: float
    omega_re: float
    omega_im: float
    wkb1_relative_error: float
    wkb3_relative_error: float
    remainder_0: float
    remainder_1: float
    remainder_2: float


def tortoise(x: float | np.ndarray) -> float | np.ndarray:
    value = np.asarray(x)
    result = value + 2.0 * np.log(value / 2.0 - 1.0)
    return float(result) if result.ndim == 0 else result


def _total_derivative(expression: sp.Expr, jets: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.Add(*[
        sp.diff(expression, jets[index]) * jets[index + 1]
        for index in range(len(jets) - 1)
    ])


@lru_cache(maxsize=1)
def _remainder_functions() -> tuple[object, ...]:
    """Compila i residui usando il getto q,...,q^(6) rispetto a x_*."""

    jets = sp.symbols("q0:7")
    epsilon = sp.symbols("epsilon", positive=True)
    q0, q1, q2, q3, q4, _, _ = jets
    u0 = sp.sqrt(q0)
    u2 = (5 * q1**2 - 4 * q0 * q2) / (32 * q0 ** sp.Rational(5, 2))
    u4 = (
        64 * q0**3 * q4
        - 448 * q0**2 * q1 * q3
        - 304 * q0**2 * q2**2
        + 1768 * q0 * q1**2 * q2
        - 1105 * q1**4
    ) / (2048 * q0 ** sp.Rational(11, 2))
    truncations = (
        u0,
        u0 + epsilon**2 * u2,
        u0 + epsilon**2 * u2 + epsilon**4 * u4,
    )
    compiled: list[object] = []
    for momentum in truncations:
        momentum_prime = _total_derivative(momentum, jets)
        momentum_second = _total_derivative(momentum_prime, jets)
        amplitude_curvature = (
            sp.Rational(3, 4) * (momentum_prime / momentum) ** 2
            - sp.Rational(1, 2) * momentum_second / momentum
        )
        remainder = q0 - (momentum**2 - epsilon**2 * amplitude_curvature)
        normalization = sp.Abs(q0) + sp.Abs(momentum) ** 2 + epsilon**2 * sp.Abs(
            amplitude_curvature
        )
        compiled.append(
            sp.lambdify(
                (epsilon, *jets), (remainder, normalization), modules="numpy"
            )
        )
    return tuple(compiled)


def integrated_remainders(
    ell: int,
    spin: int,
    omega: complex,
    points: int = 4001,
    half_widths: float = 2.25,
    window_scale: float = 1.0,
) -> tuple[float, float, float]:
    derivatives = _derivative_functions(ell, spin)
    x_peak = potential_peak(ell, spin)
    v0 = float(derivatives[0](x_peak))
    v2 = float(derivatives[2](x_peak))
    sigma = float(np.sqrt(v0 / (-v2)))
    center = float(tortoise(x_peak))

    # Si integra in x_* ma le derivate del potenziale sono gia' funzioni di r.
    xstar = np.linspace(
        center - half_widths * sigma,
        center + half_widths * sigma,
        points,
    )
    # Inversione stabile per il ristretto intervallo di barriera.
    from scipy.special import lambertw

    radius = 2.0 * (1.0 + lambertw(np.exp(xstar / 2.0 - 1.0)).real)
    epsilon = 1.0 / (ell + 0.5)
    potential_jets = [np.asarray(function(radius), dtype=float) for function in derivatives]
    q_jets: list[np.ndarray] = [
        epsilon**2 * (omega**2 - potential_jets[0])
    ]
    q_jets.extend(-epsilon**2 * item for item in potential_jets[1:])
    weight = np.exp(
        -0.5 * ((xstar - center) / (window_scale * sigma)) ** 2
    )

    indicators: list[float] = []
    for function in _remainder_functions():
        remainder, normalization = function(epsilon, *q_jets)
        numerator = np.trapezoid(weight * np.abs(remainder), xstar)
        denominator = np.trapezoid(weight * np.asarray(normalization, dtype=float), xstar)
        indicators.append(float(numerator / denominator))
    return tuple(indicators)  # type: ignore[return-value]


def benchmark_mode(
    ell: int, overtone: int, spin: int, depth: int = 3000
) -> RemainderRow:
    wkb1 = qnm_wkb(ell, overtone, spin, order=1).omega_M
    wkb3 = qnm_wkb(ell, overtone, spin, order=3).omega_M
    exact = leaver_qnm(
        ell, overtone, spin, guess=wkb3, depth=depth, tolerance=2.0e-13
    )
    remainders = integrated_remainders(ell, spin, exact)
    return RemainderRow(
        spin=spin,
        ell=ell,
        overtone=overtone,
        epsilon=float(1.0 / (ell + 0.5)),
        omega_re=float(exact.real),
        omega_im=float(exact.imag),
        wkb1_relative_error=float(abs(wkb1 - exact) / abs(exact)),
        wkb3_relative_error=float(abs(wkb3 - exact) / abs(exact)),
        remainder_0=remainders[0],
        remainder_1=remainders[1],
        remainder_2=remainders[2],
    )


def modes(max_ell: int, overtones: tuple[int, ...]) -> list[tuple[int, int, int]]:
    return [
        (spin, ell, overtone)
        for overtone in overtones
        for spin in (0, 1, 2)
        for ell in range(max(1, spin), max_ell + 1)
    ]


def statistics(
    rows: list[RemainderRow], predictor: str, error: str
) -> dict[str, float]:
    x = np.log10([getattr(row, predictor) for row in rows])
    y = np.log10([getattr(row, error) for row in rows])
    pearson = pearsonr(x, y)
    spearman = spearmanr(x, y)
    return {
        "pearson": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def save_csv(rows: list[RemainderRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-ell", type=int, default=8)
    parser.add_argument("--n", type=int, action="append", dest="overtones")
    parser.add_argument("--depth", type=int, default=3000)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "data" / "order_resolved_benchmark.csv",
    )
    args = parser.parse_args()
    overtones = tuple(args.overtones) if args.overtones else (0,)
    rows = [
        benchmark_mode(ell, overtone, spin, args.depth)
        for spin, ell, overtone in modes(args.max_ell, overtones)
    ]
    save_csv(rows, args.output)

    print(" s ell n       R0          R1          R2       err.WKB1   err.WKB3")
    for row in rows:
        print(
            f" {row.spin:d} {row.ell:3d} {row.overtone:1d}  "
            f"{row.remainder_0:10.3e} {row.remainder_1:10.3e} "
            f"{row.remainder_2:10.3e} {row.wkb1_relative_error:10.3e} "
            f"{row.wkb3_relative_error:10.3e}"
        )
    print()
    for error in ("wkb1_relative_error", "wkb3_relative_error"):
        for predictor in ("remainder_0", "remainder_1", "remainder_2"):
            result = statistics(rows, predictor, error)
            print(
                f"{predictor} vs {error}: Pearson={result['pearson']:+.4f}, "
                f"Spearman={result['spearman']:+.4f}"
            )
        print()
    print(f"CSV: {args.output}")


if __name__ == "__main__":
    main()
