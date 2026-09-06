#!/usr/bin/env python3
"""QNMs di Schwarzschild al primo e terzo ordine WKB (Iyer--Will).

Le frequenze restituite sono adimensionali: Omega = M * omega, con G=c=1.
La convenzione temporale e' exp(-i omega t), quindi un modo smorzato ha
Im(omega) < 0.
"""

from __future__ import annotations

import argparse
import cmath
from dataclasses import dataclass
from functools import lru_cache
from typing import Callable

import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar


@dataclass(frozen=True)
class WKBResult:
    spin: int
    ell: int
    overtone: int
    order: int
    r_peak_over_M: float
    omega_M: complex
    lambda2: float = 0.0
    lambda3: float = 0.0


def validate_mode(spin: int, ell: int, overtone: int) -> None:
    if spin not in (0, 1, 2):
        raise ValueError("spin deve essere 0 (scalare), 1 (EM) o 2 (assiale gravitazionale)")
    ell_min = spin
    if ell < ell_min:
        raise ValueError(f"per spin={spin} serve ell >= {ell_min}")
    if overtone < 0:
        raise ValueError("l'overtone n deve essere >= 0")


def schwarzschild_potential(x: np.ndarray | float, ell: int, spin: int) -> np.ndarray | float:
    """M^2 V_s(x), x=r/M, per campi massless s=0,1 e assiale s=2."""
    x_arr = np.asarray(x)
    f = 1.0 - 2.0 / x_arr
    value = f * (ell * (ell + 1) / x_arr**2 + 2.0 * (1.0 - spin**2) / x_arr**3)
    return float(value) if np.ndim(value) == 0 else value


@lru_cache(maxsize=None)
def _derivative_functions(ell: int, spin: int) -> tuple[Callable[[float], float], ...]:
    """V e derivate d^k V/dx_*^k, k=1,...,6; x_*=r_*/M."""
    x = sp.symbols("x", positive=True)
    f = 1 - 2 / x
    expr = f * (sp.Integer(ell * (ell + 1)) / x**2 + 2 * (1 - spin**2) / x**3)
    expressions = [expr]
    for _ in range(6):
        expressions.append(sp.factor(f * sp.diff(expressions[-1], x)))
    return tuple(sp.lambdify(x, item, modules="numpy") for item in expressions)


def potential_peak(ell: int, spin: int) -> float:
    """Posizione x=r/M del massimo esterno del potenziale."""
    validate_mode(spin, ell, 0)
    result = minimize_scalar(
        lambda x: -schwarzschild_potential(x, ell, spin),
        bounds=(2.0 + 1.0e-8, 50.0),
        method="bounded",
        options={"xatol": 1.0e-14},
    )
    if not result.success:
        raise RuntimeError(f"ricerca del massimo fallita: {result.message}")
    return float(result.x)


def qnm_wkb(ell: int, overtone: int = 0, spin: int = 2, order: int = 3) -> WKBResult:
    """Calcola Omega=M*omega al primo oppure al terzo ordine WKB.

    La forma esplicita usata e'

        Omega^2 = V0 + sqrt(-2 V0'') Lambda2
                  - i alpha sqrt(-2 V0'') (1 + Lambda3),

    con alpha=n+1/2 e derivate rispetto a x_*=r_*/M, valutate al
    massimo di V. Lambda2 e Lambda3 sono spesso chiamate Lambda e Omega
    nella letteratura, rispettivamente.
    """
    validate_mode(spin, ell, overtone)
    if order not in (1, 3):
        raise ValueError("questo modulo implementa order=1 oppure order=3")

    x0 = potential_peak(ell, spin)
    derivatives = _derivative_functions(ell, spin)
    values = np.array([float(fn(x0)) for fn in derivatives], dtype=float)
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

    omega_sq = (
        v0
        + root_curvature * lambda2
        - 1j * alpha * root_curvature * (1.0 + lambda3)
    )
    omega = cmath.sqrt(omega_sq)
    if omega.real < 0:
        omega = -omega

    return WKBResult(
        spin=spin,
        ell=ell,
        overtone=overtone,
        order=order,
        r_peak_over_M=x0,
        omega_M=omega,
        lambda2=float(lambda2),
        lambda3=float(lambda3),
    )


def physical_units(omega_M: complex, mass_solar: float) -> tuple[float, float]:
    """Restituisce frequenza oscillatoria [Hz] e tempo di damping [s]."""
    if mass_solar <= 0:
        raise ValueError("mass_solar deve essere positivo")
    gravitational_second_per_solar_mass = 4.925490947e-6
    time_scale = gravitational_second_per_solar_mass * mass_solar
    frequency_hz = omega_M.real / (2.0 * np.pi * time_scale)
    damping_s = time_scale / abs(omega_M.imag)
    return float(frequency_hz), float(damping_s)


def _format_result(result: WKBResult, mass_solar: float | None) -> str:
    w = result.omega_M
    text = (
        f"s={result.spin}, ell={result.ell}, n={result.overtone}, WKB{result.order}\n"
        f"r_peak/M = {result.r_peak_over_M:.12f}\n"
        f"M omega  = {w.real:.12f} {w.imag:+.12f} i"
    )
    if result.order == 3:
        text += f"\nLambda2 = {result.lambda2:.12e}\nLambda3 = {result.lambda3:.12e}"
    if mass_solar is not None:
        hz, tau = physical_units(w, mass_solar)
        text += f"\nf_osc = {hz:.6f} Hz\ntau   = {tau:.6e} s  (M={mass_solar:g} M_sun)"
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spin", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--ell", type=int, default=2)
    parser.add_argument("--n", type=int, default=0, dest="overtone")
    parser.add_argument("--order", choices=("1", "3", "both"), default="both")
    parser.add_argument("--mass-solar", type=float, default=None)
    args = parser.parse_args()

    orders = (1, 3) if args.order == "both" else (int(args.order),)
    results = [qnm_wkb(args.ell, args.overtone, args.spin, order) for order in orders]
    print("\n\n".join(_format_result(item, args.mass_solar) for item in results))
    if len(results) == 2:
        delta = abs(results[1].omega_M - results[0].omega_M)
        print(f"\n|Omega_WKB3 - Omega_WKB1| = {delta:.6e}")


if __name__ == "__main__":
    main()
