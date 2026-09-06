#!/usr/bin/env python3
"""Third-order barrier WKB check for massless Dirac QNMs on Schwarzschild.

This reproduces the classic Cho (2003) values for the V_+ partner.  Units M=1
are used internally.  All potential derivatives are iterated tortoise
derivatives D_x = f D_r; replacing D_x^n by f^n D_r^n would be incorrect.
"""

from __future__ import annotations

import cmath

import sympy as sp


r, Ksym = sp.symbols("r K", positive=True)
f = 1 - 2 / r
Delta = r * (r - 2)
V = Ksym * sp.sqrt(Delta) / r**4 * (Ksym * sp.sqrt(Delta) - (r - 3))


def Dx(expr: sp.Expr) -> sp.Expr:
    return sp.factor(f * sp.diff(expr, r))


Vx = [V]
for _ in range(6):
    Vx.append(Dx(Vx[-1]))


def peak_radius(K: int) -> float:
    """Physical root of Cho's quartic, mapped to r/M."""

    y = sp.symbols("y", real=True)
    polynomial = 21 * y**4 - 12 * K * y**3 - 14 * y**2 + 4 * K * y + 1
    roots = [complex(z) for z in sp.nroots(polynomial, n=40, maxsteps=200)]
    physical = [z.real for z in roots if abs(z.imag) < 1e-25 and 0 < z.real < 1]
    candidates = [2 / (1 - yy**2) for yy in physical]
    outside = [rr for rr in candidates if rr > 2]
    if not outside:
        raise RuntimeError(f"No physical peak root for K={K}: {roots}")
    # The barrier maximum is the root in the interval (2,3).
    return min(outside, key=lambda rr: abs(rr - 2.7))


def value(expr: sp.Expr, K: int, r0: float) -> float:
    return float(sp.N(expr.subs({Ksym: K, r: r0}), 40))


def wkb3_frequency(K: int, overtone: int = 0) -> tuple[float, complex]:
    r0 = peak_radius(K)
    derivatives = [value(expr, K, r0) for expr in Vx]
    V0, _, V2, V3, V4, V5, V6 = derivatives
    alpha = overtone + 0.5
    root = (-2 * V2) ** 0.5

    Lambda = (
        (V4 / V2) * (0.25 + alpha**2) / 8
        - (V3 / V2) ** 2 * (7 + 60 * alpha**2) / 288
    ) / root

    Omega_corr = (
        5 * (V3 / V2) ** 4 * (77 + 188 * alpha**2) / 6912
        - (V3**2 * V4 / V2**3) * (51 + 100 * alpha**2) / 384
        + (V4 / V2) ** 2 * (67 + 68 * alpha**2) / 2304
        + (V3 * V5 / V2**2) * (19 + 28 * alpha**2) / 288
        - (V6 / V2) * (5 + 4 * alpha**2) / 288
    ) / (-2 * V2)

    omega_sq = V0 + root * Lambda - 1j * alpha * root * (1 + Omega_corr)
    omega = cmath.sqrt(omega_sq)
    if omega.real < 0:
        omega = -omega
    if omega.imag > 0:
        omega = omega.conjugate()
    return r0, omega


def main() -> None:
    reference = {
        1: complex(0.176, -0.100),
        2: complex(0.379, -0.0965),
        3: complex(0.574, -0.0963),
        4: complex(0.767, -0.0963),
        5: complex(0.960, -0.0963),
    }
    print(" K       r_peak/M               M omega_WKB3              |delta Cho|")
    for kval in range(1, 6):
        r0, omega = wkb3_frequency(kval)
        err = abs(omega - reference[kval])
        print(f"{kval:2d}   {r0:13.9f}   {omega.real: .9f}{omega.imag:+.9f}i   {err:.3e}")


if __name__ == "__main__":
    main()
