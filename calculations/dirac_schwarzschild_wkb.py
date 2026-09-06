#!/usr/bin/env python3
"""Symbolic checks for massless Dirac WKB–Madelung on Schwarzschild.

Conventions
-----------
f = 1 - 2 M/r,  dx/dr = 1/f,  D_x = f D_r.
h = sqrt(f)/r and W = kappa h.

The supersymmetric partner potentials are

    V_sigma = W**2 + sigma D_x W,   sigma = +/-1.

For the eikonal scaling K=|kappa|, epsilon=1/K, s=sign(kappa),
omega=K Omega, W=K s h, the scalar component equations read

    epsilon**2 Z_sigma,xx
      + [Omega**2 - h**2 - sigma*s*epsilon*h_x] Z_sigma = 0.

For a real oscillatory solution written as

    Z = P**(-1/2) exp(i integral(P dx)/epsilon),

the exact scalar Madelung closure is

    q = P**2 - epsilon**2 (P**(-1/2))_xx/P**(-1/2)
      = P**2 + Q_M.

This script verifies the partner potentials, the first phase-momentum
coefficients, and the large-K peak expansion.
"""

from __future__ import annotations

import sympy as sp


r, M = sp.symbols("r M", positive=True, finite=True)
Omega = sp.symbols("Omega", finite=True)
K = sp.symbols("K", positive=True, finite=True)
epsilon = sp.symbols("epsilon", positive=True, finite=True)
sigma, sgn_kappa = sp.symbols("sigma s", real=True)

f = 1 - 2 * M / r
Delta = r * (r - 2 * M)


def Dx(expr: sp.Expr) -> sp.Expr:
    """Derivative with respect to Schwarzschild tortoise coordinate x=r*."""

    return sp.factor(f * sp.diff(expr, r))


h = sp.sqrt(f) / r
h_x = sp.factor(Dx(h))
h_x_expected = sp.sqrt(f) * (3 * M - r) / r**3

V_plus = sp.factor(K**2 * h**2 + K * h_x)
V_minus = sp.factor(K**2 * h**2 - K * h_x)
V_cho = sp.factor(K * sp.sqrt(Delta) / r**4 * (K * sp.sqrt(Delta) - (r - 3 * M)))

# Eikonal scalar equation: q_sigma = q0 + epsilon q1.
q0 = sp.factor(Omega**2 - h**2)
q1 = sp.factor(-sigma * sgn_kappa * h_x)
q0_x = sp.factor(Dx(q0))
q0_xx = sp.factor(Dx(q0_x))

# Scalar Madelung coefficient Q_M = epsilon**2 Q2 + O(epsilon**3).
Q2 = sp.factor(q0_xx / (4 * q0) - 5 * q0_x**2 / (16 * q0**2))

# P = P0 + epsilon P1 + epsilon**2 P2 + ... from
# P**2 = q0 + epsilon*q1 + epsilon**2 F[P],
# F[P] = (P**(-1/2))_xx/P**(-1/2).
P0 = sp.sqrt(q0)
P1 = sp.factor(q1 / (2 * P0))
F0 = sp.factor(3 * (Dx(P0) / P0) ** 2 / 4 - Dx(Dx(P0)) / (2 * P0))
P2 = sp.factor((F0 - P1**2) / (2 * P0))

# The same coefficient written only in q0 and q1.
P2_expected = sp.factor(
    5 * q0_x**2 / (32 * q0 ** sp.Rational(5, 2))
    - q0_xx / (8 * q0 ** sp.Rational(3, 2))
    - q1**2 / (8 * q0 ** sp.Rational(3, 2))
)

# The next Madelung coefficient is odd because the Dirac subprincipal
# symbol q1 is already O(epsilon).  It is the first mixed spin-gradient term.
P_trial = P0 + epsilon * P1 + epsilon**2 * P2
F_trial = sp.factor(
    3 * (Dx(P_trial) / P_trial) ** 2 / 4
    - Dx(Dx(P_trial)) / (2 * P_trial)
)
F1 = sp.factor(sp.diff(F_trial, epsilon).subs(epsilon, 0))
Q3 = sp.factor(-F1)
q1_x = sp.factor(Dx(q1))
q1_xx = sp.factor(Dx(q1_x))
Q3_expected = sp.factor(
    q1_xx / (4 * q0)
    - 5 * q0_x * q1_x / (8 * q0**2)
    - q1 * q0_xx / (4 * q0**2)
    + 5 * q1 * q0_x**2 / (8 * q0**3)
)

# Photon-sphere / potential-peak data.
r_ph = 3 * M
h0 = sp.simplify(h.subs(r, r_ph))
h2_xx_0 = sp.simplify(Dx(Dx(h**2)).subs(r, r_ph))
h_xx_0 = sp.simplify(Dx(h_x).subs(r, r_ph))

# Let r_peak = 3 M + c/K.  Vanishing dV_sigma/dr at leading order fixes c.
c = sp.symbols("c", real=True)
z = sp.symbols("z", positive=True)
r_trial = r_ph + c * z
V_sigma_z = h**2 / z**2 + sigma * sgn_kappa * h_x / z
dV_dr_trial = sp.diff(V_sigma_z, r).subs(r, r_trial)
leading_peak_equation = sp.simplify(sp.limit(z * dV_dr_trial, z, 0, dir="+"))
c_solution = sp.solve(sp.Eq(leading_peak_equation, 0), c)[0]

V_peak_series = sp.series(V_sigma_z.subs(r, r_trial).subs(c, c_solution), z, 0, 1).removeO()
V_peak_series = sp.factor(V_peak_series.subs(z, 1 / K))

# Leading barrier-top quantities for the usual first-order WKB QNM formula.
V0_leading = sp.factor(K**2 * h0**2)
V0_xx_leading = sp.factor(K**2 * h2_xx_0)
barrier_curvature = sp.simplify(sp.sqrt(-2 * V0_xx_leading))


def check_zero(name: str, expr: sp.Expr) -> None:
    residue = sp.simplify(expr)
    if residue != 0:
        raise AssertionError(f"{name} failed: {residue}")
    print(f"[ok] {name}")


def main() -> None:
    check_zero("D_x h", h_x - h_x_expected)
    check_zero("Cho potential equals V_plus", V_plus - V_cho)
    check_zero("P2 recursion", P2 - P2_expected)
    check_zero("Q2 = -F[P0]", Q2 + F0)
    check_zero("mixed spin-Madelung Q3", Q3 - Q3_expected)

    print("\n--- Schwarzschild Dirac data ---")
    print("f             =", f)
    print("h             =", h)
    print("D_x h         =", h_x)
    print("V_plus        =", V_plus)
    print("V_minus       =", V_minus)

    print("\n--- Eikonal WKB–Madelung hierarchy ---")
    print("q0            =", q0)
    print("q1            =", q1)
    print("D_x q0        =", q0_x)
    print("D_x^2 q0      =", q0_xx)
    print("Q2            =", Q2)
    print("Q3            =", Q3)
    print("P1            =", P1)
    print("P2            =", P2)

    print("\n--- Photon sphere and partner-potential peak ---")
    print("h(3M)         =", h0)
    print("D_x^2(h^2)|3M=", h2_xx_0)
    print("D_x^2 h|3M    =", h_xx_0)
    print("peak equation =", leading_peak_equation)
    print("c in r_peak=3M+c/K:", c_solution)
    print("V_peak through K^0:", V_peak_series)
    print("K^2 h0^2      =", V0_leading)
    print("K^2 (h^2)_xx =", V0_xx_leading)
    print("sqrt(-2V0_xx)=", barrier_curvature)


if __name__ == "__main__":
    main()
