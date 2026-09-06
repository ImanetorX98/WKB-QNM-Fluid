#!/usr/bin/env python3
"""Symbolic check of the projective WKB recursion for massless Dirac on Kerr.

After the eikonal scaling, the radial Chandrasekhar system can be written

    epsilon P_-' + i kappa(x) P_- = w(x) P_+
    epsilon P_+' - i kappa(x) P_+ = w(x) P_-

and z=P_+/P_- obeys

    epsilon z' = w + 2 i kappa z - w z**2.

This script constructs z0,z1,z2 and checks the Riccati equation coefficient
by coefficient.  The two signs label the local momentum branches.
"""

from __future__ import annotations

import sympy as sp


x = sp.symbols("x", real=True)
epsilon = sp.symbols("epsilon")
tau = sp.symbols("tau", real=True, nonzero=True)
kappa = sp.Function("kappa")(x)
w = sp.Function("w")(x)
k = sp.sqrt(kappa**2 - w**2)

z0 = sp.I * (kappa + tau * k) / w
denominator = -2 * sp.I * tau * k
z1 = sp.diff(z0, x) / denominator
z2 = (sp.diff(z1, x) + w * z1**2) / denominator


def branch_reduce(expr: sp.Expr) -> sp.Expr:
    """Use tau**2=1 and simplify the chosen WKB branch."""

    return sp.factor(sp.simplify(expr.subs(tau**2, 1)))


z_trial = z0 + epsilon * z1 + epsilon**2 * z2
residual = sp.expand(
    epsilon * sp.diff(z_trial, x)
    - w
    - 2 * sp.I * kappa * z_trial
    + w * z_trial**2
)


def coefficient(order: int) -> sp.Expr:
    return branch_reduce(sp.expand(residual).coeff(epsilon, order))


def main() -> None:
    for order in range(3):
        value = coefficient(order)
        if value != 0:
            raise AssertionError(f"Riccati coefficient epsilon^{order}: {value}")
        print(f"[ok] Riccati coefficient epsilon^{order}")

    logarithmic_derivative = sp.expand(
        (w * z_trial - sp.I * kappa) / epsilon
    )
    leading = branch_reduce(
        sp.limit(epsilon * logarithmic_derivative, epsilon, 0)
        - sp.I * tau * k
    )
    if leading != 0:
        raise AssertionError(f"Leading phase momentum: {leading}")
    print("[ok] (log P_-)' leading term = i tau k / epsilon")

    print("\nProjective hierarchy")
    print("z0 =", z0)
    print("z1 = z0'/(-2 i tau k)")
    print("z2 = (z1' + w z1**2)/(-2 i tau k)")


if __name__ == "__main__":
    main()
