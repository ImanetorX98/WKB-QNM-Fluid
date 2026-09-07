#!/usr/bin/env python3
"""Controlli simbolici per il ponte scalare Schwarzschild -> Kerr.

Parte dall'equazione radiale di Klein-Gordon separata in Boyer-Lindquist,
rimuove la derivata prima con Psi=sqrt(r^2+a^2) R e verifica il termine
lineare slow-Kerr.  Le convenzioni sono exp(-i omega t+i m phi).
"""

from __future__ import annotations

import sympy as sp


r = sp.symbols("r", positive=True)
mass = sp.symbols("M", positive=True)
a = sp.symbols("a", real=True)
omega = sp.symbols("omega")
m = sp.symbols("m", integer=True)
separation = sp.symbols("lambda")
ell_term = sp.symbols("L2", positive=True)

delta = r**2 - 2 * mass * r + a**2
h_squared = r**2 + a**2
h = sp.sqrt(h_squared)
radial_lapse = delta / h_squared
frame_rate = a / h_squared
kerr_k = h_squared * omega - a * m


def dstar(expression: sp.Expr) -> sp.Expr:
    """Derivata rispetto a r_* con dr_*/dr=(r^2+a^2)/Delta."""

    return sp.factor(radial_lapse * sp.diff(expression, r))


geometric_potential = sp.factor(dstar(dstar(h)) / h)
radial_q = sp.factor(
    (omega - m * frame_rate) ** 2
    - delta * separation / h_squared**2
    - geometric_potential
)


def assert_zero(label: str, expression: sp.Expr) -> None:
    reduced = sp.factor(sp.together(expression))
    if reduced != 0:
        raise AssertionError(f"{label}: {reduced}")
    print(f"[ok] {label}")


def check_exact_schrodinger_reduction() -> None:
    psi = sp.Function("Psi")(r)
    radial_r = psi / h
    teukolsky_radial = sp.diff(delta * sp.diff(radial_r, r), r) + (
        kerr_k**2 / delta - separation
    ) * radial_r
    schrodinger_radial = dstar(dstar(psi)) + radial_q * psi
    assert_zero(
        "exact scalar Kerr radial reduction",
        schrodinger_radial - delta * teukolsky_radial / h**3,
    )


def check_schwarzschild_limit() -> None:
    f = 1 - 2 * mass / r
    expected = omega**2 - f * (ell_term / r**2 + 2 * mass / r**3)
    assert_zero(
        "Schwarzschild scalar potential",
        radial_q.subs({a: 0, separation: ell_term}) - expected,
    )


def check_slow_kerr_limit() -> None:
    # A_lm(a omega)=ell(ell+1)+O(a^2), hence
    # lambda=A_lm+a^2 omega^2-2 a m omega.
    slow_separation = ell_term - 2 * a * m * omega
    expansion = sp.series(
        radial_q.subs(separation, slow_separation), a, 0, 2
    ).removeO()
    f = 1 - 2 * mass / r
    expected = (
        omega**2
        - f * (ell_term / r**2 + 2 * mass / r**3)
        - 4 * a * m * mass * omega / r**3
    )
    assert_zero("linear slow-Kerr frame dragging", expansion - expected)


def check_horizon_symbol() -> None:
    r_plus = sp.symbols("r_plus", positive=True)
    omega_horizon = a / (r_plus**2 + a**2)
    at_horizon = radial_q.subs(r, r_plus)
    # On Delta(r_+)=0 both the separation and rescaling potentials vanish.
    numerator = sp.together(
        at_horizon - (omega - m * omega_horizon) ** 2
    ).subs(mass, (r_plus**2 + a**2) / (2 * r_plus))
    assert_zero("co-rotating horizon symbol", numerator)


def check_slow_kerr_open_source() -> None:
    omega_re, omega_im = sp.symbols("omega_R omega_I", real=True)
    slow_q = (
        (omega_re + sp.I * omega_im) ** 2
        - (1 - 2 * mass / r) * (ell_term / r**2 + 2 * mass / r**3)
        - 4 * a * m * mass * (omega_re + sp.I * omega_im) / r**3
    )
    expected_imaginary = 2 * omega_im * (
        omega_re - 2 * a * m * mass / r**3
    )
    assert_zero(
        "slow-Kerr radial open-flow source",
        sp.im(sp.expand_complex(slow_q)) - expected_imaginary,
    )


def main() -> None:
    check_exact_schrodinger_reduction()
    check_schwarzschild_limit()
    check_slow_kerr_limit()
    check_horizon_symbol()
    check_slow_kerr_open_source()


if __name__ == "__main__":
    main()
