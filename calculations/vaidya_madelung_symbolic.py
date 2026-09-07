#!/usr/bin/env python3
"""Symbolic checks for the scalar Madelung system on ingoing Vaidya.

The script checks three identities used in ``research/covariant_formalism.md``:

1. the spherical reduction of the massless Klein--Gordon equation;
2. the real Hamilton--Jacobi/Madelung equation;
3. the imaginary characteristic continuity equation.

No field equation for M(v) is assumed.  The calculation only uses

    ds^2 = -f(v,r) dv^2 + 2 dv dr + r^2 dOmega^2

and the angular eigenvalue -ell(ell+1).
"""

from __future__ import annotations

import sympy as sp


v, r = sp.symbols("v r", real=True)
epsilon = sp.symbols("epsilon", real=True, nonzero=True)
ell = sp.symbols("ell", integer=True, nonnegative=True)

f = sp.Function("f", real=True)(v, r)
psi = sp.Function("psi", real=True)(v, r)
A = sp.Function("A", positive=True)(v, r)
S = sp.Function("S", real=True)(v, r)


def assert_zero(label: str, expression: sp.Expr) -> None:
    value = sp.simplify(sp.expand(expression))
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    print(f"[ok] {label}")


def check_radial_reduction() -> None:
    """Reduce Box Phi=0 after Phi=psi(v,r)Y_lm/r."""

    phi = psi / r
    # The angular Laplacian has already been replaced by -ell(ell+1).
    box_reduced = (
        sp.diff(r**2 * sp.diff(phi, r), v)
        + sp.diff(r**2 * (sp.diff(phi, v) + f * sp.diff(phi, r)), r)
    ) / r**2 - ell * (ell + 1) * phi / r**2

    potential = sp.diff(f, r) / r + ell * (ell + 1) / r**2
    target = (
        2 * sp.diff(psi, v, r)
        + sp.diff(f * sp.diff(psi, r), r)
        - potential * psi
    )
    assert_zero("Vaidya Klein-Gordon radial reduction", r * box_reduced - target)


def check_madelung_split() -> None:
    """Check real and imaginary parts after psi=A exp(iS/epsilon)."""

    phase = sp.exp(sp.I * S / epsilon)
    field = A * phase
    potential = sp.diff(f, r) / r + ell * (ell + 1) / r**2
    residual = (
        2 * sp.diff(field, v, r)
        + sp.diff(f * sp.diff(field, r), r)
        - potential * field
    )
    reduced = sp.expand(residual / phase)

    amplitude_part = (
        2 * sp.diff(A, v, r)
        + sp.diff(f * sp.diff(A, r), r)
        - potential * A
        - A
        * (
            2 * sp.diff(S, v) * sp.diff(S, r)
            + f * sp.diff(S, r) ** 2
        )
        / epsilon**2
    )
    phase_part = (
        2
        * (
            sp.diff(A, v) * sp.diff(S, r)
            + sp.diff(A, r) * sp.diff(S, v)
            + A * sp.diff(S, v, r)
        )
        + sp.diff(f * A * sp.diff(S, r), r)
        + f * sp.diff(A, r) * sp.diff(S, r)
    ) / epsilon

    # SymPy does not propagate the ``real=True`` assumption from a function to
    # all its unevaluated derivatives.  Checking the complete complex identity
    # is both assumption-free and stronger than calling re()/im() separately.
    assert_zero(
        "complex Madelung split into real and imaginary coefficients",
        reduced - amplitude_part - sp.I * phase_part,
    )

    q_v = -epsilon**2 * (
        2 * sp.diff(A, v, r) + sp.diff(f * sp.diff(A, r), r)
    ) / A
    hamilton_jacobi = (
        2 * sp.diff(S, v) * sp.diff(S, r)
        + f * sp.diff(S, r) ** 2
        + epsilon**2 * potential
        + q_v
    )
    assert_zero(
        "rearranged Hamilton-Jacobi equation",
        hamilton_jacobi + epsilon**2 * amplitude_part / A,
    )

    rho = A**2
    continuity = sp.diff(rho * sp.diff(S, r), v) + sp.diff(
        rho * (sp.diff(S, v) + f * sp.diff(S, r)), r
    )
    assert_zero(
        "characteristic continuity equation",
        continuity - A * epsilon * phase_part,
    )


def check_memory_normalization() -> None:
    """The logarithmic mixed curvature ignores time-only normalizations."""

    normalization = sp.Function("N", positive=True)(v)
    frozen_normalization = sp.Function("N_f", positive=True)(v)
    profile = sp.Function("B", positive=True)(v, r)
    frozen_profile = sp.Function("B_f", positive=True)(v, r)

    normalized_ratio = normalization * profile / (
        frozen_normalization * frozen_profile
    )
    raw_ratio = profile / frozen_profile
    assert_zero(
        "memory diagnostic is invariant under time-only normalizations",
        sp.diff(sp.log(normalized_ratio), v, r)
        - sp.diff(sp.log(raw_ratio), v, r),
    )


def check_kodama_energy_balance() -> None:
    """Check the null-eikonal Kodama-energy drift for f=1-2M(v)/r."""

    mass = sp.Function("M", real=True)(v)
    vaidya_f = 1 - 2 * mass / r
    orbit_metric = sp.Matrix([[-vaidya_f, 1], [1, 0]])
    w_v, w_r = sp.symbols("W_v W_r", real=True)
    wave_vector = sp.Matrix([w_v, w_r])

    # K=partial_v has constant components, so L_K g_ab=partial_v g_ab.
    lie_k_metric = orbit_metric.diff(v)
    drift = -sp.Rational(1, 2) * (
        wave_vector.T * lie_k_metric * wave_vector
    )[0]
    expected = -sp.diff(mass, v) * w_v**2 / r
    assert_zero("Kodama-energy drift in ingoing Vaidya", drift - expected)


def main() -> None:
    check_radial_reduction()
    check_madelung_split()
    check_memory_normalization()
    check_kodama_energy_balance()


if __name__ == "__main__":
    main()
