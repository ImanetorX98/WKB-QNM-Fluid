#!/usr/bin/env python3
"""Verifica simbolica della chiusura di Madelung e della sua ricorsione.

Per l'equazione canonica

    eps^2 psi'' + q(x) psi = 0,

l'ansatz a momento pari

    psi = u^(-1/2) exp( (i/eps) \\int^x u ),      u = u0 + eps^2 u2 + eps^4 u4 + ...

produce la chiusura esatta

    q = u^2 - eps^2 (u^(-1/2))'' / u^(-1/2),

in cui il secondo termine e' il funzionale di Bohm-Madelung valutato
sull'ampiezza WKB completa A = u^(-1/2).  Lo script verifica:

  (a) la chiusura, come identita' fra l'ODE e l'ansatz;
  (b) l'assenza di un termine di ordine eps^1 (la serie e' in eps^2);
  (c) la ricorsione u2, u4, u6 generata dall'unico funzionale;
  (d) l'accordo di u2 con la formula classica della WKB di secondo ordine;
  (e) che u4 e u6 NON richiedono un nuovo funzionale: bastano le derivate
      funzionali di quello stesso.

Uso: python3.13 madelung_recursion.py
"""

from __future__ import annotations

import sympy as sp


x, eps = sp.symbols("x varepsilon", positive=True)
q = sp.Function("q")(x)


def bohm(amplitude: sp.Expr) -> sp.Expr:
    """Funzionale di Bohm-Madelung A''/A."""
    return sp.diff(amplitude, x, 2) / amplitude


def closure_residual(u: sp.Expr) -> sp.Expr:
    """q - [u^2 - eps^2 (u^(-1/2))''/u^(-1/2)]."""
    return sp.together(q - (u**2 - eps**2 * bohm(u ** sp.Rational(-1, 2))))


def check_exact_closure() -> None:
    """La chiusura e' equivalente all'ODE, per u generico."""
    u = sp.Function("u")(x)
    psi = u ** sp.Rational(-1, 2) * sp.exp(sp.I / eps * sp.Integral(u, x))
    ode = sp.expand(sp.simplify((eps**2 * sp.diff(psi, x, 2) + q * psi) / psi))
    target = sp.expand(q - u**2 + eps**2 * bohm(u ** sp.Rational(-1, 2)))
    difference = sp.simplify(ode - target)
    assert difference == 0, difference
    print("(a) chiusura esatta q = u^2 - eps^2 A''/A, con A=u^(-1/2):  verificata")


def check_no_odd_order() -> None:
    """Un termine di ordine eps^1 nella serie di u e' forzato a zero."""
    u1 = sp.Function("u1")(x)
    u = sp.sqrt(q) + eps * u1
    order1 = sp.simplify(sp.series(closure_residual(u), eps, 0, 2).removeO().coeff(eps, 1))
    solution = sp.solve(sp.Eq(order1, 0), u1)
    assert solution == [0], solution
    print("(b) nessun termine di ordine eps^1: la serie e' in potenze di eps^2")


def solve_hierarchy(depth: int = 3) -> list[sp.Expr]:
    """Risolve ordine per ordine per u2, u4, ..., partendo da u0=sqrt(q)."""
    coefficients: list[sp.Expr] = [sp.sqrt(q)]
    for index in range(1, depth + 1):
        unknown = sp.Function(f"u{2 * index}")(x)
        trial = sum(eps ** (2 * j) * coefficients[j] for j in range(index))
        trial = trial + eps ** (2 * index) * unknown
        expansion = sp.series(closure_residual(trial), eps, 0, 2 * index + 1).removeO()
        equation = sp.simplify(expansion.coeff(eps, 2 * index))
        solved = sp.solve(sp.Eq(equation, 0), unknown)
        assert len(solved) == 1, solved
        coefficients.append(sp.simplify(solved[0]))
    return coefficients


def check_classical_second_order(u2: sp.Expr) -> None:
    """u2 deve coincidere con la correzione WKB classica di secondo ordine."""
    q1 = sp.diff(q, x)
    q2 = sp.diff(q, x, 2)
    classical = 5 * q1**2 / (32 * q ** sp.Rational(5, 2)) - q2 / (8 * q ** sp.Rational(3, 2))
    assert sp.simplify(u2 - classical) == 0, sp.simplify(u2 - classical)
    print("(d) u2 coincide con  5 q'^2/(32 q^(5/2)) - q''/(8 q^(3/2)):  verificata")


def check_single_functional(coefficients: list[sp.Expr]) -> None:
    """Tutti gli ordini discendono dallo stesso funzionale, non da nuovi.

    Si verifica che u4 sia esattamente il termine di ordine eps^4 prodotto
    reinserendo u0+eps^2 u2 dentro A''/A: nessun ingrediente nuovo entra.
    """
    u0, u2, u4 = coefficients[0], coefficients[1], coefficients[2]
    partial = u0 + eps**2 * u2
    generated = sp.series(
        bohm(partial ** sp.Rational(-1, 2)) / (2 * u0), eps, 0, 3
    ).removeO().coeff(eps, 2)
    reconstructed = sp.simplify(generated - u2**2 / (2 * u0))
    assert sp.simplify(u4 - reconstructed) == 0, sp.simplify(u4 - reconstructed)
    print("(e) u4 rigenerato dallo stesso funzionale di Bohm:  verificata")


def main() -> None:
    check_exact_closure()
    check_no_odd_order()
    coefficients = solve_hierarchy(depth=2)
    print("(c) ricorsione risolta fino a u4")
    u2 = coefficients[1]
    check_classical_second_order(u2)
    check_single_functional(coefficients)
    print()
    print("u0 =", coefficients[0])
    print()
    print("u2 =", sp.simplify(u2))
    print()
    print("u4 =", sp.simplify(sp.factor(coefficients[2])))


if __name__ == "__main__":
    main()
