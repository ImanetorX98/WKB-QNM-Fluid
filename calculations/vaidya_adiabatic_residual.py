#!/usr/bin/env python3
"""L'ansatz adiabatico risolve Vaidya a ordine Mdot?  E cosa resta.

DOMANDA.  L'ipotesi nulla del programma Vaidya e': se l'ampiezza dipende da v
solo attraverso M(v), allora

    Q_x = -2 eps^2 A_vr/A = -2 eps^2 Mdot (d_M d_r A0)/A0 + O(Mdot^2),

cioe' il termine misto e' determinato da Mdot e dal profilo congelato, ed e'
ridondante come lo era E_M.  Va deciso **prima** di scrivere un'evoluzione
caratteristica: e' una domanda all'ordine Mdot, e a quell'ordine si risponde
simbolicamente.

IMPOSTAZIONE.  Vaidya entrante, ds^2 = -f dv^2 + 2 dv dr + r^2 dOmega^2 con
f = 1 - 2M(v)/r.  Per Phi = psi(v,r) Y/r la Klein-Gordon massless si riduce a

    2 psi_vr + d_r(f psi_r) - U psi = 0,      U = f_r/r + l(l+1)/r^2.

Un modo congelato e' la soluzione statica a M fissato, riscritta in coordinate
entranti: con v = t + r_*, il modo e^{-i omega t} R(r) diventa

    psi_ad(v,r) = Z(r; M) exp(-i \\int^v omega(M(v')) dv'),   Z = e^{i omega r_*} R.

Sostituendo e raccogliendo in potenze di Mdot si legge il residuo.

RISULTATO (vedi output): a ordine Mdot^0 si ritrova l'equazione del modo
congelato; a ordine Mdot^1 il residuo NON e' identicamente nullo, e la sua forma
dice se l'ansatz vada corretto e come.

Uso: python3.13 vaidya_adiabatic_residual.py
"""

from __future__ import annotations

import sympy as sp


def build() -> dict[str, sp.Expr]:
    v, r = sp.symbols("v r", real=True)
    ell = sp.symbols("ell", positive=True)
    mdot = sp.symbols("Mdot")  # parametro di conteggio: M(v) = M0 + Mdot*v + ...

    mass = sp.Function("M", positive=True)(v)
    f = 1 - 2 * mass / r
    potential = sp.diff(f, r) / r + ell * (ell + 1) / r**2

    # Modo congelato: Z dipende da r e (parametricamente) da M; la fase accumula
    # omega(M(v)).  Si tiene M(v) simbolica e si espande alla fine.
    amplitude = sp.Function("Z")(r, mass)
    frequency = sp.Function("omega")(mass)
    phase = sp.Function("Phi")(v)
    psi = amplitude * sp.exp(-sp.I * phase)

    equation = 2 * sp.diff(psi, v, r) + sp.diff(f * sp.diff(psi, r), r) - potential * psi
    equation = sp.expand(sp.simplify(equation / sp.exp(-sp.I * phase)))
    # Phi' = omega(M(v))
    equation = equation.subs(sp.Derivative(phase, v), frequency)
    equation = sp.expand(sp.simplify(equation.doit()))
    return {"eq": equation, "v": v, "r": r, "M": mass, "Mdot": mdot,
            "Z": amplitude, "omega": frequency, "f": f, "U": potential}


def split_orders(data: dict[str, sp.Expr]) -> tuple[sp.Expr, sp.Expr]:
    """Separa il residuo in ordine Mdot^0 e Mdot^1."""
    equation = data["eq"]
    v, mass, mdot = data["v"], data["M"], data["Mdot"]
    # M'(v) -> Mdot,  M''(v) -> 0 (tasso costante)
    equation = equation.subs(sp.Derivative(mass, (v, 2)), 0)
    equation = equation.subs(sp.Derivative(mass, v), mdot)
    equation = sp.expand(equation)
    order0 = equation.coeff(mdot, 0)
    order1 = equation.coeff(mdot, 1)
    return sp.simplify(order0), sp.simplify(order1)


def main() -> None:
    data = build()
    order0, order1 = split_orders(data)

    print("=== Ordine Mdot^0: equazione del modo congelato ===")
    print(sp.simplify(order0))
    print()
    print("=== Ordine Mdot^1: residuo dell'ansatz adiabatico ===")
    residual = sp.simplify(order1)
    print(residual)
    print()
    if residual == 0:
        print("Il residuo e' identicamente nullo: l'ansatz adiabatico e' esatto a")
        print("ordine Mdot, e Q_x e' ridondante come previsto dall'ipotesi nulla.")
    else:
        print("Il residuo NON e' identicamente nullo.  L'ansatz adiabatico non")
        print("risolve Vaidya a ordine Mdot: serve una correzione psi_1, e Q_x")
        print("contiene un contributo non determinato dal solo profilo congelato.")
        print()
        print("Struttura del residuo, raccolta nelle derivate di Z:")
        Z, r, M = data["Z"], data["r"], data["M"]
        for derivative in (sp.Derivative(Z, r), sp.Derivative(Z, M),
                           sp.Derivative(Z, r, M), Z):
            coefficient = sp.simplify(residual.coeff(derivative))
            if coefficient != 0:
                print(f"  coefficiente di {derivative}:  {sp.simplify(coefficient)}")


if __name__ == "__main__":
    main()
