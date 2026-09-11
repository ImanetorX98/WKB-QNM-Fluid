#!/usr/bin/env python3
"""La crescita secolare r^2 della correzione adiabatica e' il tempo ritardato.

PROBLEMA (Codex, 11 settembre, `research/vaidya_forced_boundary_check_2026-09-11.md`
§5).  Scrivendo la correzione forzata come H = G y, l'equazione esterna da'
y ~ k r^2 con k = 2 i omega, quindi Mdot H/G cresce come Mdot |omega| r^2: **lo
sviluppo adiabatico non e' uniformemente piccolo a grande r**, e i limiti
adiabatico e di grande distanza non commutano.  Verificato qui in modo
indipendente: il coefficiente quadratico e' esattamente k.

DIAGNOSI.  Non e' una patologia del problema: e' un **termine secolare**, cioe'
l'espansione di Taylor di un argomento ritardato.  Nella regione esterna
l'equazione caratteristica di Vaidya si riduce a

    2 psi_vr + psi_rr = 0,

e si verifica per sostituzione diretta che **ogni** funzione della coordinata
uscente u = v - 2 r_* la risolve esattamente, qualunque sia la dipendenza di M
da v.  Percio' nella regione esterna l'adiabaticita' non e' un'approssimazione:
e' esatta, purche' la si organizzi attorno a u.

Sviluppare F(v - 2r) attorno a v a r fissato produce

    F(v-2r) = F(v) - 2r F'(v) + 2 r^2 F''(v) - ...

e il termine quadratico e' esattamente il secolare osservato.  Il conto sotto
mostra che anche il **coefficiente** coincide: con F = exp(-i int omega dv) e
omega proporzionale a 1/M si ottiene, nella parte lineare in Mdot,

    [H/G]_{r^2} = -2 i omega'(M) = +2 i omega = k .

CONSEGUENZA.  Il difetto sta nella scelta della variabile, non nell'espansione.
L'ansatz adiabatico va congelato al **tempo di emissione** v - 2r_*, non al tempo
locale v: l'onda uscente osservata a raggio r e' partita dalla fotosfera quando la
massa era diversa, e congelarla a v accumula un errore di fase quadratico nella
distanza percorsa.  Con quella variabile il secolare non c'e'.

Uso: python3.13 calculations/vaidya_retarded_uniformity.py
"""

from __future__ import annotations

import sympy as sp


def outgoing_characteristic_is_exact() -> bool:
    """Ogni F(v - 2 r) risolve 2 psi_vr + psi_rr = 0.  Nessuna ipotesi su M(v)."""
    v, r = sp.symbols("v r")
    F = sp.Function("F")
    psi = F(v - 2 * r)
    residual = sp.simplify(2 * sp.diff(psi, v, r) + sp.diff(psi, r, 2))
    return residual == 0


def secular_coefficient_from_source() -> sp.Expr:
    """Coefficiente di r^2 in H/G, dall'equazione forzata (via Codex)."""
    r, omega, transport = sp.symbols("r omega K")
    k = 2 * sp.I * omega
    mode = sp.exp(k * r)
    source = 2 * sp.diff(r * sp.diff(mode, r), r) - 2 * transport * sp.diff(mode, r)
    forcing = sp.expand(sp.simplify(source / mode))

    quad, lin = sp.symbols("a b")
    trial = quad * r**2 + lin * r
    operator = sp.diff(trial, r, 2) + 2 * sp.I * omega * sp.diff(trial, r)
    poly = sp.Poly(sp.expand(operator - forcing), r)
    return sp.simplify(sp.solve(poly.coeff_monomial(r), quad)[0])


def secular_coefficient_from_retardation() -> sp.Expr:
    """Coefficiente di r^2 nell'espansione di Taylor del ritardo, a ordine Mdot.

    psi = F(v - 2r),  F = exp(-i int omega(M(v)) dv),  omega = omega_1/M.
    Il termine 2 r^2 F''(v) contiene, lineare in Mdot, il pezzo -2 i omega'(M) Mdot.
    """
    v, r, mass, rate, omega1 = sp.symbols("v r M Mdot omega_1", positive=True)
    epsilon = sp.symbols("epsilon")
    mass_of_v = mass + rate * v
    omega_of_v = omega1 / mass_of_v
    phase = sp.integrate(omega_of_v, v)
    F = sp.exp(-sp.I * phase)

    # termine quadratico dello sviluppo di F(v-2r) attorno a v
    quadratic = sp.diff(F, v, 2) * (2 * r**2) / F
    # parte lineare in Mdot, valutata a v=0
    linear = sp.series(quadratic.subs(v, 0).rewrite(sp.exp), rate, 0, 2).removeO()
    coefficient = sp.simplify(sp.expand(linear).coeff(rate) / r**2)
    return sp.simplify(coefficient.subs(mass, 1)), epsilon


def exact_cancellation_in_far_zone():
    """Sottrazione di ritardo che azzera la sorgente forzata, zona esterna.

    Con chi = a r^2 + b r, l'annullamento simultaneo del termine lineare e della
    costante da' a = 2 i omega e b = -2K, e il residuo e' **esattamente zero**:
    a ordine Mdot la correzione forzata nella zona esterna e' interamente
    ritardo, non una correzione dinamica.
    """
    r, omega, transport = sp.symbols("r omega K")
    quad, lin = sp.symbols("a b")
    mode = sp.exp(2 * sp.I * omega * r)
    operator = lambda X: sp.diff(X, r, 2) - 2 * sp.I * omega * sp.diff(X, r)
    local = sp.expand(sp.simplify(
        (2 * sp.diff(r * sp.diff(mode, r), r) - 2 * transport * sp.diff(mode, r)) / mode))
    chi = quad * r**2 + lin * r
    retarded = sp.expand(sp.simplify(local - operator(mode * chi) / mode))
    poly = sp.Poly(retarded, r)
    solution = sp.solve([poly.coeff_monomial(r), poly.coeff_monomial(1)],
                        [quad, lin], dict=True)[0]
    return solution[quad], solution[lin], sp.simplify(retarded.subs(solution))


def main() -> None:
    print("1. Ogni F(v - 2 r_*) risolve esattamente l'equazione esterna?")
    print(f"   {outgoing_characteristic_is_exact()}")
    print("   Nessuna ipotesi su M(v): nella regione esterna l'adiabaticita'")
    print("   organizzata attorno a u = v - 2 r_* e' ESATTA, non approssimata.")
    print()

    from_source = secular_coefficient_from_source()
    print("2. Coefficiente di r^2 in H/G dall'equazione forzata:")
    print(f"   a = {from_source}")
    print(f"   coincide con k = 2 i omega: {sp.simplify(from_source - 2*sp.I*sp.Symbol('omega')) == 0}")
    print()

    from_retard, _ = secular_coefficient_from_retardation()
    omega1 = sp.Symbol("omega_1", positive=True)
    print("3. Coefficiente di r^2 dallo sviluppo del ritardo, lineare in Mdot:")
    print(f"   a = {from_retard}")
    print(f"   con omega = omega_1/M valutata a M=1, cioe' omega = omega_1:")
    print(f"   confronto con 2 i omega: {sp.simplify(from_retard - 2*sp.I*omega1)}")
    print()
    print("Se la differenza e' zero, il secolare e' interamente il ritardo,")
    print("e la cura e' congelare a v - 2 r_* invece che a v.")
    print()
    quad, lin, residue = exact_cancellation_in_far_zone()
    print("4. Sottrazione completa nella zona esterna, chi = a r^2 + b r:")
    print(f"   a = {quad}   b = {lin}")
    print(f"   residuo della sorgente forzata: {residue}")
    print("   Zero: a ordine Mdot la correzione esterna e' TUTTA ritardo.")
    print("   Il quadratico e' la deriva di massa, il lineare lo shift di frequenza,")
    print("   entrambi valutati al tempo di emissione.")


if __name__ == "__main__":
    main()
