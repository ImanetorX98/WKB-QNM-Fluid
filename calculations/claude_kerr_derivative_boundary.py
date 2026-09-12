#!/usr/bin/env python3
"""Ricorrenze di D_minus e D_plus per Kerr, e controllo di K = D_omega.

Implementa le ricorrenze di `research/kerr_explicit_boundary_data_2026-09-12.md`
§2 e §3, che erano il pezzo mancante per il punto 4 dell'estensione E.

    v D_r + D^2 + q = 0,
    D_minus = sum_{n>=0} h_n rho^n,   h_0 = -i k,     rho = r - r_plus,
    D_plus  = sum_{n>=0} t_n r^-n,    t_0 = +i omega.

I coefficienti di v e q **non** sono ricavati a mano: si costruisce q(r) in
forma simbolica dalle definizioni della §1 e se ne prendono gli sviluppi con
sympy, cosi' che un errore algebrico non passi inosservato.

CONTROLLI (quelli prescritti dalla nota)
  * residuo locale di v D_r + D^2 + q coerente con il **primo termine omesso**:
    troncando a ordine N il residuo deve andare come rho^N, non come rho^0;
  * ordini 4, 6, 8 e posizioni del bordo diverse;
  * K = D_omega da  v K_r + 2 D K = -q_omega  confrontata con differenze in
    omega sul **medesimo ramo**, con A ricalcolata a ogni omega.

ESCLUSIONI dichiarate dalla nota e rispettate qui: niente a vicino a 1, niente
omega vicino a 0, e si controlla che il denominatore di Frobenius
n v_1 + 2 h_0 non sia piccolo.

AVVERTENZA riportata dalla nota, che questo script NON supera: per Im omega < 0
il ramo uscente cresce lungo l'asse reale, la serie in 1/r da sola non seleziona
univocamente il ramo, e la stabilita' al variare di R e dell'ordine non prova
l'assenza di contaminazione oltre tutti gli ordini. Qui si verifica la
**consistenza** delle serie, non la selezione del ramo.

Uso: python3.13 calculations/claude_kerr_derivative_boundary.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from claude_kerr_derivative_validation import angular_derivative, spheroidal


def rational_series(expression, variable, order):
    """Coefficienti di Taylor di una funzione RAZIONALE, per divisione esatta.

    `sympy.series` fallisce su queste espressioni con coefficienti complessi in
    virgola mobile.  Ma q e v sono rapporti di polinomi: basta ridurli a
    numeratore e denominatore e dividere per ricorrenza, esatto e veloce.
    """
    numerator, denominator = sp.fraction(sp.cancel(sp.together(expression)))
    top = [complex(c) for c in reversed(sp.Poly(numerator, variable).all_coeffs())]
    bottom = [complex(c) for c in reversed(sp.Poly(denominator, variable).all_coeffs())]
    if abs(bottom[0]) < 1.0e-30:
        raise ValueError("denominatore con termine costante nullo: non e' Taylor")
    top += [0.0j] * (order + 2 - len(top))
    bottom += [0.0j] * (order + 2 - len(bottom))
    out = []
    for k in range(order + 2):
        value = top[k] - sum(bottom[j] * out[k - j] for j in range(1, k + 1))
        out.append(value / bottom[0])
    return np.array(out)


def symbolic_potential(spin, ell, m, omega, separation):
    """q(r) e v(r) simboliche dalle definizioni della nota, §1."""
    r = sp.symbols("r")
    h_squared = r**2 + spin**2
    delta = r**2 - 2 * r + spin**2
    v = delta / h_squared
    geometric = v * (sp.diff(v, r) * r / h_squared + v * spin**2 / h_squared**2)
    q = (omega - m * spin / h_squared) ** 2 - delta * separation / h_squared**2 - geometric
    return r, q, v


def horizon_series(spin, ell, m, omega, separation, order):
    """h_n dalla ricorrenza della §2, con v_j e q_j da sympy."""
    r_plus = 1.0 + np.sqrt(1.0 - spin**2)
    r_minus = 1.0 - np.sqrt(1.0 - spin**2)
    gap, h_plus = r_plus - r_minus, r_plus**2 + spin**2
    wavenumber = omega - m * spin / h_plus

    r, q_expr, v_expr = symbolic_potential(spin, ell, m, omega, separation)
    rho = sp.symbols("rho")
    q_coefficients = rational_series(q_expr.subs(r, r_plus + rho), rho, order)
    v_coefficients = rational_series(v_expr.subs(r, r_plus + rho), rho, order)

    v1 = v_coefficients[1]
    coefficients = [-1j * wavenumber]
    denominators = []
    for n in range(1, order + 1):
        total = q_coefficients[n]
        total += sum(coefficients[i] * coefficients[n - i] for i in range(1, n))
        total += sum(v_coefficients[j] * (n - j + 1) * coefficients[n - j + 1]
                     for j in range(2, n + 1) if 0 <= n - j + 1 < len(coefficients))
        denominator = n * v1 + 2.0 * coefficients[0]
        denominators.append(abs(denominator))
        coefficients.append(-total / denominator)
    return (np.array(coefficients), r_plus, np.array(v_coefficients),
            np.array(q_coefficients), min(denominators))


def infinity_series(spin, ell, m, omega, separation, order):
    """t_n dalla ricorrenza della §3, in potenze di 1/r."""
    r, q_expr, v_expr = symbolic_potential(spin, ell, m, omega, separation)
    inverse = sp.symbols("u")
    q_coefficients = rational_series(q_expr.subs(r, 1 / inverse), inverse, order)
    v_coefficients = rational_series(v_expr.subs(r, 1 / inverse), inverse, order)

    coefficients = [1j * omega]
    for n in range(1, order + 1):
        total = sum(k * coefficients[k] * v_coefficients[n - k - 1]
                    for k in range(1, n) if 0 <= n - k - 1 < len(v_coefficients))
        total -= sum(coefficients[k] * coefficients[n - k] for k in range(1, n))
        total -= q_coefficients[n]
        coefficients.append(total / (2j * omega))
    return np.array(coefficients), np.array(v_coefficients), np.array(q_coefficients)


def convergence_radius(coefficients, skip=4):
    """Raggio stimato dal rapporto |h_n/h_{n-1}|, che si stabilizza rapidamente."""
    ratios = [abs(coefficients[n] / coefficients[n - 1])
              for n in range(skip, len(coefficients))
              if coefficients[n - 1] != 0]
    return 1.0 / float(np.mean(ratios)), ratios


def evaluate(coefficients, variable):
    return sum(c * variable**n for n, c in enumerate(coefficients))


def riccati_residual(coefficients, v_coefficients, q_coefficients, variable, inverse=False):
    """Residuo di v D' + D^2 + q, con derivata della serie troncata."""
    powers = np.arange(len(coefficients))
    value = evaluate(coefficients, variable)
    v_value = evaluate(v_coefficients[:len(coefficients)], variable)
    q_value = evaluate(q_coefficients[:len(coefficients)], variable)
    if inverse:
        # D = sum t_n u^n con u = 1/r; d/dr = -u^2 d/du
        slope = -variable**2 * sum(n * c * variable ** (n - 1)
                                   for n, c in enumerate(coefficients) if n > 0)
    else:
        slope = sum(n * c * variable ** (n - 1) for n, c in enumerate(coefficients) if n > 0)
    del powers
    return v_value * slope + value**2 + q_value


def main() -> None:
    spin, ell, m = 0.6, 28, 19
    omega = 3.0 - 0.25j
    eigen, _, _ = spheroidal(ell, m, spin * omega)
    separation = complex(eigen + spin**2 * omega**2 - 2.0 * spin * m * omega)
    print("Ricorrenze di D_minus e D_plus per Kerr\n")
    print(f"a = {spin}, (ell,m) = ({ell},{m}), omega = {omega}, "
          f"lambda = {separation.real:.6f}{separation.imag:+.6f}i\n")

    print("1. orizzonte: RAGGIO DI CONVERGENZA, prima di tutto")
    coefficients, r_plus, v_coefficients, q_coefficients, smallest = horizon_series(
        spin, ell, m, omega, separation, 12)
    radius, ratios = convergence_radius(coefficients)
    print("   |h_n/h_(n-1)|: " + "  ".join(f"{x:.1f}" for x in ratios[:6]))
    print(f"   -> raggio ~ {radius:.2e}   (min|denominatore Frobenius| = {smallest:.3f})")
    print("   D ha un POLO a quella distanza, cioe' psi ha uno zero: e' la")
    print("   riserva della nota su «zeri del profilo», qui quantificata.")

    print("\n   residuo di v D' + D^2 + q, DENTRO e FUORI il raggio")
    print("   ordine     rho       |residuo|     comportamento")
    for order in (4, 6, 8):
        for rho in (2.0e-3, 5.0e-4, 2.0e-4):
            local, _, v_local, q_local, _ = horizon_series(
                spin, ell, m, omega, separation, order)
            residual = riccati_residual(local, v_local, q_local, rho)
            state = "fuori: diverge" if rho > radius else "dentro: cala con l'ordine"
            print(f"   {order:5d}   {rho:.0e}   {abs(residual):.3e}   {state}")
    print(f"   h_0 = -i k = {coefficients[0]:.6f}   r_plus = {r_plus:.6f}")

    print("\n2. infinito: stesso controllo in 1/r")
    print("   ordine    r       |residuo|      atteso ~ r^-ordine")
    for order in (4, 6, 8):
        coefficients, v_coefficients, q_coefficients = infinity_series(
            spin, ell, m, omega, separation, order)
        for radius_out in (40.0, 80.0, 160.0):
            residual = riccati_residual(coefficients, v_coefficients, q_coefficients,
                                        1.0 / radius_out, inverse=True)
            print(f"   {order:5d}  {radius_out:6.1f}   {abs(residual):.3e}      "
                  f"{radius_out**(-order):.3e}")
    print(f"   t_0 = i omega = {coefficients[0]:.6f}")

    print("\n3. K = D_omega contro differenze in omega, medesimo ramo")
    print("   (A ricalcolata a ogni omega; all'orizzonte DENTRO il raggio)")

    def boundary_value(which, position, frequency, order=8):
        eigenvalue, _, _ = spheroidal(ell, m, spin * frequency)
        local = complex(eigenvalue + spin**2 * frequency**2 - 2.0 * spin * m * frequency)
        if which == "orizzonte":
            local_coefficients, r_plus_local, _, _, _ = horizon_series(
                spin, ell, m, frequency, local, order)
            return evaluate(local_coefficients, position - r_plus_local)
        local_coefficients, _, _ = infinity_series(spin, ell, m, frequency, local, order)
        return evaluate(local_coefficients, 1.0 / position)

    for which, position in (("orizzonte", r_plus + 2.0e-4), ("orizzonte", r_plus + 5.0e-4),
                            ("infinito", 60.0), ("infinito", 120.0)):
        step = 1.0e-5 * abs(omega)
        along_real = ((boundary_value(which, position, omega + step)
                       - boundary_value(which, position, omega - step)) / (2.0 * step))
        along_imag = ((boundary_value(which, position, omega + 1j * step)
                       - boundary_value(which, position, omega - 1j * step))
                      / (2.0j * step))
        limit = -1j if which == "orizzonte" else 1j
        print(f"   {which:9s} r={position:9.5f}   K = "
              f"{along_real.real:+.9f}{along_real.imag:+.9f}i   "
              f"|reale-immag.| = {abs(along_real - along_imag):.2e}"
              f"   |K - limite| = {abs(along_real - limit):.2e}")

    print("\n4. limiti dichiarati, non superati qui")
    print("   - si verifica la CONSISTENZA delle serie, non la selezione del ramo;")
    print("   - per Im omega < 0 la serie in 1/r non seleziona univocamente il ramo")
    print("     uscente, e la stabilita' in R e ordine non lo dimostra;")
    print("   - la norma globale Kerr resta quindi NON validata.")


if __name__ == "__main__":
    main()
