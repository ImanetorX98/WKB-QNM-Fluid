#!/usr/bin/env python3
"""Estensione E del brief: derivata spettrale di Kerr, punti 1-4.

FORMULE DA VERIFICARE (Codex, `kerr_nonlinear_spectral_proof_2026-09-12.md`)

    A_c = -2c [int_0^pi S^2 cos^2(theta) sin(theta) dtheta]
              / [int_0^pi S^2 sin(theta) dtheta],            pairing BILINEARE

    q_omega = 2(omega - m a/H) - (Delta/H^2)[a A_c + 2 a^2 omega - 2 a m].

Nella base armonica sferica, ortonormale rispetto a int ... sin(theta) dtheta,
con S = sum_j v_j Y_j si ha

    int S^2 sin dtheta = v^T v          (trasposto, NON coniugato)
    int S^2 cos^2 sin dtheta = v^T M2 v,   M2 = (matrice cos theta)^2

quindi A_c = -2c (v^T M2 v)/(v^T v).  Il rapporto e' invariante per riscalatura,
quindi non dipende dalla normalizzazione dell'autovettore; ma **v^T v e' la
quantita' da registrare**, perche' e' lei ad annullarsi vicino a
auto-ortogonalita' e punti eccezionali, dove la formula non vale.

`numpy.linalg.eig` normalizza a v^dagger v = 1, non v^T v = 1: per c complessa i
due differiscono, ed e' esattamente l'errore che il brief vieta
(«NON sostituire S^2 con |S|^2»).

PUNTI COPERTI
  1. A_c bilineare contro differenze centrali lungo le direzioni reale e
     immaginaria del piano c, per (ell,m) = (28,19), (61,41) — entrambi con
     m/L = 2/3 esatto — e c/L = 0.2, 0.4, 0.3+0.1i.
  2. Tre passi relativi 1e-3, 1e-4, 1e-5 di max(1,|c|), e due dimensioni di
     base.  Denominatore bilineare normalizzato riportato.
  3. q_omega a tre raggi esterni all'orizzonte, per a = 0 e 0.6, differenziando
     l'INTERO potenziale e ricalcolando A a ogni omega; piu' la variante con A
     congelata, etichettata come controllo negativo.
  4. Limite a = 0: q_omega deve dare 2 omega.

Uso: python3.13 calculations/claude_kerr_derivative_validation.py
"""

from __future__ import annotations

import numpy as np

from kerr_eikonal_order_test import cosine_matrix


def spheroidal(ell, m, c, pad=60, steps=24):
    """A e autovettore, stesso tracciamento in due fasi del modulo storico."""
    ells, cosine = cosine_matrix(m, ell + pad)
    squared = cosine @ cosine
    diagonal = np.diag(ells * (ells + 1.0))
    index = ell - abs(m)

    real = float(np.real(c))
    values, vectors = np.linalg.eigh(diagonal - (real**2) * squared)
    value = complex(values[index])
    vector = vectors[:, index].astype(complex)

    imaginary = float(np.imag(c))
    if imaginary != 0.0:
        complex_diagonal = diagonal.astype(complex)
        for step in range(1, steps + 1):
            partial = real + 1j * imaginary * step / steps
            values, vectors = np.linalg.eig(complex_diagonal - (partial**2) * squared)
            overlap = np.abs(vector.conj() @ vectors)
            best = int(np.argmax(overlap))
            value, vector = complex(values[best]), vectors[:, best]
    return value, vector, squared


def angular_derivative(ell, m, c, pad=60):
    """A_c dalla formula bilineare, piu' il denominatore normalizzato."""
    _, vector, squared = spheroidal(ell, m, c, pad)
    vector = vector / np.sqrt(np.vdot(vector, vector))      # v^dagger v = 1
    bilinear = complex(vector @ vector)                      # v^T v, NON |v|^2
    weighted = complex(vector @ (squared @ vector))
    return -2.0 * c * weighted / bilinear, bilinear


def potential_pieces(radius, spin, ell, m, omega, pad=60):
    """lambda e A alla frequenza data: A viene RICALCOLATA a ogni omega."""
    eigen, _, _ = spheroidal(ell, m, spin * omega, pad)
    return eigen + spin**2 * omega**2 - 2.0 * spin * m * omega, eigen


def full_potential(radius, spin, ell, m, omega, pad=60, freeze=None):
    """q(r, omega) completo.  `freeze` congela lambda: controllo negativo."""
    h_squared = radius**2 + spin**2
    delta = radius**2 - 2.0 * radius + spin**2
    if freeze is None:
        separation, _ = potential_pieces(radius, spin, ell, m, omega, pad)
    else:
        separation = freeze
    lapse = delta / h_squared
    slope = (2.0 * radius * h_squared - 2.0 * radius * delta) / h_squared**2
    geometric = lapse * (slope * radius / h_squared + lapse * spin**2 / h_squared**2)
    return (omega - m * spin / h_squared) ** 2 - delta * separation / h_squared**2 - geometric


def main() -> None:
    print("Estensione E — derivata spettrale di Kerr\n")

    print("1-2. A_c bilineare contro differenze centrali")
    print("   ell   m     c            A_c bilineare                passo    "
          "reale       immag.     v^T v")
    for ell, m in ((28, 19), (61, 41)):
        for c in (0.2 * (ell + 0.5), 0.4 * (ell + 0.5), (0.3 + 0.1j) * (ell + 0.5)):
            analytic, bilinear = angular_derivative(ell, m, c)
            line = f"   {ell:3d}  {m:3d}  {c.real:6.2f}{c.imag:+6.2f}i  " \
                   f"{analytic.real:+.9f}{analytic.imag:+.9f}i"
            for step_ratio in (1.0e-3, 1.0e-4, 1.0e-5):
                step = step_ratio * max(1.0, abs(c))
                along_real = (spheroidal(ell, m, c + step)[0]
                              - spheroidal(ell, m, c - step)[0]) / (2.0 * step)
                along_imag = (spheroidal(ell, m, c + 1j * step)[0]
                              - spheroidal(ell, m, c - 1j * step)[0]) / (2.0j * step)
                print(f"{line}  {step_ratio:.0e}  {abs(along_real/analytic-1):.2e}"
                      f"  {abs(along_imag/analytic-1):.2e}  {abs(bilinear):.4f}")
                line = " " * len(line)

    print("\n   dimensione della base: pad 40 contro 80")
    for ell, m in ((28, 19), (61, 41)):
        c = 0.4 * (ell + 0.5)
        small, _ = angular_derivative(ell, m, c, pad=40)
        large, _ = angular_derivative(ell, m, c, pad=80)
        print(f"   ell={ell}: |A_c(40)/A_c(80) - 1| = {abs(small/large-1):.2e}")

    print("\n3. q_omega a tre raggi, differenziando l'intero potenziale")
    print("   a     r      q_omega formula            differenze centrali       "
          "scarto      A congelata (neg.)")
    for spin, ell, m in ((0.0, 28, 19), (0.6, 28, 19)):
        omega = 3.0 - 0.25j
        horizon = 1.0 + np.sqrt(max(1.0 - spin**2, 0.0))
        for radius in (1.6 * horizon, 3.0, 6.0):
            h_squared = radius**2 + spin**2
            delta = radius**2 - 2.0 * radius + spin**2
            derivative, _ = angular_derivative(ell, m, spin * omega)
            formula = (2.0 * (omega - m * spin / h_squared)
                       - delta / h_squared**2
                       * (spin * derivative + 2.0 * spin**2 * omega - 2.0 * spin * m))
            step = 1.0e-5 * max(1.0, abs(omega))
            numeric = ((full_potential(radius, spin, ell, m, omega + step)
                        - full_potential(radius, spin, ell, m, omega - step))
                       / (2.0 * step))
            frozen_value, _ = potential_pieces(radius, spin, ell, m, omega)
            frozen = ((full_potential(radius, spin, ell, m, omega + step, freeze=frozen_value)
                       - full_potential(radius, spin, ell, m, omega - step, freeze=frozen_value))
                      / (2.0 * step))
            print(f"   {spin:.1f}  {radius:5.2f}  {formula.real:+.9f}{formula.imag:+.9f}i  "
                  f"{numeric.real:+.9f}{numeric.imag:+.9f}i  {abs(numeric/formula-1):.2e}"
                  f"   {abs(frozen/formula-1):.2e}")

    print("\n4. limite a = 0: q_omega deve dare 2 omega")
    omega = 3.0 - 0.25j
    for radius in (3.0, 6.0, 12.0):
        h_squared, delta = radius**2, radius**2 - 2.0 * radius
        derivative, _ = angular_derivative(28, 19, 0.0)
        formula = 2.0 * omega - delta / h_squared**2 * (0.0 * derivative + 0.0 - 0.0)
        print(f"   r={radius:5.1f}   q_omega = {formula.real:+.12f}{formula.imag:+.12f}i"
              f"   2 omega = {(2*omega).real:+.12f}{(2*omega).imag:+.12f}i"
              f"   scarto {abs(formula - 2*omega):.2e}")
    print("\n   Nota: a=0 annulla il coefficiente di A_c, quindi il limite e'")
    print("   esatto per costruzione.  Il controllo utile e' il punto 3, dove")
    print("   A_c entra con peso non nullo.")


if __name__ == "__main__":
    main()
