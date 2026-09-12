#!/usr/bin/env python3
"""Gate A del brief del 12 settembre: secondo e terzo metodo spettrale.

Il brief chiede esplicitamente che «cambiare linguaggio/CAS allo stesso shooting
non basta come indipendenza di metodo».  Qui ci sono **tre** metodi con
discretizzazioni diverse:

  1. `shooting`    — integrazione da sinistra, residuo uscente a destra.
                     E' il metodo storico, riportato per confronto.
  2. `matching`    — integrazione dai DUE lati con le rispettive condizioni
                     uscenti, e annullamento del Wronskiano normalizzato in un
                     punto interno.  Propaga gli errori in modo diverso dal (1).
  3. `collocation` — collocazione di Chebyshev, problema agli autovalori
                     **quadratico** in omega (le condizioni al bordo dipendono da
                     omega), linearizzato in forma compagna e risolto in blocco.
                     Nessuna iterazione, nessun valore iniziale: restituisce
                     l'intero spettro discretizzato in un colpo.

Il (3) e' il controllo che conta: non e' uno shooting, non parte da un guess, e
**classifica lo spettro** invece di inseguire una radice.  Il brief chiede di non
chiamare «fondamentale» la risonanza finche' lo spettro non e' classificato.

Convenzioni identiche al resto del compito: V0 = 10 bump(x), perturbazione
eta bump(x; 2.0, 0.4), dominio [-1, 2.4], exp(-i omega t), bordi uscenti
psi'(a) = -i omega psi(a) e psi'(c) = +i omega psi(c).

Uso: python3.13 calculations/claude_compact_spectrum.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eig
from scipy.optimize import root

LEFT, RIGHT = -1.0, 2.4
HEIGHT = 10.0
GUESS = 3.423295488 - 0.607387592j


def bump(x, center=0.0, width=1.0):
    y = (np.asarray(x, dtype=float) - center) / width
    out = np.zeros_like(y, dtype=float)
    inside = np.abs(y) < 1.0
    out[inside] = np.exp(1.0 - 1.0 / (1.0 - y[inside] ** 2))
    return out


def potential(x, eta=0.0):
    return HEIGHT * bump(x) + eta * bump(x, 2.0, 0.4)


# --------------------------------------------------------------------------
# metodo 1: shooting da sinistra
# --------------------------------------------------------------------------
def shooting(eta=0.0, guess=GUESS, tol=1.0e-14):
    def march(omega):
        return solve_ivp(lambda x, y: [y[1], (potential(x, eta) - omega**2) * y[0]],
                         (LEFT, RIGHT), [1.0 + 0j, -1j * omega], method="DOP853",
                         rtol=tol, atol=tol / 100.0, max_step=0.035)

    def residual(v):
        omega = complex(*v)
        out = march(omega)
        defect = out.y[1, -1] - 1j * omega * out.y[0, -1]
        return [defect.real, defect.imag]

    found = root(residual, [guess.real, guess.imag], tol=1.0e-13)
    return complex(*found.x), float(np.linalg.norm(residual(found.x)))


# --------------------------------------------------------------------------
# metodo 2: matching dai due lati, Wronskiano interno
# --------------------------------------------------------------------------
def matching(eta=0.0, guess=GUESS, tol=1.0e-14, meeting=0.7):
    """Wronskiano normalizzato in `meeting` fra i due rami uscenti."""
    def march(omega, start, stop, initial):
        return solve_ivp(lambda x, y: [y[1], (potential(x, eta) - omega**2) * y[0]],
                         (start, stop), initial, method="DOP853",
                         rtol=tol, atol=tol / 100.0, max_step=0.035)

    def residual(v):
        omega = complex(*v)
        left = march(omega, LEFT, meeting, [1.0 + 0j, -1j * omega])
        right = march(omega, RIGHT, meeting, [1.0 + 0j, 1j * omega])
        a, da = left.y[0, -1], left.y[1, -1]
        b, db = right.y[0, -1], right.y[1, -1]
        # Wronskiano normalizzato: insensibile alla scala delle due soluzioni
        wronskian = (a * db - da * b) / (abs(a) * abs(b) + abs(da) * abs(b) / abs(omega))
        return [wronskian.real, wronskian.imag]

    found = root(residual, [guess.real, guess.imag], tol=1.0e-13)
    return complex(*found.x), float(np.linalg.norm(residual(found.x)))


# --------------------------------------------------------------------------
# metodo 3: collocazione di Chebyshev, autovalori quadratici in blocco
# --------------------------------------------------------------------------
def chebyshev(size):
    """Nodi di Chebyshev-Lobatto su [-1,1] e matrice di differenziazione."""
    index = np.arange(size + 1)
    nodes = np.cos(np.pi * index / size)
    weights = np.ones(size + 1)
    weights[0] = weights[-1] = 2.0
    weights *= (-1.0) ** index
    difference = nodes[:, None] - nodes[None, :]
    matrix = np.outer(weights, 1.0 / weights) / (difference + np.eye(size + 1))
    matrix -= np.diag(matrix.sum(axis=1))
    return nodes, matrix


def collocation(eta=0.0, size=180):
    """Spettro completo dal problema quadratico A0 + omega A1 + omega^2 A2."""
    nodes, derivative = chebyshev(size)
    half = 0.5 * (RIGHT - LEFT)
    grid = 0.5 * (RIGHT + LEFT) + half * nodes        # nodi decrescenti: grid[0]=RIGHT
    derivative = derivative / half
    second = derivative @ derivative
    count = size + 1

    a0 = second - np.diag(potential(grid, eta))
    a1 = np.zeros((count, count), dtype=complex)
    a2 = np.eye(count, dtype=complex)

    # bordo destro (grid[0] = RIGHT): psi' - i omega psi = 0
    a0[0, :] = derivative[0, :]
    a1[0, :] = 0.0
    a1[0, 0] = -1j
    a2[0, :] = 0.0
    # bordo sinistro (grid[-1] = LEFT): psi' + i omega psi = 0
    a0[-1, :] = derivative[-1, :]
    a1[-1, :] = 0.0
    a1[-1, -1] = +1j
    a2[-1, :] = 0.0

    # linearizzazione compagna:  [[A0, A1],[0, I]] u = omega [[0, -A2],[I, 0]] u
    zero, identity = np.zeros((count, count), dtype=complex), np.eye(count, dtype=complex)
    left = np.block([[a0, a1], [zero, identity]])
    right = np.block([[zero, -a2], [identity, zero]])
    values = eig(left, right, right=False)
    values = values[np.isfinite(values)]
    return values


def classify(eta=0.0, size=180, window=(0.5, 12.0), damping=3.0):
    """Risonanze nella finestra, ordinate per smorzamento crescente."""
    values = collocation(eta, size)
    keep = [v for v in values
            if window[0] < v.real < window[1] and -damping < v.imag < -1.0e-6]
    return sorted(keep, key=lambda v: -v.imag)


def main() -> None:
    print("Gate A — tre metodi indipendenti per la stessa risonanza\n")

    print("1. Classificazione dello spettro con la collocazione (nessun guess)")
    for size in (140, 180, 220):
        found = classify(size=size)
        print(f"   N={size}:  {len(found)} risonanze con Im>-3, Re in (0.5,12)")
        for value in found[:4]:
            print(f"        {value.real:+.9f}{value.imag:+.9f}i")
    spectrum = classify(size=220)
    least_damped = spectrum[0]
    print(f"\n   la meno smorzata e' {least_damped.real:.9f}{least_damped.imag:+.9f}i")
    print("   (e' quella del brief; 'fondamentale' resta una classificazione,")
    print("    non un'etichetta: dipende dalla finestra scelta)")

    print("\n2. Accordo fra i tre metodi, a eta = 0, +/-0.001, +/-0.01")
    print("   eta      shooting                    matching                    "
          "collocazione            |sh-ma|    |sh-co|")
    drifts = {}
    for eta in (0.0, 1.0e-3, -1.0e-3, 1.0e-2, -1.0e-2):
        one, _ = shooting(eta)
        two, _ = matching(eta)
        three = classify(eta, size=220)[0]
        drifts[eta] = one
        print(f"   {eta:+.3f}  {one.real:.9f}{one.imag:+.9f}i  "
              f"{two.real:.9f}{two.imag:+.9f}i  {three.real:.6f}{three.imag:+.6f}i  "
              f"{abs(one-two):.2e}  {abs(one-three):.2e}")

    print("\n3. Le variazioni in eta sono risolte molto sopra l'errore fra metodi?")
    base = drifts[0.0]
    for eta in (1.0e-3, 1.0e-2):
        change = abs(drifts[eta] - base)
        one, _ = shooting(eta)
        two, _ = matching(eta)
        gap = abs(one - two)
        print(f"   eta={eta:.0e}:  |Delta omega| = {change:.3e}   errore fra metodi "
              f"= {gap:.1e}   rapporto = {change/max(gap,1e-16):.1e}")


if __name__ == "__main__":
    main()
