#!/usr/bin/env python3
"""Gate C del brief del 12 settembre: fattorizzazione e rango della risposta.

TESI DA FALSIFICARE (Codex, `analytic_core_response_2026-09-12.md` §3).  Per
perturbazioni **esterne** — nulle fra il bordo sinistro e la finestra J — la
risposta linearizzata di qualunque diagnostico interno passa interamente per la
frequenza.  Esplicitamente, con d = d_R + i d_I,

    partial_eta Q_M = d_R f1(x) + d_I f2(x),
    f1 = 2[Re omega - p Im k],      f2 = -2[Im omega + p Re k],

dove p = Im z e k risolve k' + 2 z k = -2 omega, k(a) = -i.  Le funzioni f1, f2
**non dipendono dalla perturbazione**: solo i due coefficienti reali lo fanno.

Ne segue che la matrice delle risposte, una colonna per perturbazione, ha rango
reale **al piu' due**, per quante perturbazioni esterne si vogliano.  E' una
predizione forte e falsificabile: basta una terza direzione singolare sopra il
rumore per smentirla.

PROTOCOLLO (i punti del brief)
  * famiglia b_{c,w} con c = 1.6, 2.0, 2.4 e w = 0.2, 0.4 — tutte esterne a
    J = [-0.8, 0.8] e alla barriera centrale [-1, 1];
  * dominio esteso a [-1, 3.2] per contenere il supporto fino a 2.8, con il
    bordo libero **fuori** da ogni bump;
  * per ciascuna, d e h = partial_eta z per differenze centrali, e confronto
    h(x) = d k(x) con metrica L2 relativa piu' errore assoluto;
  * SVD delle colonne REALI partial_eta Q_M campionate sullo stesso J, con la
    terza singolare confrontata con il rumore **stimato**, non con una soglia
    arbitraria;
  * predizione di ciascuna colonna con la base {f1, f2} **senza rifittare d**;
  * controllo negativo con perturbazione INTERNA (centro 0.3, larghezza 0.1),
    dove l'ipotesi non vale e h = d k non deve reggere.

Uso: python3.13 calculations/claude_compact_rank.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

LEFT, RIGHT = -1.0, 3.2
WINDOW = (-0.8, 0.8)
HEIGHT = 10.0
GUESS = 3.423295488 - 0.607387592j
STEP = 1.0e-4


def bump(x, center=0.0, width=1.0):
    y = (np.asarray(x, dtype=float) - center) / width
    out = np.zeros_like(y, dtype=float)
    inside = np.abs(y) < 1.0
    out[inside] = np.exp(1.0 - 1.0 / (1.0 - y[inside] ** 2))
    return out


def make_potential(center, width, eta):
    def potential(x):
        return HEIGHT * bump(x) + eta * bump(x, center, width)
    return potential


def solve(center, width, eta, guess=GUESS, tol=1.0e-14):
    potential = make_potential(center, width, eta)

    def march(omega, dense=False):
        return solve_ivp(lambda x, y: [y[1], (potential(x) - omega**2) * y[0]],
                         (LEFT, RIGHT), [1.0 + 0j, -1j * omega], method="DOP853",
                         rtol=tol, atol=tol / 100.0, max_step=0.035,
                         dense_output=dense)

    def residual(v):
        omega = complex(*v)
        out = march(omega)
        defect = out.y[1, -1] - 1j * omega * out.y[0, -1]
        return [defect.real, defect.imag]

    found = root(residual, [guess.real, guess.imag], tol=1.0e-13)
    omega = complex(*found.x)
    return omega, march(omega, True)


def logarithmic(solution, grid):
    values = solution.sol(grid)
    return values[1] / values[0]


def kernel(solution, omega, grid, tol=1.0e-13):
    """k dall'ODE k' + 2 z k = -2 omega, k(a) = -i."""
    def rhs(x, y):
        values = solution.sol(np.array([x]))
        return [-2.0 * (values[1][0] / values[0][0]) * y[0] - 2.0 * omega]

    out = solve_ivp(rhs, (LEFT, WINDOW[1] + 0.05), [-1j], method="DOP853",
                    rtol=tol, atol=tol / 100.0, dense_output=True)
    if not out.success:
        raise RuntimeError(out.message)
    return out.sol(grid)[0]


def response(center, width, grid, step=STEP):
    """d, h = partial_eta z, e partial_eta Q_M, per differenze centrali."""
    omega0, base = solve(center, width, 0.0)
    plus_omega, plus = solve(center, width, +step, guess=omega0)
    minus_omega, minus = solve(center, width, -step, guess=omega0)

    drift = (plus_omega - minus_omega) / (2.0 * step)
    slope = (logarithmic(plus, grid) - logarithmic(minus, grid)) / (2.0 * step)

    def quantum(solution, omega, eta):
        potential = make_potential(center, width, eta)
        return (omega * omega - potential(grid)).real - logarithmic(solution, grid).imag ** 2

    quantum_dot = (quantum(plus, plus_omega, +step)
                   - quantum(minus, minus_omega, -step)) / (2.0 * step)
    return {"omega": omega0, "solution": base, "d": drift,
            "h": slope, "dQ": quantum_dot}


def basis(solution, omega, grid):
    """f1, f2: le due funzioni reali della §3, indipendenti dalla perturbazione."""
    p = logarithmic(solution, grid).imag
    k = kernel(solution, omega, grid)
    return 2.0 * (omega.real - p * k.imag), -2.0 * (omega.imag + p * k.real)


def relative(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def main() -> None:
    grid = np.linspace(*WINDOW, 801)
    print("Gate C — fattorizzazione e rango della risposta\n")

    omega0, base_solution = solve(2.0, 0.4, 0.0)
    print(f"dominio [{LEFT}, {RIGHT}],  omega(0) = {omega0.real:.12f}{omega0.imag:+.12f}i")
    peak = HEIGHT * bump(np.linspace(-1, 1, 20001)).max()
    print(f"massimo dominante della barriera centrale = {peak:.10f} (invariato: le "
          f"perturbazioni sono fuori da [-1,1])")

    family = [(c, w) for c in (1.6, 2.0, 2.4) for w in (0.2, 0.4)]
    print("\n1. h(x) = d k(x) per ciascuna perturbazione esterna")
    print("   centro  largh.   d                                  L2 rel.    err. ass. max")
    columns, drifts = [], []
    for center, width in family:
        out = response(center, width, grid)
        predicted = out["d"] * kernel(out["solution"], out["omega"], grid)
        print(f"   {center:5.1f}  {width:5.2f}   {out['d'].real:+.9f}{out['d'].imag:+.9f}i"
              f"   {relative(out['h'], predicted):.2e}   "
              f"{np.abs(out['h'] - predicted).max():.2e}")
        columns.append(out["dQ"])
        drifts.append(out["d"])

    print("\n2. rango della matrice delle risposte reali (6 colonne su J)")
    matrix = np.array(columns).T
    singular = np.linalg.svd(matrix, compute_uv=False)
    print("   valori singolari: " + "  ".join(f"{s:.3e}" for s in singular))
    print(f"   rapporto sigma3/sigma1 = {singular[2] / singular[0]:.2e}")

    # rumore stimato: rifare una colonna con passo eta diverso e misurarne lo scarto
    reference = response(2.0, 0.4, grid, step=STEP)["dQ"]
    coarse = response(2.0, 0.4, grid, step=3.0 * STEP)["dQ"]
    noise = float(np.linalg.norm(reference - coarse))
    print(f"   rumore stimato (stessa colonna, passo eta x3) = {noise:.3e}")
    print(f"   sigma3 = {singular[2]:.3e}   ->  "
          f"{'compatibile col rumore' if singular[2] < 5 * noise else 'SOPRA il rumore'}")

    print("\n3. predizione con la base {f1, f2}, senza rifittare d")
    first, second = basis(base_solution, omega0, grid)
    print("   centro  largh.   L2 rel. della predizione d_R f1 + d_I f2")
    for (center, width), column, drift in zip(family, columns, drifts):
        predicted = drift.real * first + drift.imag * second
        print(f"   {center:5.1f}  {width:5.2f}   {relative(column, predicted):.2e}")

    print("\n4. controllo negativo: perturbazione INTERNA (centro 0.3, largh. 0.1)")
    inner = response(0.3, 0.1, grid)
    predicted_inner = inner["d"] * kernel(inner["solution"], inner["omega"], grid)
    print(f"   d = {inner['d'].real:+.9f}{inner['d'].imag:+.9f}i")
    print(f"   L2 rel. di h contro d k = {relative(inner['h'], predicted_inner):.2e}"
          f"   <- deve essere GRANDE: l'ipotesi non vale")
    predicted_column = inner["d"].real * first + inner["d"].imag * second
    print(f"   L2 rel. della colonna contro la base = "
          f"{relative(inner['dQ'], predicted_column):.2e}")


if __name__ == "__main__":
    main()
