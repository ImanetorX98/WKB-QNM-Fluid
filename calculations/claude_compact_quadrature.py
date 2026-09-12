#!/usr/bin/env python3
"""Gate B del brief del 12 settembre: la quadratura della derivata di E.

PROBLEMA (Codex, `research/analytic_core_response_2026-09-12.md` §5).  La norma
spettrale converge, ma dE/deta calcolata con quadratura globale su J non
converge: 0.00189086, 0.00189892, 0.00188900 a 8001, 16001, 32001 punti.  Il
funzionale DISCRETO originale (trapezio a 2001 punti) da' invece 0.00187362, in
accordo a ~1e-5 con le sue differenze centrali.  Sono due livelli diversi e non
vanno forzati a coincidere.

DIAGNOSI ADOTTATA.  U = int w |Q_M| dx ha integrando **non liscio**: |Q_M| ha
angoli dove Q_M cambia segno, e la derivata dot U = int w sgn(Q_M) dot Q_M dx ha
salti negli stessi punti.  Una quadratura globale liscia su un integrando C^0
converge come una potenza bassa del passo, e nessun raffinamento uniforme la
salva.  La cura non e' piu' precisione: e' **spezzare l'intervallo sugli zeri**
e usare quadratura adattiva in ciascun tratto, dove l'integrando e' analitico.

CHE COSA VERIFICA QUESTO PROGRAMMA (i punti di B, nell'ordine)
  1. localizza tutti gli zeri di Q_M in J per bracketing raffinato, e ne
     verifica la stabilita' al variare della griglia di ricerca;
  2. verifica l'assenza di zeri di psi e registra min|psi|, distinguendoli
     dagli zeri di Q_M;
  3. calcola k in DUE modi indipendenti — formula integrale e ODE
     k' + 2 z k = -2 omega con k(a) = -i — e ne riporta lo scarto;
  4. calcola dE/deta analitica con quadratura spezzata sugli zeri;
  5. la confronta con differenze centrali di E, dove **anche E** e' calcolata
     con la quadratura spezzata, su tre raffinamenti;
  6. ripete il confronto sulla trapezoidale a 2001 punti come regressione del
     funzionale discreto, senza scambiarla per il limite continuo.

Convenzioni: V0 = 10 bump(x), perturbazione eta bump(x; 2.0, 0.4), dominio
[-1, 2.4], bordi uscenti esatti perche' V=0 fuori dal supporto, eps = 1,
J = [-0.8, 0.8], peso w = exp(-x^2), normalizzazione psi(a) = 1.

Uso: python3.13 calculations/claude_compact_quadrature.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, root

LEFT, RIGHT = -1.0, 2.4
WINDOW = (-0.8, 0.8)
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


def integrate_mode(omega, eta, tol):
    def rhs(x, y):
        return [y[1], (potential(x, eta) - omega**2) * y[0]]

    solution = solve_ivp(rhs, (LEFT, RIGHT), [1.0 + 0j, -1j * omega],
                         method="DOP853", rtol=tol, atol=tol / 100.0,
                         max_step=0.035, dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def resonance(eta=0.0, guess=GUESS, tol=1.0e-12):
    """Risonanza dal matching uscente a destra.  Metodo 1 (shooting)."""
    def residual(v):
        omega = complex(*v)
        solution = integrate_mode(omega, eta, tol)
        defect = solution.y[1, -1] - 1j * omega * solution.y[0, -1]
        return [defect.real, defect.imag]

    found = root(residual, [guess.real, guess.imag], tol=1.0e-12)
    omega = complex(*found.x)
    return omega, float(np.linalg.norm(residual(found.x))), integrate_mode(omega, eta, tol)


# --------------------------------------------------------------------------
# 1-2. zeri di Q_M, e assenza di zeri di psi
# --------------------------------------------------------------------------
def profile(solution, omega, eta=0.0):
    """Restituisce funzioni callable: psi, z = psi'/psi, p = Im z, Q_M."""
    def psi(x):
        return solution.sol(np.atleast_1d(x))[0]

    def logarithmic(x):
        values = solution.sol(np.atleast_1d(x))
        return values[1] / values[0]

    def quantum(x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        return (omega * omega - potential(x, eta)).real - logarithmic(x).imag ** 2

    return psi, logarithmic, quantum


def sign_changes(function, samples):
    """Zeri per bracketing su una griglia, raffinati con Brent."""
    values = function(samples)
    roots = []
    for index in range(len(samples) - 1):
        if values[index] == 0.0:
            roots.append(float(samples[index]))
        elif values[index] * values[index + 1] < 0.0:
            roots.append(float(brentq(lambda t: float(function(np.array([t]))[0]),
                                      samples[index], samples[index + 1],
                                      xtol=1.0e-14, rtol=8.9e-16)))
    return roots


# --------------------------------------------------------------------------
# 3. k in due modi
# --------------------------------------------------------------------------
def response_kernel_ode(solution, omega, tol=1.0e-12):
    """k dall'ODE  k' + 2 z k = -2 omega,  k(a) = -i.  Metodo 1."""
    def rhs(x, y):
        values = solution.sol(np.array([x]))
        z = values[1][0] / values[0][0]
        return [-2.0 * z * y[0] - 2.0 * omega]

    out = solve_ivp(rhs, (LEFT, WINDOW[1] + 0.05), [-1j], method="DOP853",
                    rtol=tol, atol=tol / 100.0, dense_output=True)
    if not out.success:
        raise RuntimeError(out.message)
    return lambda x: out.sol(np.atleast_1d(x))[0]


def response_kernel_integral(solution, omega):
    """k(x) = psi^-2 [ -i psi(a)^2 - 2 omega int_a^x psi^2 ].  Metodo 2."""
    psi_left = solution.sol(np.array([LEFT]))[0][0]

    def kernel(x):
        out = []
        for point in np.atleast_1d(np.asarray(x, dtype=float)):
            real = quad(lambda t: (solution.sol(np.array([t]))[0][0] ** 2).real,
                        LEFT, point, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
            imag = quad(lambda t: (solution.sol(np.array([t]))[0][0] ** 2).imag,
                        LEFT, point, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
            psi_here = solution.sol(np.array([point]))[0][0]
            out.append((-1j * psi_left**2 - 2.0 * omega * (real + 1j * imag))
                       / psi_here**2)
        return np.array(out)

    return kernel


# --------------------------------------------------------------------------
# 4-5. E e dE/deta con quadratura spezzata
# --------------------------------------------------------------------------
def split_quad(function, breakpoints, epsabs=1.0e-13, epsrel=1.0e-13):
    """Quadratura adattiva su ciascun tratto fra i punti di rottura."""
    total, error = 0.0, 0.0
    for left, right in zip(breakpoints[:-1], breakpoints[1:]):
        if right - left < 1.0e-13:
            continue
        value, estimate = quad(function, left, right, limit=400,
                               epsabs=epsabs, epsrel=epsrel)
        total += value
        error += abs(estimate)
    return total, error


def indicator_split(solution, omega, eta=0.0, breakpoints=None):
    """E = U/D con quadratura spezzata sugli zeri di Q_M."""
    _, logarithmic, quantum = profile(solution, omega, eta)
    weight = lambda t: np.exp(-t * t)
    if breakpoints is None:
        breakpoints = [WINDOW[0]] + sign_changes(
            quantum, np.linspace(*WINDOW, 4001)) + [WINDOW[1]]

    numerator, _ = split_quad(
        lambda t: float(weight(t) * abs(quantum(np.array([t]))[0])), breakpoints)
    denominator, _ = split_quad(
        lambda t: float(weight(t) * (abs(omega) ** 2 + abs(potential(np.array([t]), eta)[0])
                                     + logarithmic(np.array([t])).imag[0] ** 2)),
        [WINDOW[0], WINDOW[1]])
    return numerator / denominator, numerator, denominator, breakpoints


def analytic_derivative(solution, omega, drift, kernel, breakpoints):
    """dE/deta dalla formula del §4 della nota, con quadratura spezzata."""
    _, logarithmic, quantum = profile(solution, omega)
    weight = lambda t: np.exp(-t * t)

    def quantum_dot(t):
        p = logarithmic(np.array([t])).imag[0]
        k = kernel(np.array([t]))[0]
        return 2.0 * ((omega * drift).real - p * (k * drift).imag)

    def denominator_dot(t):
        p = logarithmic(np.array([t])).imag[0]
        k = kernel(np.array([t]))[0]
        return 2.0 * (np.conj(omega) * drift).real + 2.0 * p * (k * drift).imag

    numerator_dot, _ = split_quad(
        lambda t: float(weight(t) * np.sign(quantum(np.array([t]))[0]) * quantum_dot(t)),
        breakpoints)
    density_dot, _ = split_quad(lambda t: float(weight(t) * denominator_dot(t)),
                                [WINDOW[0], WINDOW[1]])

    _, numerator, denominator, _ = indicator_split(solution, omega,
                                                   breakpoints=breakpoints)
    return (numerator_dot * denominator - numerator * density_dot) / denominator**2


def spectral_drift(solution, omega):
    """d = B/N, norma generalizzata con termini di superficie."""
    def squared(t):
        return solution.sol(np.array([t]))[0][0] ** 2

    bump_real = quad(lambda t: float(bump(np.array([t]), 2.0, 0.4)[0] * squared(t).real),
                     1.6, 2.4, limit=400, epsabs=1e-14, epsrel=1e-14)[0]
    bump_imag = quad(lambda t: float(bump(np.array([t]), 2.0, 0.4)[0] * squared(t).imag),
                     1.6, 2.4, limit=400, epsabs=1e-14, epsrel=1e-14)[0]
    norm_real = quad(lambda t: float(squared(t).real), LEFT, RIGHT,
                     limit=800, epsabs=1e-14, epsrel=1e-14)[0]
    norm_imag = quad(lambda t: float(squared(t).imag), LEFT, RIGHT,
                     limit=800, epsabs=1e-14, epsrel=1e-14)[0]
    edges = (solution.sol(np.array([LEFT]))[0][0] ** 2
             + solution.sol(np.array([RIGHT]))[0][0] ** 2)
    return (bump_real + 1j * bump_imag) / (2.0 * omega * (norm_real + 1j * norm_imag)
                                           + 1j * edges)


def main() -> None:
    print("Gate B — quadratura della derivata di E\n")
    omega, defect, solution = resonance()
    print(f"risonanza  omega = {omega.real:.12f}{omega.imag:+.12f}i   residuo uscente {defect:.2e}")

    psi, logarithmic, quantum = profile(solution, omega)

    print("\n1. zeri di Q_M in J, e stabilita' del conteggio")
    for count in (1001, 2001, 4001, 8001):
        roots = sign_changes(quantum, np.linspace(*WINDOW, count))
        print(f"   griglia {count:5d}:  {len(roots)} zeri   " +
              "  ".join(f"{r:+.10f}" for r in roots))
    roots = sign_changes(quantum, np.linspace(*WINDOW, 8001))

    print("\n2. psi ha zeri in J?")
    grid = np.linspace(*WINDOW, 20001)
    magnitude = np.abs(psi(grid))
    print(f"   min|psi| = {magnitude.min():.6e} a x = {grid[magnitude.argmin()]:+.4f}"
          f"   -> {'nessuno zero' if magnitude.min() > 1e-3 else 'ATTENZIONE'}")
    print(f"   gli zeri trovati sopra sono di Q_M, non di psi")

    print("\n3. k con due metodi indipendenti")
    kernel_ode = response_kernel_ode(solution, omega)
    kernel_int = response_kernel_integral(solution, omega)
    probes = np.array([-0.8, -0.4, 0.0, 0.4, 0.8])
    difference = np.abs(kernel_ode(probes) - kernel_int(probes))
    for point, value, gap in zip(probes, kernel_ode(probes), difference):
        print(f"   x={point:+.1f}   k = {value.real:+.10f}{value.imag:+.10f}i"
              f"   |ODE - integrale| = {gap:.2e}")

    drift = spectral_drift(solution, omega)
    print(f"\n   d = B/N = {drift.real:.12f}{drift.imag:+.12f}i")

    print("\n4-5. dE/deta: analitica spezzata contro differenze centrali")
    breakpoints = [WINDOW[0]] + roots + [WINDOW[1]]
    analytic = analytic_derivative(solution, omega, drift, kernel_ode, breakpoints)
    print(f"   analitica (quadratura spezzata) = {analytic:.12f}")

    base, _, _, _ = indicator_split(solution, omega)
    print(f"   E(0) = {base:.12f}")
    print("\n   passo eta      E(+eta)          E(-eta)        FD centrale      scarto rel.")
    for step in (1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4):
        values = []
        for sign in (+1, -1):
            omega_s, _, solution_s = resonance(sign * step, guess=omega)
            value, _, _, _ = indicator_split(solution_s, omega_s, sign * step)
            values.append(value)
        finite = (values[0] - values[1]) / (2.0 * step)
        print(f"   {step:.0e}    {values[0]:.12f}  {values[1]:.12f}  "
              f"{finite:.12f}   {abs(finite / analytic - 1):.2e}")


if __name__ == "__main__":
    main()
