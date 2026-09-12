#!/usr/bin/env python3
"""Gate D del brief del 12 settembre: utilita' predittiva, separata dalla struttura.

LA DOMANDA.  Il gate C mostra che il diagnostico interno **risponde** alle
perturbazioni esterne, con rango due.  Rispondere non e' prevedere.  Qui si
chiede: E predice l'errore della WKB3?

IL TEST E' NETTO PER COSTRUZIONE.  Le perturbazioni della famiglia sono nulle su
[-1, 1], quindi il getto centrale di V e' **invariato a tutti gli ordini**:
Lambda_2, Lambda_3 e la frequenza WKB3 non si muovono di una cifra.  Ma la
frequenza vera si muove.  Dunque l'errore della WKB3 cambia **solo** perche'
cambia il bersaglio, e si puo' chiedere se E lo segua.

Se E rispondesse in modo proporzionale all'errore, sarebbe un predittore; se
risponde senza seguirlo, e' una diagnosi a posteriori della soluzione, non una
previsione.  Il brief impone anche di dichiarare training e test **prima** di
qualunque fit: qui non si fitta nulla, si confrontano andamenti.

DEFINIZIONI (il brief chiede di dichiararle esattamente)
  * peso w(x) = exp(-x^2) su J = [-0.8, 0.8];
  * rapporto integrale  E_int = int_J w |Q_M| dx / int_J w [|omega|^2+|V0|+p^2] dx,
    con quadratura spezzata sugli zeri di Q_M (gate B);
  * rapporto delle mediane pesate E_med = med_w(|Q_M|) / med_w(|omega|^2+|V0|+p^2),
    dove med_w(g) e' il valore g(x*) nel punto x* in cui la somma cumulata dei
    pesi ordinata per g crescente raggiunge meta' del peso totale.  Quantile
    discreto, nessuna interpolazione, griglia uniforme di 4001 punti su J.
  * le due NON soddisfano la stessa formula di derivata: la mediana non e'
    differenziabile in presenza di plateaux (nota di Codex §4), e qui e' usata
    solo come valore, mai derivata.

Uso: python3.13 calculations/claude_compact_predictive.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, root

LEFT, RIGHT = -1.0, 3.2
WINDOW = (-0.8, 0.8)
HEIGHT = 10.0
GUESS = 3.423295488 - 0.607387592j


def bump(x, center=0.0, width=1.0):
    y = (np.asarray(x, dtype=float) - center) / width
    out = np.zeros_like(y, dtype=float)
    inside = np.abs(y) < 1.0
    out[inside] = np.exp(1.0 - 1.0 / (1.0 - y[inside] ** 2))
    return out


def make_potential(center, width, eta):
    return lambda x: HEIGHT * bump(x) + eta * bump(x, center, width)


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


# --------------------------------------------------------------------------
# WKB3 dal getto centrale
# --------------------------------------------------------------------------
def central_jet():
    """V, V'', ..., V^(6) del bump centrale in x=0, in forma esatta."""
    x = sp.symbols("x")
    expression = HEIGHT * sp.exp(1 - 1 / (1 - x**2))
    return [float(sp.diff(expression, x, k).subs(x, 0)) for k in range(7)]


def wkb3(overtone=0):
    values = central_jet()
    v0, _, v2, v3, v4, v5, v6 = values
    alpha = overtone + 0.5
    root_curvature = np.sqrt(-2.0 * v2)
    lambda2 = (0.125 * (v4 / v2) * (0.25 + alpha**2)
               - (v3 / v2) ** 2 * (7.0 + 60.0 * alpha**2) / 288.0) / root_curvature
    lambda3 = 1.0 / (-2.0 * v2) * (
        5.0 * (v3 / v2) ** 4 * (77.0 + 188.0 * alpha**2) / 6912.0
        - (v3**2 * v4 / v2**3) * (51.0 + 100.0 * alpha**2) / 384.0
        + (v4 / v2) ** 2 * (67.0 + 68.0 * alpha**2) / 2304.0
        + (v3 * v5 / v2**2) * (19.0 + 28.0 * alpha**2) / 288.0
        - (v6 / v2) * (5.0 + 4.0 * alpha**2) / 288.0)
    squared = v0 + root_curvature * lambda2 - 1j * alpha * root_curvature * (1 + lambda3)
    omega = np.sqrt(squared)
    return (omega if omega.real > 0 else -omega), lambda2, lambda3


# --------------------------------------------------------------------------
# i due diagnostici
# --------------------------------------------------------------------------
def pieces(solution, omega, potential):
    def logarithmic(t):
        values = solution.sol(np.atleast_1d(t))
        return values[1] / values[0]

    def quantum(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        return (omega * omega - potential(t)).real - logarithmic(t).imag ** 2

    def density(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        return abs(omega) ** 2 + np.abs(potential(t)) + logarithmic(t).imag ** 2

    return quantum, density


def integral_ratio(solution, omega, potential):
    quantum, density = pieces(solution, omega, potential)
    samples = np.linspace(*WINDOW, 4001)
    values = quantum(samples)
    breaks = [WINDOW[0]]
    for index in range(len(samples) - 1):
        if values[index] * values[index + 1] < 0.0:
            breaks.append(float(brentq(lambda t: float(quantum(np.array([t]))[0]),
                                       samples[index], samples[index + 1],
                                       xtol=1e-14, rtol=8.9e-16)))
    breaks.append(WINDOW[1])

    def piecewise(function):
        total = 0.0
        for a, b in zip(breaks[:-1], breaks[1:]):
            if b - a < 1e-13:
                continue
            total += quad(function, a, b, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
        return total

    weight = lambda t: float(np.exp(-t * t))
    upper = piecewise(lambda t: weight(t) * abs(quantum(np.array([t]))[0]))
    lower = piecewise(lambda t: weight(t) * density(np.array([t]))[0])
    return upper / lower, len(breaks) - 2


def weighted_median(values, weights):
    """Quantile discreto: ordina per valore, cumula i pesi, prende la meta'."""
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, 0.5 * cumulative[-1])])


def median_ratio(solution, omega, potential, count=4001):
    quantum, density = pieces(solution, omega, potential)
    grid = np.linspace(*WINDOW, count)
    weights = np.exp(-grid * grid)
    return (weighted_median(np.abs(quantum(grid)), weights)
            / weighted_median(density(grid), weights))


def main() -> None:
    print("Gate D — il diagnostico predice l'errore della WKB3?\n")

    approximate, lambda2, lambda3 = wkb3()
    print("1. getto centrale, invariato per costruzione (le bump sono fuori da [-1,1])")
    print(f"   Lambda2 = {lambda2:.12f}   Lambda3 = {lambda3:.12f}")
    print(f"   omega_WKB3 = {approximate.real:.12f}{approximate.imag:+.12f}i")

    family = [(1.6, 0.2), (1.6, 0.4), (2.0, 0.2), (2.0, 0.4), (2.4, 0.2), (2.4, 0.4)]
    etas = (-0.05, -0.02, 0.0, 0.02, 0.05)

    print("\n2. errore spettrale indipendente contro i due diagnostici")
    print("   centro largh.   eta     omega esatta                  err.WKB3    "
          "E_integrale    E_mediana   zeri")
    rows = []
    for center, width in family:
        guess = GUESS
        for eta in etas:
            omega, solution = solve(center, width, eta, guess=guess)
            guess = omega
            potential = make_potential(center, width, eta)
            error = abs(omega - approximate) / abs(omega)
            integral, zeros = integral_ratio(solution, omega, potential)
            median = median_ratio(solution, omega, potential)
            rows.append((center, width, eta, omega, error, integral, median))
            print(f"   {center:5.1f} {width:5.2f}  {eta:+.2f}  "
                  f"{omega.real:.9f}{omega.imag:+.9f}i  {error:.6e}  "
                  f"{integral:.9f}  {median:.9f}   {zeros}")

    print("\n3. E segue l'errore?  (nessun fit: si confrontano gli andamenti)")
    error = np.array([r[4] for r in rows])
    integral = np.array([r[5] for r in rows])
    median = np.array([r[6] for r in rows])
    for name, series in (("E_integrale", integral), ("E_mediana", median)):
        centred_e = error - error.mean()
        centred_s = series - series.mean()
        correlation = float(centred_e @ centred_s /
                            np.sqrt((centred_e @ centred_e) * (centred_s @ centred_s)))
        spread = float(np.ptp(series) / series.mean())
        print(f"   {name:12s}  correlazione di Pearson con l'errore = {correlation:+.4f}"
              f"   escursione relativa = {spread:.2e}")
    print(f"   {'errore WKB3':12s}  escursione relativa = {float(np.ptp(error)/error.mean()):.2e}")

    print("\n4. lettura")
    print("   Il getto e' identico in tutte le righe, quindi omega_WKB3 e' una")
    print("   costante e l'errore varia SOLO perche' varia il bersaglio.")
    print("   Se E fosse un predittore dell'errore la correlazione sarebbe ~1 e le")
    print("   escursioni comparabili.  I numeri sopra decidono.")


if __name__ == "__main__":
    main()
