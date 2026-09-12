#!/usr/bin/env python3
"""Controllo G del brief: secondo ordine e zeri mobili del diagnostico.

TESI DA VERIFICARE (`research/second_order_core_response_2026-09-12.md` §4).
Per U(eta) = int_J w |Q(x,eta)| dx la prima derivata non ha termini agli zeri,
perche' |Q| vi si annulla, ma la **seconda** si', perche' l'integrando della
prima ha un salto:

    U1 = int_J w sgn(Q) Q1 dx,
    U2 = int_J w sgn(Q) Q2 dx  +  2 sum_j w(x_j) Q1(x_j)^2 / |Q_x(x_j)| .

Il termine aggiuntivo e' non negativo.  Se si omette, la seconda derivata e'
sbagliata — e il controesempio esatto del §4 lo mostra senza alcuna
approssimazione QNM.

STRUTTURA (i sei punti del brief)
  1. Q = x - eta, w = 1, J = [-1,1]:  U = 1 + eta^2, dunque U2 = 2 mentre
     Q2 = 0.  **L'intero contributo viene dallo zero mobile.**  Omettendolo si
     ottiene 0.
  2. d = omega' ed e = omega'' sul benchmark esterno (centro 2, largh. 0.4).
     Sono DERIVATE: omega = omega0 + eta d + eta^2 e/2, non coefficienti.
  3. j' + 2 z j = -2 - 2 k^2, j(a) = 0, piu' la forma integrale come controllo.
  4. Q1, Q2, zeri semplici di Q0 in J, Q_x agli zeri da derivate controllate
     (non differenze grossolane), e x_j' = -Q1/Q_x seguendo gli stessi zeri.
  5. U2 ed E2, contro [E(eta) - 2E(0) + E(-eta)]/eta^2 con quadratura spezzata
     sugli zeri di CIASCUN eta.  Termine integrale e contributo degli zeri
     riportati separatamente.
  6. Controllo negativo: la stessa cosa omettendo il termine degli zeri.

Benchmark compatto: q = omega^2 - V0, quindi q_omega = 2 omega, q_omegaomega = 2,
D_minus = -i omega, k(a) = -i, j(a) = 0.  Per Kerr non si puo' porre
q_omegaomega = 2: servirebbe la seconda derivata dell'autovalore sferoidale.

Uso: python3.13 calculations/claude_compact_second_order_response.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, root

LEFT, RIGHT = -1.0, 2.4
WINDOW = (-0.8, 0.8)
HEIGHT = 10.0
GUESS = 3.423295488 - 0.607387592j
CENTER, WIDTH = 2.0, 0.4


def bump(x, center=0.0, width=1.0):
    y = (np.asarray(x, dtype=float) - center) / width
    out = np.zeros_like(y, dtype=float)
    inside = np.abs(y) < 1.0
    out[inside] = np.exp(1.0 - 1.0 / (1.0 - y[inside] ** 2))
    return out


def potential(x, eta=0.0):
    return HEIGHT * bump(x) + eta * bump(x, CENTER, WIDTH)


def march(omega, eta=0.0, tol=1.0e-14):
    return solve_ivp(lambda x, y: [y[1], (potential(x, eta) - omega**2) * y[0]],
                     (LEFT, RIGHT), [1.0 + 0j, -1j * omega], method="DOP853",
                     rtol=tol, atol=tol / 100.0, max_step=0.035, dense_output=True)


def resonance(eta=0.0, guess=GUESS, tol=1.0e-14):
    def real_pair(v):
        omega = complex(*v)
        out = march(omega, eta, tol)
        defect = out.y[1, -1] - 1j * omega * out.y[0, -1]
        return [defect.real, defect.imag]

    found = root(real_pair, [guess.real, guess.imag], tol=1.0e-13)
    omega = complex(*found.x)
    return omega, march(omega, eta, tol)


# --------------------------------------------------------------------------
# 1. controesempio esatto
# --------------------------------------------------------------------------
def moving_zero_toy():
    """U(eta) = int_-1^1 |x - eta| dx = 1 + eta^2.  Q1 = -1, Q2 = 0, Q_x = 1."""
    exact = 2.0
    integral_term = 0.0                       # Q2 = 0 identicamente
    zero_term = 2.0 * 1.0 * (-1.0) ** 2 / 1.0  # 2 w(x_j) Q1^2 / |Q_x|
    finite = None
    step = 1.0e-3
    values = [quad(lambda t: abs(t - s), -1.0, s)[0] + quad(lambda t: abs(t - s), s, 1.0)[0]
              for s in (step, 0.0, -step)]
    finite = (values[0] - 2.0 * values[1] + values[2]) / step**2
    return exact, integral_term, zero_term, finite


# --------------------------------------------------------------------------
# 2. d ed e
# --------------------------------------------------------------------------
def spectral_derivatives(steps=(1.0e-2, 3.0e-3, 1.0e-3), tol=1.0e-14):
    base, _ = resonance(tol=tol)
    rows = []
    for step in steps:
        plus, _ = resonance(+step, guess=base, tol=tol)
        minus, _ = resonance(-step, guess=base, tol=tol)
        first = (plus - minus) / (2.0 * step)
        second = (plus - 2.0 * base + minus) / step**2
        rows.append((step, first, second))
    return base, rows


# --------------------------------------------------------------------------
# 3. k e j
# --------------------------------------------------------------------------
def kernels(solution, omega, tol=1.0e-13):
    """k e j dalle rispettive ODE, integrate insieme al profilo."""
    def rhs(x, y):
        values = solution.sol(np.array([x]))
        z = values[1][0] / values[0][0]
        k, j = y
        return [-2.0 * z * k - 2.0 * omega, -2.0 * z * j - 2.0 - 2.0 * k * k]

    out = solve_ivp(rhs, (LEFT, WINDOW[1] + 0.05), [-1j, 0.0 + 0j], method="DOP853",
                    rtol=tol, atol=tol / 100.0, dense_output=True)
    if not out.success:
        raise RuntimeError(out.message)
    return (lambda x: out.sol(np.atleast_1d(x))[0],
            lambda x: out.sol(np.atleast_1d(x))[1])


def kernel_integral_check(solution, omega, kernel, points):
    """j(x) = psi^-2 [ - int_a^x (2 + 2 k^2) psi^2 dt ],  con j(a)=0."""
    out = []
    for point in points:
        real = quad(lambda t: float(((2.0 + 2.0 * kernel(np.array([t]))[0] ** 2)
                                     * solution.sol(np.array([t]))[0][0] ** 2).real),
                    LEFT, point, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        imag = quad(lambda t: float(((2.0 + 2.0 * kernel(np.array([t]))[0] ** 2)
                                     * solution.sol(np.array([t]))[0][0] ** 2).imag),
                    LEFT, point, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        out.append(-(real + 1j * imag) / solution.sol(np.array([point]))[0][0] ** 2)
    return np.array(out)


# --------------------------------------------------------------------------
# 4-5. Q1, Q2, zeri, U2, E2
# --------------------------------------------------------------------------
def field(solution, omega, eta=0.0):
    def logarithmic(t):
        values = solution.sol(np.atleast_1d(t))
        return values[1] / values[0]

    def quantum(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        return (omega * omega - potential(t, eta)).real - logarithmic(t).imag ** 2

    def density(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        return abs(omega) ** 2 + np.abs(potential(t, eta)) + logarithmic(t).imag ** 2

    return logarithmic, quantum, density


def quantum_slope(solution, omega, t):
    """Q_x da derivate controllate: p_x = -Im(z^2) - Im(omega^2), Q_x = -V0' - 2 p p_x."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    values = solution.sol(t)
    z = values[1] / values[0]
    p = z.imag
    slope_p = -(z * z).imag - (omega * omega).imag
    step = 1.0e-6
    slope_v = (HEIGHT * bump(t + step) - HEIGHT * bump(t - step)) / (2.0 * step)
    return -slope_v - 2.0 * p * slope_p


def zeros_of(function, samples):
    values = function(samples)
    out = []
    for index in range(len(samples) - 1):
        if values[index] * values[index + 1] < 0.0:
            out.append(float(brentq(lambda t: float(function(np.array([t]))[0]),
                                    samples[index], samples[index + 1],
                                    xtol=1e-14, rtol=8.9e-16)))
    return out


def split_integral(function, breaks):
    total = 0.0
    for a, b in zip(breaks[:-1], breaks[1:]):
        if b - a < 1e-13:
            continue
        total += quad(function, a, b, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
    return total


def indicator(solution, omega, eta=0.0):
    _, quantum, density = field(solution, omega, eta)
    breaks = [WINDOW[0]] + zeros_of(quantum, np.linspace(*WINDOW, 4001)) + [WINDOW[1]]
    weight = lambda t: float(np.exp(-t * t))
    upper = split_integral(lambda t: weight(t) * abs(quantum(np.array([t]))[0]), breaks)
    lower = split_integral(lambda t: weight(t) * density(np.array([t]))[0],
                           [WINDOW[0], WINDOW[1]])
    return upper / lower, upper, lower, breaks


def main() -> None:
    print("Controllo G — secondo ordine e zeri mobili\n")

    print("1. controesempio esatto  Q = x - eta,  w = 1,  J = [-1,1]")
    exact, integral_term, zero_term, finite = moving_zero_toy()
    print(f"   U(eta) = 1 + eta^2  ->  U2 esatta = {exact}")
    print(f"   termine integrale (Q2 = 0)        = {integral_term}")
    print(f"   termine degli zeri 2 w Q1^2/|Q_x| = {zero_term}")
    print(f"   somma = {integral_term + zero_term}   differenze centrali = {finite:.10f}")
    print(f"   OMETTENDO il termine degli zeri si otterrebbe {integral_term}: sbagliato.")

    print("\n2. d = omega' ed e = omega'' (derivate, non coefficienti)")
    omega, rows = spectral_derivatives()
    print(f"   omega0 = {omega.real:.12f}{omega.imag:+.12f}i")
    for step, first, second in rows:
        print(f"   passo {step:.0e}:  d = {first.real:+.10f}{first.imag:+.10f}i"
              f"   e = {second.real:+.8f}{second.imag:+.8f}i")
    drift = rows[-1][1]
    curvature = rows[-1][2]

    _, solution = resonance()
    kernel, second_kernel = kernels(solution, omega)

    print("\n3. j dalle due vie")
    probes = np.array([-0.8, -0.4, 0.0, 0.4, 0.8])
    from_ode = second_kernel(probes)
    from_integral = kernel_integral_check(solution, omega, kernel, probes)
    for point, value, other in zip(probes, from_ode, from_integral):
        print(f"   x={point:+.1f}   j = {value.real:+.9f}{value.imag:+.9f}i"
              f"   |ODE - integrale| = {abs(value - other):.2e}")

    print("\n4. zeri di Q0 in J, pendenza, e loro moto")
    _, quantum, density = field(solution, omega)
    roots = zeros_of(quantum, np.linspace(*WINDOW, 4001))
    slopes = quantum_slope(solution, omega, np.array(roots))
    p_at = field(solution, omega)[0](np.array(roots)).imag
    k_at = kernel(np.array(roots))
    q1_at = 2.0 * (omega * drift).real - 2.0 * p_at * (k_at * drift).imag
    print("      x_j            Q_x           Q1(x_j)      -Q1/Q_x   x_j' misurata")
    step = 1.0e-3
    moved = []
    for sign in (+1, -1):
        omega_s, solution_s = resonance(sign * step, guess=omega)
        _, quantum_s, _ = field(solution_s, omega_s, sign * step)
        moved.append(zeros_of(quantum_s, np.linspace(*WINDOW, 4001)))
    for index, (root, slope, q1) in enumerate(zip(roots, slopes, q1_at)):
        measured = (moved[0][index] - moved[1][index]) / (2.0 * step)
        print(f"   {root:+.10f}  {slope:+.8f}  {q1:+.8f}  {-q1/slope:+.8f}  {measured:+.8f}"
              f"   scarto {abs(measured + q1/slope):.2e}")

    print("\n5. U2 ed E2")
    weight = lambda t: float(np.exp(-t * t))
    value, upper, lower, breaks = indicator(solution, omega)
    p_all = lambda t: field(solution, omega)[0](t).imag
    p1 = lambda t: (kernel(t) * drift).imag
    p2 = lambda t: (kernel(t) * curvature + second_kernel(t) * drift**2).imag
    q1 = lambda t: 2.0 * (omega * drift).real - 2.0 * p_all(t) * p1(t)
    q2 = lambda t: (2.0 * (omega * curvature + drift**2).real
                    - 2.0 * p1(t) ** 2 - 2.0 * p_all(t) * p2(t))

    u1 = split_integral(lambda t: weight(t) * np.sign(quantum(np.array([t]))[0])
                        * q1(np.array([t]))[0], breaks)
    u2_integral = split_integral(lambda t: weight(t) * np.sign(quantum(np.array([t]))[0])
                                 * q2(np.array([t]))[0], breaks)
    u2_zeros = sum(2.0 * weight(root) * q1(np.array([root]))[0] ** 2 / abs(slope)
                   for root, slope in zip(roots, slopes))
    d1 = split_integral(lambda t: weight(t) * (2.0 * (np.conj(omega) * drift).real
                                               + 2.0 * p_all(np.array([t]))[0] * p1(np.array([t]))[0]),
                        [WINDOW[0], WINDOW[1]])
    d2 = split_integral(lambda t: weight(t) * (2.0 * abs(drift) ** 2
                                               + 2.0 * (np.conj(omega) * curvature).real
                                               + 2.0 * p1(np.array([t]))[0] ** 2
                                               + 2.0 * p_all(np.array([t]))[0] * p2(np.array([t]))[0]),
                        [WINDOW[0], WINDOW[1]])
    e1 = (u1 - value * d1) / lower
    e2_full = ((u2_integral + u2_zeros) - value * d2 - 2.0 * e1 * d1) / lower
    e2_naked = (u2_integral - value * d2 - 2.0 * e1 * d1) / lower
    print(f"   U1 = {u1:+.10f}      D1 = {d1:+.10f}      E1 = {e1:+.10f}")
    print(f"   U2, termine integrale     = {u2_integral:+.10f}")
    print(f"   U2, contributo degli zeri = {u2_zeros:+.10f}   <- non negativo")
    print(f"   U2 totale                 = {u2_integral + u2_zeros:+.10f}")
    print(f"   D2 = {d2:+.10f}")
    print(f"   E2 con il termine degli zeri   = {e2_full:+.10f}")
    print(f"   E2 SENZA (controllo negativo)  = {e2_naked:+.10f}")

    print("\n   contro differenze centrali seconde di E, quadratura spezzata a ogni eta")
    print("   passo eta      E2 numerica        scarto con E2      scarto senza zeri")
    for step in (1.0e-2, 3.0e-3, 1.0e-3):
        values = []
        for sign in (+1, 0, -1):
            omega_s, solution_s = resonance(sign * step, guess=omega)
            values.append(indicator(solution_s, omega_s, sign * step)[0])
        numeric = (values[0] - 2.0 * values[1] + values[2]) / step**2
        print(f"   {step:.0e}     {numeric:+.10f}     {abs(numeric/e2_full-1):.3e}"
              f"        {abs(numeric/e2_naked-1):.3e}")

    print("\n6. lettura del controllo negativo")
    print(f"   il termine degli zeri vale {u2_zeros:.6f} su un U2 totale di "
          f"{u2_integral + u2_zeros:.6f}")
    print(f"   ometterlo sposta E2 da {e2_full:.6f} a {e2_naked:.6f}, "
          f"cioe' del {abs(e2_naked/e2_full-1)*100:.1f}%")


if __name__ == "__main__":
    main()
