#!/usr/bin/env python3
"""Controllo F del brief: residuo, norma generalizzata e matching.

TESI DA VERIFICARE (`research/matching_invariance_proof_2026-09-12.md`).
Con normalizzazione psi(a)=1 e residuo **non normalizzato**

    F(omega) = psi'(b) - D_plus(omega) psi(b),

al QNM vale  N = -psi(b) F_omega,  dove N e' il denominatore della norma
generalizzata,  N = 2 omega int_a^b psi^2 dx + i[psi(a)^2 + psi(b)^2].

La normalizzazione conta: un Wronskiano **normalizzato** porta fattori diversi,
e il confronto va fatto sul residuo nudo.

F_omega e' calcolata in due modi:
  * **equazione variazionale**: phi = partial_omega psi risolve
    phi'' = (V - omega^2) phi - 2 omega psi con phi(a)=0, phi'(a)=-i, perche'
    psi(a)=1 e' fissata e psi'(a) = -i omega;
  * **differenze centrali in omega**, lungo le direzioni reale e immaginaria,
    con tre passi.

Si verifica poi che la stima locale delta omega = -delta F / F_omega descriva
l'effetto di un raffinamento numerico: il residuo da solo **non** e' l'errore
spettrale.

Infine, sullo spostamento dei bordi: si distingue l'invarianza del **rapporto**
B/N — che regge perche' il fattore di normalizzazione si cancella — da quella di
B e N **separatamente**, che richiede la stessa normalizzazione globale e non
vale se si rinormalizza a ogni dominio.

Uso: python3.13 calculations/claude_compact_matching.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp
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


def march(omega, eta=0.0, left=LEFT, right=RIGHT, tol=1.0e-14, dense=False):
    """psi con psi(a)=1, psi'(a)=-i omega."""
    return solve_ivp(lambda x, y: [y[1], (potential(x, eta) - omega**2) * y[0]],
                     (left, right), [1.0 + 0j, -1j * omega], method="DOP853",
                     rtol=tol, atol=tol / 100.0, max_step=0.035, dense_output=dense)


def residual(omega, eta=0.0, left=LEFT, right=RIGHT, tol=1.0e-14):
    """F = psi'(b) - i omega psi(b),  non normalizzato."""
    out = march(omega, eta, left, right, tol)
    return out.y[1, -1] - 1j * omega * out.y[0, -1], out.y[0, -1]


def resonance(eta=0.0, guess=GUESS, left=LEFT, right=RIGHT, tol=1.0e-14):
    def real_pair(v):
        value, _ = residual(complex(*v), eta, left, right, tol)
        return [value.real, value.imag]

    found = root(real_pair, [guess.real, guess.imag], tol=1.0e-13)
    return complex(*found.x)


def residual_derivative_variational(omega, eta=0.0, left=LEFT, right=RIGHT,
                                    tol=1.0e-14):
    """F_omega dall'equazione variazionale, senza differenze finite."""
    def rhs(x, y):
        psi, dpsi, phi, dphi = y
        curvature = potential(x, eta) - omega**2
        return [dpsi, curvature * psi, dphi, curvature * phi - 2.0 * omega * psi]

    out = solve_ivp(rhs, (left, right), [1.0 + 0j, -1j * omega, 0.0 + 0j, -1j],
                    method="DOP853", rtol=tol, atol=tol / 100.0, max_step=0.035)
    if not out.success:
        raise RuntimeError(out.message)
    psi, phi, dphi = out.y[0, -1], out.y[2, -1], out.y[3, -1]
    return dphi - 1j * omega * phi - 1j * psi, psi


def generalized_norm(omega, eta=0.0, left=LEFT, right=RIGHT, tol=1.0e-14):
    """N = 2 omega int psi^2 + i[psi(a)^2 + psi(b)^2],  stessa normalizzazione."""
    out = march(omega, eta, left, right, tol, dense=True)
    square = lambda t: out.sol(np.array([t]))[0][0] ** 2
    real = quad(lambda t: float(square(t).real), left, right, limit=900,
                epsabs=1e-14, epsrel=1e-14)[0]
    imag = quad(lambda t: float(square(t).imag), left, right, limit=900,
                epsabs=1e-14, epsrel=1e-14)[0]
    edges = square(left) + square(right)
    return 2.0 * omega * (real + 1j * imag) + 1j * edges, out


def source_overlap(out, left=LEFT, right=RIGHT):
    """B = int b psi^2, stessa normalizzazione psi(a)=1."""
    square = lambda t: out.sol(np.array([t]))[0][0] ** 2
    lo, hi = max(left, 1.6), min(right, 2.4)
    real = quad(lambda t: float(bump(np.array([t]), 2.0, 0.4)[0] * square(t).real),
                lo, hi, limit=400, epsabs=1e-14, epsrel=1e-14)[0]
    imag = quad(lambda t: float(bump(np.array([t]), 2.0, 0.4)[0] * square(t).imag),
                lo, hi, limit=400, epsabs=1e-14, epsrel=1e-14)[0]
    return real + 1j * imag


def main() -> None:
    print("Controllo F — residuo, norma e matching\n")
    omega = resonance()
    print(f"omega = {omega.real:.12f}{omega.imag:+.12f}i")

    print("\n1. F_omega: equazione variazionale contro differenze centrali")
    variational, psi_right = residual_derivative_variational(omega)
    print(f"   variazionale      {variational.real:+.12f}{variational.imag:+.12f}i")
    for step in (1.0e-5, 1.0e-6, 1.0e-7):
        along_real = ((residual(omega + step)[0] - residual(omega - step)[0])
                      / (2.0 * step))
        along_imag = ((residual(omega + 1j * step)[0] - residual(omega - 1j * step)[0])
                      / (2.0j * step))
        print(f"   passo {step:.0e}   reale {along_real.real:+.12f}{along_real.imag:+.12f}i"
              f"   |diff| {abs(along_real - variational):.2e}"
              f"   |imm.-var| {abs(along_imag - variational):.2e}")

    print("\n2. N = -psi(b) F_omega  al QNM?")
    norm, solution = generalized_norm(omega)
    predicted = -psi_right * variational
    print(f"   N diretta         {norm.real:+.12f}{norm.imag:+.12f}i")
    print(f"   -psi(b) F_omega   {predicted.real:+.12f}{predicted.imag:+.12f}i")
    print(f"   scarto relativo   {abs(predicted / norm - 1):.3e}")
    print(f"   psi(b) = {psi_right.real:+.9f}{psi_right.imag:+.9f}i   "
          "(il fattore che un Wronskiano normalizzato nasconderebbe)")

    print("\n3. delta omega = -delta F / F_omega contro un raffinamento vero")
    coarse = resonance(tol=1.0e-9)
    fine = resonance(tol=1.0e-14)
    value_at_coarse, _ = residual(coarse, tol=1.0e-14)
    estimate = -value_at_coarse / variational
    print(f"   omega a tol 1e-9   {coarse.real:.12f}{coarse.imag:+.12f}i")
    print(f"   omega a tol 1e-14  {fine.real:.12f}{fine.imag:+.12f}i")
    print(f"   spostamento vero   {abs(fine - coarse):.3e}")
    print(f"   stima -deltaF/F_om {abs(estimate):.3e}")
    print(f"   |residuo| nudo     {abs(value_at_coarse):.3e}   <- da solo NON e'"
          " l'errore spettrale")

    print("\n4. spostamento dei bordi: il rapporto B/N contro B e N separati")
    print("   [a, b]         B                      N                      B/N")
    reference = None
    for left, right in ((-1.0, 2.4), (-1.5, 2.4), (-1.0, 3.0), (-2.0, 3.5)):
        omega_local = resonance(left=left, right=right)
        norm_local, out_local = generalized_norm(omega_local, left=left, right=right)
        overlap = source_overlap(out_local, left, right)
        ratio = overlap / norm_local
        if reference is None:
            reference = (overlap, norm_local, ratio)
        print(f"   [{left:5.1f},{right:4.1f}]  {abs(overlap):.10f}  {abs(norm_local):.10f}  "
              f"{ratio.real:+.12f}{ratio.imag:+.12f}i")
    print("\n   scarti dal primo dominio:")
    for left, right in ((-1.5, 2.4), (-1.0, 3.0), (-2.0, 3.5)):
        omega_local = resonance(left=left, right=right)
        norm_local, out_local = generalized_norm(omega_local, left=left, right=right)
        overlap = source_overlap(out_local, left, right)
        print(f"   [{left:5.1f},{right:4.1f}]  |B/B0-1| = {abs(overlap/reference[0]-1):.2e}"
              f"   |N/N0-1| = {abs(norm_local/reference[1]-1):.2e}"
              f"   |(B/N)/(B0/N0)-1| = {abs((overlap/norm_local)/reference[2]-1):.2e}")


if __name__ == "__main__":
    main()
