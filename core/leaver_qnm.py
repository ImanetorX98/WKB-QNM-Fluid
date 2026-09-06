#!/usr/bin/env python3
"""QNM di Schwarzschild col metodo delle frazioni continue di Leaver (1985).

A differenza della WKB di barriera, che usa solo il getto del potenziale nel
massimo, questo metodo impone entrambe le condizioni al contorno dentro
l'ansatz e non approssima la barriera: la precisione e' limitata solo dalla
profondita' della frazione continua.

Ansatz (unita' M=1, orizzonte in r=2, z=(r-2)/r):

    psi = e^{i omega r} (r-2)^{-2 i omega} r^{4 i omega} sum_n a_n z^n,

da cui la ricorrenza a tre termini alpha_n a_(n+1) + beta_n a_n + gamma_n a_(n-1) = 0
con (derivata simbolicamente, non trascritta):

    alpha_n = (n+1)(n+1-4 i omega)
    beta_n  = -l(l+1) - 2n^2 - 2n + 16 i n omega + 32 omega^2 + 8 i omega + s^2 - 1
    gamma_n = (n - 4 i omega)^2 - s^2

La condizione di QNM e' che gli a_n siano la soluzione minimale, equivalente
all'annullarsi della frazione continua.  Per l'overtone n si usa l'equazione
invertita n volte: le radici ci sono tutte in ogni inversione, ma solo nella
p-esima la p-esima radice e' numericamente stabile.

Uso: python3.13 leaver_qnm.py --ell 2 --spin 2 --n 0
"""

from __future__ import annotations

import argparse
import cmath

# Il default e' volutamente profondo: la convergenza peggiora quando
# |Im omega| cresce rispetto a |Re omega|, cioe' proprio agli overtone alti.
DEFAULT_DEPTH = 3000


def coefficients(n: int, omega: complex, ell: int, spin: int) -> tuple[complex, complex, complex]:
    """(alpha_n, beta_n, gamma_n) della ricorrenza a tre termini."""
    large_l = ell * (ell + 1)
    alpha = (n + 1) * (n + 1 - 4j * omega)
    beta = (
        -large_l
        - 2 * n**2
        - 2 * n
        + 16j * n * omega
        + 32 * omega**2
        + 8j * omega
        + spin**2
        - 1
    )
    gamma = (n - 4j * omega) ** 2 - spin**2
    return alpha, beta, gamma


def _tail(
    start: int, omega: complex, ell: int, spin: int, depth: int
) -> complex:
    """Coda u_start = a_(start+1)/a_start della soluzione minimale.

    Calcolata all'indietro da `depth`; il valore iniziale nullo e' la scelta
    conservativa (troncamento netto).  La convergenza si controlla variando
    `depth`, come raccomandato da Leaver.
    """
    u = 0.0 + 0.0j
    for n in range(depth, start, -1):
        alpha, beta, gamma = coefficients(n, omega, ell, spin)
        denominator = beta + alpha * u
        if denominator == 0:
            return complex("nan")
        u = -gamma / denominator
    return u


def continued_fraction(
    omega: complex, ell: int, spin: int, overtone: int = 0, depth: int = DEFAULT_DEPTH
) -> complex:
    """Residuo dell'equazione invertita `overtone` volte; zero sui QNM."""
    alpha0, beta0, _ = coefficients(0, omega, ell, spin)
    if alpha0 == 0:
        return complex("nan")

    # Testa: u_0 dalla condizione a_(-1)=0, poi ricorrenza in avanti.
    head = -beta0 / alpha0
    for n in range(1, overtone):
        alpha, beta, gamma = coefficients(n, omega, ell, spin)
        if head == 0:
            return complex("nan")
        head = -(beta + gamma / head) / alpha

    alpha, beta, gamma = coefficients(overtone, omega, ell, spin)
    residual = beta + alpha * _tail(overtone, omega, ell, spin, depth)
    if overtone > 0:
        if head == 0:
            return complex("nan")
        residual += gamma / head
    return residual


def _muller(
    function, guesses: tuple[complex, complex, complex], tolerance: float, max_steps: int
) -> complex:
    """Metodo di Muller: root finder complesso senza derivate."""
    x0, x1, x2 = guesses
    f0, f1, f2 = function(x0), function(x1), function(x2)
    for _ in range(max_steps):
        h0, h1 = x1 - x0, x2 - x1
        if h0 == 0 or h1 == 0:
            break
        d0, d1 = (f1 - f0) / h0, (f2 - f1) / h1
        if h1 + h0 == 0:
            break
        a = (d1 - d0) / (h1 + h0)
        b = a * h1 + d1
        radicand = b * b - 4 * a * f2
        denominator = max(
            (b + cmath.sqrt(radicand), b - cmath.sqrt(radicand)), key=abs
        )
        if denominator == 0:
            break
        step = -2 * f2 / denominator
        x0, x1, x2 = x1, x2, x2 + step
        f0, f1, f2 = f1, f2, function(x2)
        if abs(step) < tolerance:
            break
    return x2


def leaver_qnm(
    ell: int,
    overtone: int = 0,
    spin: int = 2,
    guess: complex | None = None,
    depth: int = DEFAULT_DEPTH,
    tolerance: float = 1.0e-12,
    max_steps: int = 200,
) -> complex:
    """M*omega del QNM, per ricerca di radice sull'equazione invertita."""
    if guess is None:
        from schwarzschild_wkb import qnm_wkb

        guess = qnm_wkb(ell, overtone, spin, 3).omega_M

    def residual(omega: complex) -> complex:
        return continued_fraction(omega, ell, spin, overtone, depth)

    scale = 1.0e-3 * abs(guess)
    root = _muller(
        residual,
        (guess * (1 + 1e-3), guess * (1 - 1e-3), guess + 1j * scale),
        tolerance,
        max_steps,
    )
    if root.real < 0:
        root = -root.conjugate()
    return root


def convergence_report(
    ell: int, overtone: int, spin: int, depths: tuple[int, ...]
) -> list[tuple[int, complex]]:
    """Stessa radice a profondita' crescenti: il test di convergenza di Leaver."""
    return [(d, leaver_qnm(ell, overtone, spin, depth=d)) for d in depths]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ell", type=int, default=2)
    parser.add_argument("--spin", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--n", type=int, default=0, dest="overtone")
    parser.add_argument("--depth", type=int, default=DEFAULT_DEPTH)
    parser.add_argument("--convergence", action="store_true")
    args = parser.parse_args()

    if args.convergence:
        print(" profondita        M omega")
        for depth, omega in convergence_report(
            args.ell, args.overtone, args.spin, (200, 500, 1000, 2000, 4000)
        ):
            print(f"  {depth:8d}   {omega.real:.10f} {omega.imag:+.10f} i")
        return

    omega = leaver_qnm(args.ell, args.overtone, args.spin, depth=args.depth)
    print(f"s={args.spin}, ell={args.ell}, n={args.overtone}")
    print(f"M omega = {omega.real:.10f} {omega.imag:+.10f} i")


if __name__ == "__main__":
    main()
