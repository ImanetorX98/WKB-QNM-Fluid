#!/usr/bin/env python3
"""Numeratore della condizione di solvibilita' di Vaidya, senza range dinamico.

MOTIVO.  Il numeratore N(L_+) = int_a^{L_+} I dr - F(L_+) non era indipendente
da L_+: derivava del 35% fra L_+=40 e L_+=90.  La diagnosi corrente e' che
NESSUNA delle cause sospettate era quella vera.  In ordine di scoperta:

  1. quadratura disallineata di una cella (`integrand[i:j]` esclude j, `F[j]-F[i]`
     no) — errore complesso 1.0e-3, diagnosi di Codex, vedi
     `repro_numerator_discrepancy.py`;
  2. `np.gradient` per S = 2(rZ')' — differenze finite del secondo ordine,
     errore (2 om h)^2/6 = 2e-7 costante.  `source_integrand` gia' evitava il
     problema ricavando S dall'ODE;
  3. trapezio invece di Simpson — 1e-7, invisibile finche' (1) e (2) dominavano;
  4. **la condizione al bordo dell'orizzonte**, che era il limite vero.

Il punto (4) merita una riga.  psi = e^{-i om r_*} al primo ordine, imposta a
r = 2 + delta, sbaglia di O(delta) e semina il modo ENTRANTE all'infinito.  Che
sia lui si vede dal tasso: |g/C - 1| decade esattamente come |e^{-2 i om r_*}|.
Che non sia la frequenza si vede perturbando omega di 1e-9 senza alcun effetto.
Ridurre delta migliora la BC ma peggiora la cancellazione vicino all'orizzonte
(f = 1-2/r diventa 1e-10): l'ottimo era un misero 2.4e-4 a delta = 1e-10.

CURA — questo modulo.  Due riscritture che tolgono l'esponenziale invece di
comprarlo con le cifre, piu' Frobenius al posto della BC troncata:

    dentro   psi = e^{-i om r_*} h(r)   =>  **Z = h**,        h e' O(1)
    fuori    psi = e^{+i om r_*} u(r) g(r)                     g -> C

Nessuna delle due variabili cresce, quindi la doppia precisione basta e mpmath
non serve.  La BC di Frobenius (due termini, errore O(delta^3)) permette
delta = 1e-5, dove non c'e' cancellazione.

RISULTATO.  N = 40.2951 + 4.1600 i, con dispersione in L_+ di 7.6e-6 — cinque
ordini meglio del 10% precedente.  Stabile a 1e-5 rispetto a delta, punto di
raccordo, risoluzione e numero di termini della serie.

AVVERTENZA.  Il limite inferiore a = 4 resta un taglio arbitrario: il bordo
interno della regolarizzazione di Leung non e' trattato qui.  N != 0 dice che la
proiezione risonante non si annulla; non e' ancora un coefficiente fisico.

Uso: python3.13 calculations/vaidya_numerator_factored.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from asymptotic_series import Series, antiderivative_series, integrand_series  # noqa: E402
from leaver_qnm import leaver_qnm  # noqa: E402
from outgoing_asymptotics import series_coefficients  # noqa: E402
from vaidya_solvability import tortoise  # noqa: E402

ELL, SPIN, TERMS = 2, 0, 20


class Asymptotic:
    """u, u', u'' esatte dai coefficienti, piu' i coefficienti dell'ODE per g."""

    def __init__(self, ell: int, spin: int, omega: complex, terms: int) -> None:
        self.ell, self.spin, self.omega = ell, spin, omega
        self.u = Series(0, series_coefficients(ell, spin, omega, terms))
        self.u1 = self.u.derivative()
        self.u2 = self.u1.derivative()

    def at(self, radius):
        r = np.atleast_1d(np.asarray(radius, dtype=float))
        u, u1, u2 = (s.evaluate(r) for s in (self.u, self.u1, self.u2))
        f, df = 1.0 - 2.0 / r, 2.0 / r**2
        w = u1 / u
        dw = u2 / u - w * w
        potential = self.ell * (self.ell + 1) / r**2 + 2.0 * (1 - self.spin**2) / r**3
        k = 2j * self.omega
        defect = (df * w + f * dw) / f + k * w / f + w * w - potential / f
        drift = df / f + k / f + 2.0 * w
        return u, u1, u2, f, df, w, defect, drift


def horizon_start(ell: int, spin: int, omega: complex, delta: float):
    """h(2+delta), h'(2+delta) da due termini di Frobenius: errore O(delta^3)."""
    u0 = ell * (ell + 1) / 4.0 + 2.0 * (1 - spin**2) / 8.0
    u1 = -ell * (ell + 1) / 4.0 - 3.0 * (1 - spin**2) / 8.0
    c1 = u0 / (0.5 - 2j * omega)
    c2 = (c1 * (0.5 + u0) + u1) / (2.0 - 4j * omega)
    return 1.0 + c1 * delta + c2 * delta**2, c1 + 2.0 * c2 * delta


def inner_solve(ell, spin, omega, delta, lower, r_match, points):
    """h su [2+delta, r_match] con `lower` esattamente sulla griglia.

    Z = h: l'integrando e' I = e^{-2 i om r_*} h * 2 (h' + r h'').
    """

    def curvature(r, h, dh):
        f = 1.0 - 2.0 / r
        potential = ell * (ell + 1) / r**2 + 2.0 * (1 - spin**2) / r**3
        return (potential * h - (2.0 / r**2 - 2j * omega) * dh) / f

    def rhs(r, y):
        return np.array([y[1], curvature(r, y[0], y[1])], dtype=complex)

    start = 2.0 + delta
    grid = np.concatenate(
        [np.linspace(start, lower, points // 4, endpoint=False),
         np.linspace(lower, r_match, points)]
    )
    solution = solve_ivp(
        rhs,
        (start, r_match),
        np.array(horizon_start(ell, spin, omega, delta), dtype=complex),
        t_eval=grid, method="DOP853", rtol=1.0e-13, atol=1.0e-16,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    h, dh = solution.y
    integrand = (np.exp(-2j * omega * tortoise(grid)) * h
                 * 2.0 * (dh + grid * curvature(grid, h, dh)))
    return grid, h, dh, integrand


def outer_solve(asy: Asymptotic, r_match: float, r_far: float, g0, dg0, points: int):
    """g'' + drift g' + defect g = 0, integrata verso l'esterno.  g -> C."""

    def rhs(r, y):
        _, _, _, _, _, _, defect, drift = asy.at(r)
        return np.array([y[1], -drift[0] * y[1] - defect[0] * y[0]], dtype=complex)

    grid = np.linspace(r_match, r_far, points)
    solution = solve_ivp(rhs, (r_match, r_far), np.array([g0, dg0], dtype=complex),
                         t_eval=grid, method="DOP853", rtol=1.0e-13, atol=1.0e-16)
    if not solution.success:
        raise RuntimeError(solution.message)
    return grid, solution.y[0], solution.y[1]


def outer_integrand(asy: Asymptotic, r, g, dg):
    """I = mu Z S nella regione esterna, da g e dalle serie esatte."""
    u, u1, u2, f, df, _, defect, drift = asy.at(r)
    k = 2j * asy.omega
    ddg = -drift * dg - defect * g
    bracket = k * u * g / f + u1 * g + u * dg
    dbracket = (k * ((u1 * g + u * dg) / f - u * g * df / f**2)
                + u2 * g + 2.0 * u1 * dg + u * ddg)
    return (2.0 * u * g * np.exp(k * tortoise(r))
            * (bracket + r * (k * bracket / f + dbracket)))


def numerator(
    ell: int = ELL,
    spin: int = SPIN,
    overtone: int = 0,
    delta: float = 1.0e-5,
    lower: float = 4.0,
    r_match: float = 25.0,
    r_far: float = 95.0,
    inner_points: int = 80001,
    outer_points: int = 200001,
    terms: int = TERMS,
    uppers: tuple[float, ...] = (40.0, 50.0, 60.0, 70.0, 80.0, 90.0),
) -> dict:
    """N(L_+) su piu' L_+, piu' la dispersione che ne misura l'indipendenza."""
    omega = leaver_qnm(ell, overtone, spin)
    asy = Asymptotic(ell, spin, omega, terms)
    r_in, h, dh, inner = inner_solve(ell, spin, omega, delta, lower, r_match, inner_points)

    u, _, _, f, _, w, _, _ = asy.at(r_match)
    phase = np.exp(-2j * omega * tortoise(r_match))
    g0 = phase * h[-1] / u[0]
    dpsi = -1j * omega * h[-1] / f[0] + dh[-1]
    dg0 = phase * (dpsi - h[-1] * (1j * omega / f[0] + w[0])) / u[0]

    r_out, g, dg = outer_solve(asy, r_match, r_far, g0, dg0, outer_points)
    outer = outer_integrand(asy, r_out, g, dg)
    constant = g[np.searchsorted(r_out, r_far - 10.0)]
    primitive = (constant**2
                 * antiderivative_series(integrand_series(ell, spin, omega, terms),
                                         omega, terms, terms).evaluate(r_out)
                 * np.exp(2j * omega * tortoise(r_out)))

    start = int(np.searchsorted(r_in, lower))
    base = simpson(inner[start:], x=r_in[start:])
    values = {}
    for upper in uppers:
        j = int(np.searchsorted(r_out, upper))
        values[upper] = complex(base + simpson(outer[: j + 1], x=r_out[: j + 1]) - primitive[j])

    reference = values[max(values)]
    spread = max(abs(v - reference) for v in values.values()) / abs(reference)
    return {
        "omega": omega,
        "C": complex(constant),
        "N": reference,
        "values": values,
        "spread": float(spread),
        "flatness": float(max(abs(g[np.searchsorted(r_out, r_far - 15.0)] / constant - 1.0), 0.0)),
    }


def main() -> None:
    base = numerator()
    print("Numeratore della condizione di solvibilita', formulazione fattorizzata")
    print(f"  M*omega = {base['omega'].real:.9f}{base['omega'].imag:+.9f}i")
    print(f"  C = g(inf) = {base['C'].real:+.9f}{base['C'].imag:+.9f}i")
    print(f"  N = {base['N'].real:+.8f}{base['N'].imag:+.8f}i")
    print()
    print("=== N(L_+): l'indipendenza dal taglio esterno ===")
    for upper, value in base["values"].items():
        print(f"   L+ = {upper:4.0f}    {value.real:+.8f}{value.imag:+.8f}i")
    print(f"   dispersione = {base['spread']:.2e}")
    print()
    print("=== robustezza: ogni parametro separatamente ===")
    checks = (
        ("delta", "delta", (1e-3, 1e-4, 1e-5, 1e-6, 1e-8)),
        ("r_match", "r_match", (15.0, 20.0, 25.0, 30.0, 35.0)),
        ("inner_points", "punti interni", (40001, 80001, 160001)),
        ("terms", "termini serie", (12, 16, 20, 24)),
    )
    for key, label, options in checks:
        results = [numerator(**{key: option}) for option in options]
        centre = np.mean([row["N"] for row in results])
        drift = max(abs(row["N"] - centre) for row in results) / abs(centre)
        print(f"  {label:<15} N = {results[-1]['N'].real:.6f}"
              f"{results[-1]['N'].imag:+.6f}i   variazione = {drift:.1e}")


if __name__ == "__main__":
    main()
