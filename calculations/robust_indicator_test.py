#!/usr/bin/env python3
"""Passo 1: l'indicatore di Madelung sopravvive alla robustificazione sui nodi?

Due difetti separati sono stati misurati nell'indicatore integrale

    E_M = int w |Q_M| / int w (|eps omega|^2 + eps^2 |V| + |P|^2):

  (a) a n=0 e' degenere, cioe' una funzione del solo eps (esattamente su
      Poschl-Teller, entro il 3% su Schwarzschild s=0);
  (b) a n=1 sembra portare informazione, ma l'integrale e' dominato all'89%
      dal picco di |Q_M| al nodo, e la dominanza cresce con ell.

Se (b) e' un artefatto, sostituendo l'integrale con una statistica robusta
al polo la dipendenza da ell deve collassare su quella di n=0.  Si provano
tre varianti:

    E_med      mediana pesata di |Q_M| su mediana pesata del denominatore
    E_trim     integrale troncato al 90-esimo percentile di |Q_M|
    E_excise   integrale escluso un intorno +-0.25 sigma del minimo di |psi|

Il secondo test e' quello richiesto dal protocollo dell'audit (§5.5):
validazione a **eps fissato** su una deformazione di forma indipendente.
A ell fissato lo spin s=0,1,2 cambia la barriera tramite 2(1-s^2)/x^3
lasciando eps=1/(ell+1/2) identico, e le frequenze esatte restano
disponibili da Leaver.  Se l'indicatore e' utile deve ordinare gli errori
WKB dentro ciascun gruppo a ell fissato.

Uso: python3.13 robust_indicator_test.py
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "core"))

from leaver_qnm import leaver_qnm  # noqa: E402
from schwarzschild_wkb import (  # noqa: E402
    _derivative_functions,
    potential_peak,
    qnm_wkb,
    schwarzschild_potential,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))

from static_madelung_benchmark import (  # noqa: E402
    frobenius_log_derivative,
    inverse_tortoise,
    tortoise,
)


@dataclass(frozen=True)
class Indicators:
    integral: float
    median: float
    trimmed: float
    excised: float
    node_depth: float
    node_fraction: float


def profile(
    ell: int,
    overtone: int,
    spin: int,
    omega: complex,
    points: int = 4001,
    half_widths: float = 2.25,
) -> dict[str, np.ndarray | float]:
    """Profilo attraverso la barriera, con ampiezza esposta."""
    x_peak = potential_peak(ell, spin)
    derivatives = _derivative_functions(ell, spin)
    v0 = float(derivatives[0](x_peak))
    v2 = float(derivatives[2](x_peak))
    sigma = float(np.sqrt(v0 / (-v2)))
    center = float(tortoise(x_peak))
    xstar = np.linspace(center - half_widths * sigma, center + half_widths * sigma, points)
    x_left = float(inverse_tortoise(xstar[0]))

    def rhs(_t: float, state: np.ndarray) -> np.ndarray:
        psi, chi, radius = state
        radius_real = float(radius.real)
        f = 1.0 - 2.0 / radius_real
        return np.array(
            [chi, (schwarzschild_potential(radius_real, ell, spin) - omega**2) * psi, f],
            dtype=complex,
        )

    solution = solve_ivp(
        rhs,
        (float(xstar[0]), float(xstar[-1])),
        np.array(
            [1.0 + 0.0j, frobenius_log_derivative(x_left, omega, ell, spin), x_left],
            dtype=complex,
        ),
        t_eval=xstar,
        method="DOP853",
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    psi = solution.y[0]
    radius = solution.y[2].real
    epsilon = 1.0 / (ell + 0.5)
    amplitude = np.abs(psi)
    momentum = epsilon * np.imag(solution.y[1] / psi)
    potential = np.asarray(schwarzschild_potential(radius, ell, spin))

    step = float(xstar[1] - xstar[0])
    curvature = np.gradient(np.gradient(amplitude, step), step)
    q_madelung = -(epsilon**2) * curvature / amplitude
    density = abs(epsilon * omega) ** 2 + epsilon**2 * np.abs(potential) + momentum**2
    weight = np.exp(-0.5 * ((xstar - center) / sigma) ** 2)

    return {
        "xstar": xstar,
        "amplitude": amplitude,
        "q_madelung": q_madelung,
        "density": density,
        "weight": weight,
        "sigma": sigma,
        "epsilon": epsilon,
    }


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values_sorted = values[order]
    cumulative = np.cumsum(weights[order])
    return float(values_sorted[np.searchsorted(cumulative, 0.5 * cumulative[-1])])


def indicators(data: dict[str, np.ndarray | float]) -> Indicators:
    xstar = np.asarray(data["xstar"])
    amplitude = np.asarray(data["amplitude"])
    magnitude = np.abs(np.asarray(data["q_madelung"]))
    density = np.asarray(data["density"])
    weight = np.asarray(data["weight"])
    sigma = float(data["sigma"])
    interior = slice(3, -3)

    x_in = xstar[interior]
    mag = magnitude[interior]
    den = density[interior]
    wt = weight[interior]
    amp = amplitude[interior]

    integral = float(np.trapezoid(wt * mag, x_in) / np.trapezoid(wt * den, x_in))
    median = weighted_median(mag, wt) / weighted_median(den, wt)

    cutoff = np.quantile(mag, 0.90)
    keep = mag <= cutoff
    trimmed = float(
        np.trapezoid(wt[keep] * mag[keep], x_in[keep])
        / np.trapezoid(wt[keep] * den[keep], x_in[keep])
    )

    node = int(np.argmin(amp))
    outside = np.abs(x_in - x_in[node]) > 0.25 * sigma
    excised = float(
        np.trapezoid(wt[outside] * mag[outside], x_in[outside])
        / np.trapezoid(wt[outside] * den[outside], x_in[outside])
    )

    total = np.trapezoid(wt * mag, x_in)
    near = ~outside
    fraction = float(np.trapezoid(wt[near] * mag[near], x_in[near]) / total)
    return Indicators(
        integral=integral,
        median=median,
        trimmed=trimmed,
        excised=excised,
        node_depth=float(amp.min() / amp.max()),
        node_fraction=fraction,
    )


def degeneracy_test() -> None:
    """(a) La dipendenza da ell a n=1 collassa con le statistiche robuste?"""
    print("=== Test 1: la dipendenza da ell a n=1 e' artefatto del nodo? ===")
    print("Spread relativo di E/eps^2 su ell=3..8 (basso = degenere = nessuna informazione)")
    print()
    print("  s  n |  integrale   mediana   troncato   escisso | frazione nodo")
    for spin in (0, 2):
        for overtone in (0, 1):
            rows = []
            fractions = []
            for ell in range(max(spin, 3), 9):
                omega = leaver_qnm(ell, overtone, spin)
                values = indicators(profile(ell, overtone, spin, omega))
                epsilon = 1.0 / (ell + 0.5)
                rows.append(
                    (
                        values.integral / epsilon**2,
                        values.median / epsilon**2,
                        values.trimmed / epsilon**2,
                        values.excised / epsilon**2,
                    )
                )
                fractions.append(values.node_fraction)
            array = np.array(rows)
            spread = (array.max(axis=0) - array.min(axis=0)) / array.mean(axis=0)
            print(
                f"  {spin}  {overtone} |  {spread[0]:8.1%}  {spread[1]:8.1%}  "
                f"{spread[2]:8.1%}  {spread[3]:8.1%} | {np.mean(fractions):8.1%}"
            )
    print()


def fixed_epsilon_test() -> None:
    """(b) A eps fissato, l'indicatore ordina gli errori WKB fra spin diversi?"""
    print("=== Test 2: eps FISSATO, forma variata via spin (protocollo audit §5.5) ===")
    print("A ell fissato: eps identico, barriera diversa, omega esatte da Leaver.")
    print()
    print("  ell  s |   err.WKB3    |  E integrale    E mediana    E troncato")
    concordant = {"integral": 0, "median": 0, "trimmed": 0}
    total_pairs = 0
    for ell in range(3, 9):
        group = []
        for spin in (0, 1, 2):
            exact = leaver_qnm(ell, 0, spin)
            approximate = qnm_wkb(ell, 0, spin, 3).omega_M
            error = abs(exact - approximate) / abs(exact)
            values = indicators(profile(ell, 0, spin, exact))
            group.append((spin, error, values))
            print(
                f"  {ell:3d}  {spin} | {error:.6e} | {values.integral:.6e} "
                f"{values.median:.6e} {values.trimmed:.6e}"
            )
        for (_, error_a, a), (_, error_b, b) in combinations(group, 2):
            total_pairs += 1
            for key, va, vb in (
                ("integral", a.integral, b.integral),
                ("median", a.median, b.median),
                ("trimmed", a.trimmed, b.trimmed),
            ):
                if (error_a - error_b) * (va - vb) > 0:
                    concordant[key] += 1
        print()
    print(f"Coppie confrontate a eps fissato: {total_pairs}")
    for key, count in concordant.items():
        print(f"  concordanza {key:10}: {count}/{total_pairs} = {count / total_pairs:.1%}")
    print("(50% = indistinguibile dal caso; serve nettamente di piu' per essere utile)")


def main() -> None:
    degeneracy_test()
    fixed_epsilon_test()


if __name__ == "__main__":
    main()
