#!/usr/bin/env python3
"""Il funzionale di Madelung sul profilo radiale di Kerr.

Completa il §9 del manoscritto: la misura dell'ordinamento e' stata fatta sul
potenziale, qui si monta anche Q_M = -eps^2 A''/A sull'ampiezza integrata, in
modo da avere per Kerr la stessa tabella che si ha per Dirac.

FREQUENZA.  Non serve Leaver per Kerr.  Nel regime eikonale la condizione WKB1
di barriera, in forma riscalata, e'

    qhat(r0)/sqrt(2 qhat''(r0)) = i eps (n+1/2),        '' = d/dr_*

e a eps->0 si riduce alla radice doppia qhat(r0)=qhat'(r0)=0 gia' risolta in
`kerr_radial_order_profile.py`.  Perturbando Omega = Omega0 + delta attorno a
quella soluzione,

    delta = i eps (n+1/2) sqrt(2 qhat'') / (d qhat/d Omega),

per cui omega = L*Omega0 + i(n+1/2) sqrt(2 qhat'')/(d qhat/d Omega): la parte
immaginaria e' indipendente da L, come deve essere.  Il controllo a=0 deve dare
il valore noto omega = L/(3 sqrt 3) - i(n+1/2)/(3 sqrt 3).

CONDIZIONE AL CONTORNO.  All'orizzonte l'onda entrante e' co-rotante:
Psi ~ exp(-i k r_*) con k = omega - m Omega_H e Omega_H = a/(r_+^2+a^2).

Uso: python3.13 kerr_madelung_profile.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from kerr_eikonal_order_test import spheroidal_eigenvalue  # noqa: E402
from kerr_radial_order_profile import (  # noqa: E402
    eikonal_solution,
    horizon,
    leading_potential,
)


def tortoise_derivative(r: np.ndarray | float, spin: float) -> np.ndarray | float:
    """dr/dr_* = Delta/H."""
    return (r**2 - 2.0 * r + spin**2) / (r**2 + spin**2)


def eikonal_frequency(
    spin: float, mu: float, overtone: int, large_l: float
) -> tuple[complex, dict[str, float]]:
    """omega = L*Omega0 + i(n+1/2) sqrt(2 qhat'')/(d qhat/d Omega)."""
    solution = eikonal_solution(spin, mu)
    r0 = solution["r_peak"]
    omega0 = solution["omega_hat"]
    abar0 = solution["Abar0"]

    step = 1.0e-5

    def q_hat(r: float, omega: float) -> float:
        return float(leading_potential(r, spin, mu, omega, abar0))

    # d^2 qhat/dr_*^2 al picco, con d/dr_* = (Delta/H) d/dr
    def first_star(r: float) -> float:
        return tortoise_derivative(r, spin) * (
            q_hat(r + step, omega0) - q_hat(r - step, omega0)
        ) / (2.0 * step)

    second_star = tortoise_derivative(r0, spin) * (
        first_star(r0 + step) - first_star(r0 - step)
    ) / (2.0 * step)

    d_omega = (q_hat(r0, omega0 + step) - q_hat(r0, omega0 - step)) / (2.0 * step)
    damping = np.sqrt(2.0 * second_star + 0j) / d_omega
    omega = large_l * omega0 + 1j * (overtone + 0.5) * damping
    if omega.imag > 0:
        omega = omega.conjugate()
    return complex(omega), {**solution, "q_second_star": float(second_star)}


def kerr_profile(
    ell: int,
    spin: float,
    mu: float,
    overtone: int = 0,
    points: int = 24000,
    r_max: float = 60.0,
    horizon_offset: float = 1.0e-6,
    omega: complex | None = None,
) -> dict[str, np.ndarray | complex | float]:
    """Integra Psi in r_* uniforme, dall'orizzonte fino a r_max.

    La finestra di misura e' poi presa **lontano dai turning point**, come per
    Schwarzschild e Dirac: al massimo di barriera P->0 e la gerarchia si rompe
    (§7 del manoscritto), quindi misurare li' sarebbe l'errore che quel
    paragrafo documenta.

    `horizon_offset` va tenuto piccolo.  La condizione entrante e' esatta solo
    per Delta->0, e a ell alto Q_M e' cosi' piccolo che il ramo riflesso
    introdotto a bordo finito lo domina: la pendenza misurata passa da 0.57 a
    2.03 spostando il bordo da 1e-3 a 1e-6.  Sotto 1e-7 l'integratore cede,
    perche' r_* diverge logaritmicamente.
    """
    large_l = ell + 0.5
    m = int(round(mu * large_l))
    eikonal, solution = eikonal_frequency(spin, mu, overtone, large_l)
    if omega is None:
        # La frequenza eikonale ha errore ~1e-3 per a!=0, e Q_M amplifica di
        # ~1e2: per il confronto con la previsione va passata quella
        # autoconsistente di `kerr_wkb3_selfconsistent.py`.
        omega = eikonal
    epsilon = 1.0 / large_l
    r0 = solution["r_peak"]
    r_plus = horizon(spin)

    def rhs(_t: float, state: np.ndarray) -> np.ndarray:
        psi, chi, radius = state
        r = float(radius.real)
        h_squared = r**2 + spin**2
        delta = r**2 - 2.0 * r + spin**2
        eigenvalue = separation
        h = np.sqrt(h_squared)
        # (1/h) d^2 h/dr_*^2 valutato analiticamente
        dh = (delta / h_squared) * (r / h)
        d2h = (delta / h_squared) * (
            (2.0 * r - 2.0) / h_squared * (r / h)
            - 2.0 * r * delta / h_squared**2 * (r / h)
            + delta / h_squared * (1.0 / h - r**2 / h**3)
        )
        geometric = d2h / h
        q = (omega - m * spin / h_squared) ** 2 - delta * eigenvalue / h_squared**2 - geometric
        return np.array([chi, -q * psi, delta / h_squared], dtype=complex)

    # L'autovalore va valutato a c = a*omega COMPLESSA e tenuto complesso: per
    # a!=0 troncarne la parte immaginaria rende omega non piu' una frequenza di
    # QNM del potenziale usato, e la soluzione integrata acquista ammistione di
    # ramo entrante che una previsione a ramo puro non puo' riprodurre.
    eigen = complex(spheroidal_eigenvalue(ell, m, spin * omega))
    separation = eigen + spin**2 * omega**2 - 2.0 * spin * m * omega

    # Estensione in r_* stimata dalla mappa dr_*/dr = H/Delta, integrata
    # numericamente una volta sola per fissare la griglia.
    r_start = r_plus + horizon_offset * (r0 - r_plus)
    sample = np.linspace(r_start, r_max, 40000)
    rstar_span = float(np.trapezoid((sample**2 + spin**2) / (sample**2 - 2.0 * sample + spin**2), sample))
    horizon_frequency = spin / (r_plus**2 + spin**2)
    wavenumber = omega - m * horizon_frequency
    grid = np.linspace(0.0, rstar_span, points)

    psi0 = 1.0 + 0.0j
    solution_ivp = solve_ivp(
        rhs,
        (grid[0], grid[-1]),
        np.array([psi0, -1j * wavenumber * psi0, r_start], dtype=complex),
        t_eval=grid,
        method="DOP853",
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    if not solution_ivp.success:
        raise RuntimeError(solution_ivp.message)

    psi = solution_ivp.y[0]
    r = solution_ivp.y[2].real
    amplitude = np.abs(psi)
    step = float(grid[1] - grid[0])
    curvature = np.gradient(np.gradient(amplitude, step), step)
    q_madelung = -(epsilon**2) * curvature / amplitude

    h_squared = r**2 + spin**2
    delta = r**2 - 2.0 * r + spin**2
    rotation = (np.real(separation) / large_l**2 - solution["Abar0"]) * delta / h_squared**2
    # (per la parte di rotazione conta la parte reale: e' quella che entra in Re qhat)

    return {
        "rstar": grid,
        "r": r,
        "amplitude": amplitude,
        "q_madelung": q_madelung,
        "rotation": rotation,
        "omega": omega,
        "epsilon": epsilon,
        "r_peak": r0,
        "rstar_span": rstar_span,
        "r_max": r_max,
        "omega_eikonal": eikonal,
    }


def scaling(
    spin: float, mu: float, ells: tuple[int, ...] = (20, 30, 45, 70, 100)
) -> list[dict[str, float]]:
    rows = []
    for ell in ells:
        data = kerr_profile(ell, spin, mu)
        r = np.asarray(data["r"])
        # Stessa convenzione di Schwarzschild e Dirac: finestra lontana dai
        # turning point, non al massimo di barriera.
        mask = (r > 20.0) & (r < 50.0)
        interior = np.zeros_like(mask)
        interior[6:-6] = True
        mask = mask & interior
        rows.append(
            {
                "ell": float(ell),
                "epsilon": float(data["epsilon"]),
                "rotation": float(np.max(np.abs(np.asarray(data["rotation"])[mask]))),
                "madelung": float(np.max(np.abs(np.asarray(data["q_madelung"])[mask]))),
            }
        )
        rows[-1]["ratio"] = rows[-1]["rotation"] / rows[-1]["madelung"]
    return rows


def slopes(rows: list[dict[str, float]]) -> dict[str, float]:
    log_eps = np.log([row["epsilon"] for row in rows])
    return {
        key: float(np.polyfit(log_eps, np.log([row[key] for row in rows]), 1)[0])
        for key in ("rotation", "madelung")
    }


def main() -> None:
    print("=== Controllo: la frequenza eikonale riproduce Schwarzschild? ===")
    omega, _ = eikonal_frequency(0.0, 0.5, 0, 100.5)
    expected = 100.5 / (3 * np.sqrt(3)) - 0.5j / (3 * np.sqrt(3))
    print(f"  calcolata {omega.real:.8f}{omega.imag:+.8f}i")
    print(f"  attesa    {expected.real:.8f}{expected.imag:+.8f}i")
    print(f"  scarto    {abs(omega - expected):.2e}")
    print()

    for spin, mu in ((0.0, 0.5), (0.6, 0.5), (0.9, 0.5), (0.9, 0.9)):
        rows = scaling(spin, mu)
        exponents = slopes(rows)
        print(f"=== a={spin}, mu={mu} ===")
        print("   ell    eps      max|rotazione|   max|Q_M|     rapporto")
        for row in rows:
            print(
                f"  {row['ell']:4.0f}  {row['epsilon']:7.5f}   {row['rotation']:.6e}   "
                f"{row['madelung']:.6e}  {row['ratio']:8.3f}"
            )
        print(
            f"  pendenze in eps: rotazione = {exponents['rotation']:.4f}, "
            f"Q_M = {exponents['madelung']:.4f}"
        )
        print()


if __name__ == "__main__":
    main()
