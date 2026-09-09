#!/usr/bin/env python3
"""RIPRODUTTORE MINIMO della discrepanza 1.65e-4 sul numeratore di Vaidya.

QUESTO E' IL PROGRAMMA CHE GENERA L'ERRORE.  Eseguirlo stampa la discrepanza e
i controlli che l'hanno gia' esclusa da quattro cause.

--------------------------------------------------------------------------
IL PROBLEMA IN UNA RIGA
--------------------------------------------------------------------------
L'integrale definito dell'integrando calcolato numericamente e la differenza
della sua antiderivata analitica dovrebbero coincidere.  Coincidono a
0.99983, cioe' con un errore RELATIVO COSTANTE di 1.65e-4, mentre ogni
ingrediente e' verificato a 1e-7 o meglio.

Perche' e' bloccante: un errore relativo costante su un integrale che cresce
esponenzialmente produce un residuo che cresce esponenzialmente.  E' questo a
impedire l'indipendenza da L+ del numeratore della condizione di solvibilita'
(Leung et al. 1998, eq. 2.15), e quindi a bloccare il calcolo della memoria di
Vaidya.

--------------------------------------------------------------------------
GLI OGGETTI
--------------------------------------------------------------------------
Modo QNM di Schwarzschild, s=0, n=0, M=1, omega da Leaver (esatta a 10 cifre).

  Z    = e^{i omega r_*} R          ampiezza in coordinate entranti
  S    = 2 (r dZ/dr)'               sorgente non adiabatica, = -2 d_r d_M Z
  I    = mu Z S,   mu = e^{-2 i omega r_*}     integrando del numeratore

Nella regione asintotica Z = C e^{2 i omega r_*} u(r) con u = sum a_k r^-k, e

  I = C^2 W(r) e^{2 i omega r_*},   W = 2u[(2i omega/f) P + P'],
                                    P = r(2i omega u/f + u')

L'antiderivata analitica e' F = C^2 V(r) e^{2 i omega r_*}, con V dalla
ricorsione V <- (f/2i omega)(W - V'), derivate ESATTE su serie troncate in 1/r.

--------------------------------------------------------------------------
CAUSE GIA' ESCLUSE — non rifarle
--------------------------------------------------------------------------
1. QUADRATURA.  Trapezio e Simpson danno lo stesso rapporto a 7 cifre.
2. METODO DI DERIVATA.  G'' per differenze finite e G'' ricavata dall'ODE
   danno lo stesso integrando (a meno del fattore 2 di bookkeeping fra
   `source_integrand` e la vecchia `projection`, che moltiplica per 0.5).
3. COSTANTE DI NORMALIZZAZIONE C.  Convergente a 12 cifre; il rapporto e'
   insensibile al troncamento di u (identico a 8 cifre per 8..20 termini).
4. ALGEBRA DELLE SERIE.  f*(1/f)=1 esatto; la ricorsione soddisfa
   V' + (2i omega/f)V = W con residuo 1.8e-15; e su griglia fine indipendente
   l'antiderivata riproduce int W e^{2i omega r_*} dr a 1.000000000.

--------------------------------------------------------------------------
SOSPETTO RESIDUO, NON VERIFICATO
--------------------------------------------------------------------------
La serie di V e' ASINTOTICA: il residuo 1.8e-15 misura l'algebra, non l'errore
di troncamento a r finito.  A r ~ 55 la coda trascurata potrebbe valere 1e-4.
Contro questa ipotesi: variando il troncamento fra 8 e 20 termini il rapporto
non si muove di una cifra.  Va capito se `Series.truncate` stia davvero
cambiando il numero di termini efficaci o se stia tagliando altrove.

Uso: python3.13 repro_numerator_discrepancy.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from asymptotic_series import Series, antiderivative_series, integrand_series  # noqa: E402
from leaver_qnm import leaver_qnm  # noqa: E402
from outgoing_asymptotics import series_coefficients  # noqa: E402
from vaidya_solvability import radial_mode  # noqa: E402

ELL, SPIN = 2, 0


def setup(points: int = 100001, x_max: float = 90.0):
    omega = leaver_qnm(ELL, 0, SPIN)
    data = radial_mode(ELL, SPIN, omega, 2.0001, x_max, points)
    r, xstar = data["x"], data["xstar"]
    lapse = 1.0 - 2.0 / r

    amplitude = np.exp(1j * omega * xstar) * data["R"]
    derivative = np.exp(1j * omega * xstar) * (1j * omega * data["R"] + data["dR"]) / lapse
    source = 2.0 * np.gradient(r * derivative, r)
    integrand = np.exp(-2j * omega * xstar) * amplitude * source

    series = Series(0, series_coefficients(ELL, SPIN, omega, 20))
    constant = (amplitude / (np.exp(2j * omega * xstar) * series.evaluate(r)))[
        np.searchsorted(r, 70.0)
    ]
    momentum = antiderivative_series(integrand_series(ELL, SPIN, omega, 20), omega, 14, 20)
    antiderivative = constant**2 * momentum.evaluate(r) * np.exp(2j * omega * xstar)
    return omega, r, xstar, integrand, antiderivative, constant


def main() -> None:
    omega, r, xstar, integrand, antiderivative, constant = setup()
    print(f"Schwarzschild s=0 n=0 ell={ELL},  M*omega = {omega.real:.9f}{omega.imag:+.9f}i")
    print(f"costante di normalizzazione C = {constant.real:+.9f}{constant.imag:+.9f}i")
    print()

    print("=== LA DISCREPANZA ===")
    print("integrale numerico / differenza dell'antiderivata, su intervalli diversi")
    print("   da    a     trapezio        Simpson")
    for low, high in ((40.0, 50.0), (50.0, 60.0), (60.0, 70.0), (70.0, 80.0)):
        i, j = np.searchsorted(r, low), np.searchsorted(r, high)
        analytic = antiderivative[j] - antiderivative[i]
        trap = np.trapezoid(integrand[i:j], r[i:j]) / analytic
        simp = simpson(integrand[i:j], x=r[i:j]) / analytic
        print(f"  {low:4.0f} {high:4.0f}   {abs(trap):.9f}   {abs(simp):.9f}")
    print()
    print("  Il rapporto e' COSTANTE a 0.99983: errore relativo 1.65e-4.")
    print("  Se fosse quadratura, trapezio e Simpson differirebbero.  Non lo fanno.")
    print()

    print("=== CONTROLLO: l'integrando coincide con la serie? ===")
    series = Series(0, series_coefficients(ELL, SPIN, omega, 20))
    predicted = (
        constant**2
        * integrand_series(ELL, SPIN, omega, 20).evaluate(r)
        * np.exp(2j * omega * xstar)
    )
    print("      r     rapporto COMPLESSO numerico/serie")
    for radius in (40.0, 55.0, 70.0, 85.0):
        i = np.searchsorted(r, radius)
        ratio = integrand[i] / predicted[i]
        print(f"   {radius:5.0f}    {ratio.real:.9f}{ratio.imag:+.9f}i")
    print()
    print("  Coincidono a 2e-7 COME NUMERI COMPLESSI, non solo in modulo.")
    print("  Quindi l'integrando e' giusto e l'antiderivata e' esatta per la serie,")
    print("  eppure i loro integrali definiti differiscono di 1.65e-4.")
    print()

    print("=== CONTROLLO: il troncamento sposta il risultato? ===")
    i, j = np.searchsorted(r, 50.0), np.searchsorted(r, 60.0)
    numeric = np.trapezoid(integrand[i:j], r[i:j])
    print("   termini   rapporto")
    for terms in (8, 12, 16, 20):
        v_local = antiderivative_series(
            integrand_series(ELL, SPIN, omega, terms), omega, 12, terms
        )
        f_local = constant**2 * v_local.evaluate(r) * np.exp(2j * omega * xstar)
        print(f"     {terms:3d}     {abs(numeric / (f_local[j] - f_local[i])):.9f}")
    print()
    print("  Insensibile: il troncamento della serie non spiega la discrepanza,")
    print("  oppure `Series.truncate` non sta cambiando i termini efficaci.")
    print("  QUESTO E' IL PUNTO DA CUI RIPARTIRE.")


if __name__ == "__main__":
    main()
