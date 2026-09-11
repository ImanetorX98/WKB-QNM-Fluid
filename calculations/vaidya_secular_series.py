#!/usr/bin/env python3
"""Struttura asintotica esatta della correzione forzata H/G, senza adattamenti.

MOTIVO.  La prima analisi della crescita secolare procedeva per adattamento di
H/G a un polinomio in r_*.  Un adattamento su un intervallo finito non distingue
r^2 da r_*^2 e assorbe un eventuale logaritmo nei coefficienti polinomiali: ne
usciva un coefficiente quadratico basso dello 0.5% e la conclusione, sbagliata,
che H/G fosse un polinomio.  Qui la stessa quantita' e' ricavata **in forma
chiusa** dall'equazione asintotica.

--------------------------------------------------------------------------
RIDUZIONE
--------------------------------------------------------------------------
Con D = d/dr_* si ha  f L X = [D^2 - 2 i om D - f U] X, e ponendo H = G y
la parte omogenea cade:

    D^2 y + (2 i om + 2 f u'/u) D y = f (S - 2K G_r)/G .

E' del **primo ordine in p = D y**, quindi si inverte per serie senza integrare
nulla numericamente.  Con g1 = G'/G = 2 i om/f + u'/u,

    Phi = 2 f (g1 + r g1' + r g1^2) - 2K (2 i om + f u'/u),
    p   = [Phi - f p'] / (2 i om + 2 f u'/u),   iterata,

e infine y = int p dr_* = int (p/f) dr, termine a termine.

--------------------------------------------------------------------------
RISULTATO
--------------------------------------------------------------------------
I coefficienti dell'ODE sono serie pure in 1/r, quindi p/f lo e'; integrando,
la forma asintotica e' **forzata**:

    H/G = 2 i om r^2 + c_1 r + c_log ln r + sum_k c_k r^-k .

Due letture, entrambe correzioni a quanto scritto prima:

* Il quadratico e' in **r**, non in r_*.  Un termine r_*^2 produrrebbe r ln r,
  che qui non compare: i coefficienti non contengono logaritmi, quindi
  l'integrazione ne genera al piu' uno.
* **c_log non e' nullo.**  H/G non e' un polinomio, e l'affermazione contraria
  era un artefatto dell'adattamento.

Il coefficiente quadratico e' invece confermato, e stavolta senza adattarlo:
vale **esattamente 2 i om** (scarto zero macchina), cioe' -2 i om'(M) per om
proporzionale a 1/M.  E' il coefficiente di deriva di massa, e resta la firma
della retrodatazione.

--------------------------------------------------------------------------
FORMA CHIUSA, E UNIVERSALITA'
--------------------------------------------------------------------------
I rapporti c_1/(2 i om) e c_log/(2 i om) valgono **8 e 24 esatti** per ogni modo
provato — ell = 2,3,4 e s = 0,1,2 — mentre da 1/r in poi dipendono dal modo.  I
tre termini di testa sono quindi geometrici, non del modo, e si riconoscono:

    H/G = 4 i om  int (r/f) dr_*  +  O(1/r),

    int r dr/f^2 = [16M^3 + (2M-r)(24M^2 ln(r-2M) + 8M r + r^2)] / (2(2M-r))
                 = r^2/2 + 4M r + 12 M^2 ln r + O(1/r)   a grande r.

Moltiplicando per 4 i om si ottengono esattamente 2 i om (r^2 + 8r + 24 ln r) a
M=1.  Il secolare di testa e' dunque il trasporto di fase lungo il raggio
uscente pesato da (dr_*/dr)^2, e non conosce ne' il multipolo ne' lo spin.

Uso: python3.13 calculations/vaidya_secular_series.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "core"))

from asymptotic_series import Series, inverse_lapse, lapse  # noqa: E402
from leaver_qnm import leaver_qnm  # noqa: E402
from outgoing_asymptotics import series_coefficients  # noqa: E402

TERMS = 26


def _unit(terms: int) -> Series:
    coeffs = np.zeros(terms, dtype=complex)
    coeffs[0] = 1.0
    return Series(0, coeffs)


def _divide(numer: Series, denom: Series, terms: int) -> Series:
    """Divisione di serie, con termine di testa del denominatore non nullo."""
    top, bottom = numer.truncate(terms).coeffs, denom.truncate(terms).coeffs
    out = np.zeros(terms, dtype=complex)
    for n in range(terms):
        out[n] = (top[n] - sum(out[m] * bottom[n - m] for m in range(n))) / bottom[0]
    return Series(numer.offset - denom.offset, out)


def secular_structure(ell: int = 2, spin: int = 0, transport: complex = 0.0,
                      terms: int = TERMS, iterations: int = 8) -> dict:
    """Coefficienti asintotici esatti di H/G, senza alcun adattamento."""
    omega = leaver_qnm(ell, 0, spin)
    wave = 2j * omega

    profile = Series(0, series_coefficients(ell, spin, omega, terms))
    logarithmic = _divide(profile.derivative(), profile, terms)      # u'/u
    slope = lapse(terms)
    inverse = inverse_lapse(terms)

    def product(*items: Series) -> Series:
        out = items[0]
        for item in items[1:]:
            out = (out * item).truncate(terms)
        return out

    def times_radius(series: Series) -> Series:
        return Series(series.offset - 1, series.coeffs)

    gradient = (wave * inverse + logarithmic).truncate(terms)        # G'/G
    forcing = (2.0 * product(slope, (gradient + times_radius(gradient.derivative())
                                     + times_radius(product(gradient, gradient))
                                     ).truncate(terms))
               + (-2.0 * transport) * (wave * _unit(terms)
                                       + product(slope, logarithmic)).truncate(terms)
               ).truncate(terms)
    damping = (wave * _unit(terms) + 2.0 * product(slope, logarithmic)).truncate(terms)

    momentum = _divide(forcing, damping, terms)
    for _ in range(iterations):
        momentum = _divide(
            (forcing + (-1.0) * product(slope, momentum.derivative())).truncate(terms),
            damping, terms)

    integrand = product(momentum, inverse)      # p/f, da integrare in dr
    coefficients = {}
    for index in range(terms):
        power = integrand.offset + index        # termine c r^-power
        value = complex(integrand.coeffs[index])
        if power == 1:
            coefficients["log"] = value
        else:
            coefficients[1 - power] = value / (1 - power)
    return {"omega": omega, "quadratic": coefficients[2], "linear": coefficients[1],
            "log": coefficients["log"], "coefficients": coefficients,
            "ratios": {"linear": coefficients[1] / coefficients[2],
                       "log": coefficients["log"] / coefficients[2]}}


def main() -> None:
    print("Struttura asintotica di H/G, ricavata in forma chiusa")
    print("(nessun adattamento: i coefficienti escono dall'inversione della serie)")
    print()
    for ell, spin in ((2, 0), (3, 0), (2, 2)):
        row = secular_structure(ell, spin)
        omega = row["omega"]
        print(f"  ell={ell}, s={spin}:  M*omega = {omega.real:.6f}{omega.imag:+.6f}i")
        print(f"     r^2   {row['quadratic'].real:+.9f}{row['quadratic'].imag:+.9f}i")
        print(f"     2i*om {(2j*omega).real:+.9f}{(2j*omega).imag:+.9f}i"
              f"     scarto {abs(row['quadratic'] - 2j*omega):.2e}")
        print(f"     r^1   {row['linear'].real:+.9f}{row['linear'].imag:+.9f}i")
        print(f"     ln r  {row['log'].real:+.9f}{row['log'].imag:+.9f}i"
              f"   <- non nullo: H/G non e' un polinomio")
        print(f"     rapporti a 2i*om:  r^1 {row['ratios']['linear'].real:.6f}"
              f"   ln r {row['ratios']['log'].real:.6f}   <- 8 e 24, universali")
        print()


if __name__ == "__main__":
    main()
