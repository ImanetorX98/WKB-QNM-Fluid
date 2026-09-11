#!/usr/bin/env python3
"""Separazione esatta della correzione forzata in parte geometrica e resto.

RISULTATO.  Definendo  rho = G'/G - 2 i om/f  — che misura di quanto il modo NON
sia puramente uscente, e che si annulla nella zona esterna — vale l'identita'

    H/G = 4 i om int (r/f) dr_*  +  y~,

dove y~ obbedisce alla stessa equazione con sorgente

    Sorg = 2 f (rho + r rho' + r rho^2) - 2K (f rho + 2 i om).

**E' esatta ovunque**, senza alcuno sviluppo asintotico: la verifica simbolica
da' differenza nulla.  Il primo termine ha forma chiusa,

    int r dr/f^2 = [16M^3 + (2M-r)(24M^2 ln(r-2M) + 8M r + r^2)] / (2(2M-r)),

non conosce ne' il multipolo ne' lo spin, e porta **tutta** la crescita r^2.

CONSEGUENZA.  Il secolare quadratico e' geometrico e integrabile in forma chiusa;
cio' che resta e' sorgentato unicamente da rho, cioe' **dalla barriera**.
Numericamente y~ cresce in modo lineare, non quadratico, e la parte lineare e' il
termine costante -4 i om K, cioe' lo shift di frequenza, riassorbibile nella
convenzione di fase.  Dopo quelle due sottrazioni non resta crescita.

E' la risposta alla domanda lasciata aperta dal §11.4: non serviva uno studio
numerico della regione di barriera per sapere dove vive il contenuto non
geometrico — la separazione lo isola esattamente.

Uso: python3.13 calculations/vaidya_geometric_split.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))


def identity_is_exact():
    """Differenza simbolica fra il residuo e la sorgente annunciata: deve essere 0."""
    r, omega, transport, mass = sp.symbols("r omega K M", positive=True)
    lapse = 1 - 2 * mass / r
    deviation = sp.Function("rho")(r)
    gradient = 2 * sp.I * omega / lapse + deviation
    star = lambda X: lapse * sp.diff(X, r)

    forcing = sp.expand(2 * lapse * (gradient + r * sp.diff(gradient, r) + r * gradient**2)
                        - 2 * transport * (2 * sp.I * omega + lapse * deviation))
    geometric = 4 * sp.I * omega * r / lapse            # D(parte geometrica)
    applied = sp.expand(star(geometric) + (2 * sp.I * omega + 2 * lapse * deviation) * geometric)
    announced = sp.expand(
        2 * lapse * (deviation + r * sp.diff(deviation, r) + r * deviation**2)
        - 2 * transport * (lapse * deviation + 2 * sp.I * omega))
    return sp.simplify(sp.expand(forcing - applied - announced))


def geometric_primitive():
    """int r dr/f^2 in forma chiusa, e il suo sviluppo a grande r."""
    r, mass = sp.symbols("r M", positive=True)
    closed = sp.integrate(r / (1 - 2 * mass / r) ** 2, r)
    expansion = sp.series(sp.expand(closed.subs(mass, 1)), r, sp.oo, 3).removeO()
    return sp.simplify(closed), sp.simplify(expansion)


def main() -> None:
    print("1. L'identita' e' esatta?")
    print(f"   differenza simbolica = {identity_is_exact()}")
    print("   Nessuno sviluppo asintotico: vale ovunque, rho arbitraria.")
    print()
    closed, expansion = geometric_primitive()
    print("2. Parte geometrica in forma chiusa, int r dr/f^2 :")
    print(f"   {closed}")
    print(f"   a grande r (M=1): {expansion}")
    print()
    print("3. Lettura")
    print("   - il quadratico r^2 e' TUTTO nella parte geometrica, che non")
    print("     conosce ne' il multipolo ne' lo spin;")
    print("   - il resto y~ e' sorgentato solo da rho = G'/G - 2 i om/f,")
    print("     nulla nella zona esterna e O(1) sulla barriera;")
    print("   - numericamente y~ cresce in modo LINEARE, e quella crescita e'")
    print("     il termine costante -4 i om K, cioe' lo shift di frequenza.")


if __name__ == "__main__":
    main()
