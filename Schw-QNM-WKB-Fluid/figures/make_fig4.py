#!/usr/bin/env python3
"""Figura 4: la condizione di solvibilita' di Vaidya, ai due bordi (§11).

(a) Il numeratore regolarizzato N(L_+) contro il taglio esterno, normalizzato al
    proprio valore: cinque modi fondamentali, cinque rette piatte a 1e-7.  E' la
    verifica che la regolarizzazione funziona — un taglio arbitrario che non
    lascia traccia non e' piu' arbitrario.  Nessuna curva passa per zero: la
    proiezione risonante non si annulla.

(b) L'esponente al bordo interno, Re(-4 i omega M) = 4 M Im(omega), contro
    l'overtone.  La soglia -1 separa singolarita' integrabile da non integrabile.
    **Tutti i fondamentali stanno sopra**: per n=0 il bordo interno non chiede
    regolarizzazione, e il limite esiste gia'.  Da n=1 in poi si scende sotto, e
    la forma chiusa del §11.3 diventa continuazione analitica.

Il pannello (b) e' il motivo per cui la figura esiste: e' l'osservazione che non
ci aspettavamo, ed e' leggibile in un colpo d'occhio.

Uso: python3.13 make_fig4.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
sys.path.insert(0, str(ROOT / "calculations"))
sys.path.insert(0, str(ROOT / "core"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from leaver_qnm import leaver_qnm  # noqa: E402
from vaidya_numerator_factored import numerator  # noqa: E402

BLUE, VERMILLION, GREEN, ORANGE, PURPLE = (
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7")
MODES = ((2, 0, BLUE, "o"), (3, 0, VERMILLION, "s"), (2, 2, GREEN, "^"),
         (3, 2, ORANGE, "D"), (2, 1, PURPLE, "v"))
CUTS = (40.0, 50.0, 60.0, 70.0, 80.0)


def main() -> None:
    figure, (left, right) = plt.subplots(1, 2, figsize=(7.6, 3.5), constrained_layout=True)

    for ell, spin, color, marker in MODES:
        row = numerator(ell=ell, spin=spin, uppers=CUTS)
        cuts = np.array(CUTS)
        values = np.array([row["values"][c] for c in CUTS])
        reference = values.mean()          # media, non un punto: nessun taglio privilegiato
        left.plot(cuts, np.abs(values / reference - 1.0), marker=marker, color=color,
                  linewidth=1.6, markersize=5.5, markeredgecolor="white",
                  markeredgewidth=0.7, zorder=3,
                  label=rf"$\ell={ell}$, $s={spin}$:  $|N|={abs(row['N']):.1f}$")

    left.set_yscale("log")
    left.set_ylim(1e-10, 3e-5)
    left.axhspan(1e-10, 1e-6, color="0.94", zorder=0)
    left.annotate("scarto sotto $10^{-6}$", xy=(41, 3e-7), fontsize=8.5, color="0.35")
    left.set_xlabel(r"taglio esterno $L_+/M$")
    left.set_ylabel(r"$\left|\,N(L_+)/\bar N - 1\right|$")
    left.set_title("(a)  il taglio esterno non lascia traccia", fontsize=10, loc="left")
    left.legend(frameon=False, fontsize=7.0, loc="upper left", ncol=2,
                handlelength=1.2, columnspacing=0.9)

    overtones = np.arange(0, 4)
    for ell, spin, color, marker in MODES:
        exponent = [4.0 * leaver_qnm(ell, n, spin).imag for n in overtones]
        right.plot(overtones, exponent, marker=marker, color=color, linewidth=1.6,
                   markersize=6.0, markeredgecolor="white", markeredgewidth=0.7,
                   zorder=3, label=rf"$\ell={ell}$, $s={spin}$")

    right.axhline(-1.0, color="0.35", linewidth=1.4, linestyle=(0, (5, 3)), zorder=2)
    right.annotate("soglia di integrabilità", xy=(3.05, -1.0), xytext=(-4, 5),
                   textcoords="offset points", ha="right", color="0.3", fontsize=8.5)
    right.axhspan(-1.0, 0.4, color="0.93", zorder=0)
    right.annotate("integrabile:\nnessuna regolarizzazione", xy=(0.05, -0.45),
                   fontsize=8.0, color="0.35", va="center")
    right.set_xticks(overtones)
    right.set_xlabel(r"overtone $n$")
    right.set_ylabel(r"$4M\,\mathrm{Im}\,\omega$")
    right.set_ylim(-4.2, 0.4)
    right.set_title("(b)  l'esponente al bordo interno", fontsize=10, loc="left")
    right.legend(frameon=False, fontsize=7.5, loc="lower left", handlelength=1.2)

    for axes in (left, right):
        axes.grid(True, which="major", linewidth=0.5, color="0.92", zorder=0)
        axes.set_axisbelow(True)
        for side in ("top", "right"):
            axes.spines[side].set_visible(False)

    output = Path(__file__).resolve().parent / "fig4_vaidya"
    figure.savefig(output.with_suffix(".pdf"))
    figure.savefig(output.with_suffix(".png"), dpi=200)
    plt.close(figure)

    print("esponente al bordo interno, 4 M Im(omega):")
    for ell, spin, _, _ in MODES:
        vals = "  ".join(f"{4.0*leaver_qnm(ell, n, spin).imag:+.3f}" for n in overtones)
        print(f"  ell={ell} s={spin}:  {vals}")
    print(f"\nscritto: {output}.pdf e .png")


if __name__ == "__main__":
    main()
