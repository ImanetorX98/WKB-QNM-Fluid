#!/usr/bin/env python3
"""Figura 3: la ridondanza del diagnostico di Madelung (§6).

Il funzionale d'ampiezza E_M contro |Lambda_3|, la correzione di Iyer-Will, per
sei multipoli e tre spin.  A overtone fissato i punti cadono su una retta di
pendenza 1 in log-log, cioe' E_M = c(n) |Lambda_3| con c indipendente da ell e
dallo spin.  Il pannello destro mostra il residuo attorno alla costante.

Non e' una correlazione: e' una **proporzionalita'**.  Lambda_3 costa sei
derivate del potenziale nel picco; E_M costa la frequenza esatta, dati iniziali
di Frobenius e l'integrazione dell'ODE attraverso la barriera.  Stessa
informazione, costo incomparabile.

Colori Okabe-Ito; lo spin e' codificato dal marcatore, non dal colore.

Uso: python3.13 make_fig3.py
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
from robust_indicator_test import indicators, profile  # noqa: E402
from schwarzschild_wkb import qnm_wkb  # noqa: E402

BLUE, VERMILLION = "#0072B2", "#D55E00"
MARKERS = {0: "o", 1: "s", 2: "^"}


def gather(overtone: int, ells=range(3, 9), spins=(0, 2)):
    rows = []
    for ell in ells:
        for spin in spins:
            omega = leaver_qnm(ell, overtone, spin)
            lam = abs(qnm_wkb(ell, overtone, spin, 3).lambda3)
            value = indicators(profile(ell, overtone, spin, omega)).median
            rows.append({"ell": ell, "spin": spin, "lambda3": lam, "E": value})
    return rows


def main() -> None:
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(7.6, 3.5), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]})

    summary = {}
    for overtone, color in ((0, BLUE), (1, VERMILLION)):
        rows = gather(overtone)
        lam = np.array([r["lambda3"] for r in rows])
        val = np.array([r["E"] for r in rows])
        ratio = val / lam
        constant = float(np.mean(ratio))
        spread = float(np.std(ratio) / constant)
        summary[overtone] = (constant, spread, rows, ratio)

        for spin in sorted({r["spin"] for r in rows}):
            sel = [i for i, r in enumerate(rows) if r["spin"] == spin]
            left.plot(lam[sel], val[sel], linestyle="none", marker=MARKERS[spin],
                      color=color, markersize=6.0, markeredgecolor="white",
                      markeredgewidth=0.8, zorder=3,
                      label=f"$n={overtone}$, $s={spin}$")
        grid = np.array([lam.min() * 0.7, lam.max() * 1.4])
        left.plot(grid, constant * grid, color=color, linewidth=1.2,
                  linestyle=(0, (5, 4)), zorder=2)
        left.annotate(rf"$c={constant:.3f}$", xy=(grid[1], constant * grid[1]),
                      xytext=(-6, -12), textcoords="offset points", ha="right",
                      color=color, fontsize=9)

        right.plot(np.arange(len(ratio)), ratio / constant - 1.0,
                   linestyle="none", marker="o", color=color, markersize=5.0,
                   markeredgecolor="white", markeredgewidth=0.7, zorder=3,
                   label=rf"$n={overtone}$ — dispersione {spread:.2%}")

    left.set_xscale("log"); left.set_yscale("log")
    left.set_xlabel(r"$|\Lambda_3|$  (Iyer–Will)")
    left.set_ylabel(r"$E_M$  (diagnostico di Madelung)")
    left.set_title("(a)  proporzionalità", fontsize=10, loc="left")
    left.legend(frameon=False, fontsize=8, loc="upper left", handlelength=1.0,
                ncol=2, columnspacing=1.0)

    right.axhline(0.0, color="0.6", linewidth=1.0, zorder=1)
    right.set_xlabel(r"modo  ($\ell=3\ldots8$, $s=0,2$)")
    right.set_ylabel(r"$E_M/(c\,|\Lambda_3|) - 1$")
    right.set_title("(b)  residuo attorno alla costante", fontsize=10, loc="left")
    right.legend(frameon=False, fontsize=8.5, loc="upper right", handlelength=1.0)

    for axes in (left, right):
        axes.grid(True, which="major", linewidth=0.5, color="0.92", zorder=0)
        axes.set_axisbelow(True)
        for side in ("top", "right"):
            axes.spines[side].set_visible(False)

    output = Path(__file__).resolve().parent / "fig3_ridondanza"
    figure.savefig(output.with_suffix(".pdf"))
    figure.savefig(output.with_suffix(".png"), dpi=200)
    plt.close(figure)

    for overtone, (constant, spread, rows, _) in summary.items():
        print(f"n={overtone}:  E_M/|Lambda_3| = {constant:.4f}  "
              f"dispersione {spread:.2%}  su {len(rows)} modi")
    print(f"scritto: {output}.pdf e .png")


if __name__ == "__main__":
    main()
