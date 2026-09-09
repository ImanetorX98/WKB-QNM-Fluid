#!/usr/bin/env python3
"""Figura 1: il criterio di ordinamento, in un pannello.

Log-log del primo termine subprincipale contro eps, per i tre casi.  Ogni serie
e' riscalata in modo che la sua retta di fit passi per 1 in eps=0.1: le ampiezze
assolute hanno significati fisici diversi e non sono confrontabili, mentre le
pendenze lo sono, e con questa scelta la pendenza e' l'unica differenza visibile.
Le due rette grigie hanno pendenza 1 e 2.

Colori: Okabe-Ito (blu / vermiglio / verde-bluastro), validati per CVD
(peggior coppia adiacente Delta-E 11.0 in deuteranopia).  Identita' codificata
anche dal marcatore e da etichette dirette: mai dal solo colore.

Uso: python3.13 make_fig1.py
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
sys.path.insert(0, str(ROOT / "Schw-QNM-WKB-Fluid" / "verification"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, VERMILLION, GREEN = "#0072B2", "#D55E00", "#009E73"


def gather() -> list[dict]:
    from dirac_madelung_profile import hierarchy_scaling
    from kerr_radial_order_profile import order_measurement
    from scalar_eikonal_scaling import scaling_table

    scalar = scaling_table(ell_values=(4, 8, 16, 32), spin=2)
    kerr = order_measurement(0.6, 0.5, ells=(40, 60, 80, 120, 160))
    dirac = hierarchy_scaling(kappa_values=(2.0, 4.0, 8.0, 16.0), region="far")

    return [
        {
            "label": "Schwarzschild, $s=2$",
            "detail": r"$\varepsilon^2 v_2$",
            "eps": [row["epsilon"] for row in scalar],
            "val": [row["subleading_max"] for row in scalar],
            "color": BLUE,
            "marker": "o",
        },
        {
            "label": "Kerr, $a=0.6$",
            "detail": r"$-\varepsilon\,\Delta A_1/H^2$",
            "eps": list(kerr["epsilons"]),
            "val": list(kerr["deviations"]),
            "color": VERMILLION,
            "marker": "s",
        },
        {
            "label": "Dirac, $\\tau=+1$",
            "detail": r"$\tau\varepsilon h'$",
            "eps": [row["epsilon"] for row in dirac],
            "val": [row["spin_max"] for row in dirac],
            "color": GREEN,
            "marker": "^",
        },
    ]


def main() -> None:
    series = gather()
    figure, axes = plt.subplots(figsize=(5.4, 4.2), constrained_layout=True)
    anchor = 0.1  # tutte le serie riportate a 1 in eps=0.1: resta solo la pendenza

    for entry in series:
        eps = np.array(entry["eps"], dtype=float)
        val = np.array(entry["val"], dtype=float)
        order = np.argsort(eps)
        eps, val = eps[order], val[order]
        slope, intercept = np.polyfit(np.log(eps), np.log(val), 1)
        scale = np.exp(slope * np.log(anchor) + intercept)
        axes.plot(
            eps,
            val / scale,
            marker=entry["marker"],
            color=entry["color"],
            linewidth=2.0,
            markersize=6.5,
            markeredgecolor="white",
            markeredgewidth=0.8,
            label=f"{entry['label']}, {entry['detail']} — {slope:.2f}",
            zorder=3,
        )
        entry["slope"] = float(slope)

    span = np.array([5.0e-3, 6.0e-1])
    for exponent in (1, 2):
        axes.plot(span, (span / anchor) ** exponent, linestyle=(0, (5, 4)),
                  color="0.6", linewidth=1.0, zorder=1)
        axes.annotate(
            rf"$\varepsilon^{{{exponent}}}$",
            xy=(span[0], (span[0] / anchor) ** exponent),
            xytext=(4, -2 if exponent == 1 else 4),
            textcoords="offset points",
            color="0.45", fontsize=9.5,
        )

    axes.axvline(anchor, color="0.85", linewidth=0.8, zorder=0)
    axes.set_xscale("log"); axes.set_yscale("log")
    axes.set_xlim(4.5e-3, 7.0e-1)
    axes.set_xlabel(r"$\varepsilon = 1/L$")
    axes.set_ylabel("primo termine subprincipale (riscalato)")
    axes.grid(True, which="major", linewidth=0.5, color="0.92", zorder=0)
    axes.set_axisbelow(True)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    axes.legend(frameon=False, fontsize=8.5, loc="lower right",
                handlelength=1.6, borderaxespad=0.6)

    output = Path(__file__).resolve().parent / "fig1_ordinamento"
    figure.savefig(output.with_suffix(".pdf"))
    figure.savefig(output.with_suffix(".png"), dpi=200)
    plt.close(figure)

    print("pendenze misurate:")
    for entry in series:
        print(f"  {entry['label']:24} {entry['detail']:26} {entry['slope']:.4f}")
    print(f"\nscritto: {output}.pdf e .png")


if __name__ == "__main__":
    main()
