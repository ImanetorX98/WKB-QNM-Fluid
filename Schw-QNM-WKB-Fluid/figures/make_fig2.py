#!/usr/bin/env python3
"""Figura 2: l'ostruzione di condizionamento (§5).

Due pannelli sullo stesso oggetto, il potenziale quantistico previsto
analiticamente, valutati dove Q_M e' noto in forma chiusa (Poschl-Teller, n=0).

(a) A frequenza **esatta** la previsione e' ottima e converge: l'errore va come
    eps^2, e includere u_2 guadagna un ordine.  Non e' un metodo che non
    funziona.
(b) Perturbando la frequenza di una frazione relativa delta, l'errore su Q_M
    cresce come ~10^2 * delta.  La retta grigia e' l'assenza di amplificazione.
    I dati stanno due decadi sopra.

Il messaggio e' la giustapposizione: **il funzionale e' accurato e inservibile
per la stessa ragione.**  Dipende dalla frequenza cosi' fortemente che per
calcolarlo a tre cifre ne servono cinque di cio' che si vorrebbe prevedere.

Colori Okabe-Ito, identita' codificata anche da marcatore e etichetta diretta.

Uso: python3.13 make_fig2.py
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

from madelung_wkb_prediction import poschl_teller_validation, predict  # noqa: E402

BLUE, VERMILLION, GREEN = "#0072B2", "#D55E00", "#009E73"


def sensitivity(scale: float = 50.0, points: int = 40001,
                deltas: tuple[float, ...] = (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2)):
    """Errore su Q_M quando la frequenza e' sbagliata di una frazione relativa."""
    from poschl_teller_madelung_benchmark import exact_frequency, exact_wavefunction

    epsilon = 1.0 / scale
    y = np.linspace(-3.0, 3.0, points)
    step = float(y[1] - y[0])
    psi, _ = exact_wavefunction(y, scale, 0)
    amplitude = np.abs(psi)
    exact = -(epsilon**2) * np.gradient(np.gradient(amplitude, step), step) / amplitude
    omega = exact_frequency(scale, 0) / scale
    window = (y > 1.0) & (y < 2.5)

    rows = []
    for delta in deltas:
        q = (omega * (1.0 + delta)) ** 2 - 1.0 / np.cosh(y) ** 2
        predicted = predict(q, epsilon, step, 2)
        rows.append((delta, float(np.median(
            np.abs(predicted - exact)[window] / np.abs(exact[window])))))
    return rows


def main() -> None:
    figure, (left, right) = plt.subplots(1, 2, figsize=(7.6, 3.5), constrained_layout=True)

    # --- (a) convergenza a frequenza esatta -------------------------------
    rows = poschl_teller_validation()
    eps = np.array([1.0 / r["L"] for r in rows])
    for key, color, marker, label in (
        ("order0", BLUE, "o", r"$u=\sqrt{q}$"),
        ("order2", VERMILLION, "s", r"$u=\sqrt{q}+\varepsilon^2u_2$"),
    ):
        val = np.array([r[key] for r in rows])
        slope = np.polyfit(np.log(eps), np.log(val), 1)[0]
        left.plot(eps, val, marker=marker, color=color, linewidth=2.0, markersize=6.0,
                  markeredgecolor="white", markeredgewidth=0.8,
                  label=f"{label} — pendenza {slope:.2f}", zorder=3)
    left.set_xscale("log"); left.set_yscale("log")
    left.set_xlabel(r"$\varepsilon = 1/L$")
    left.set_ylabel(r"errore relativo su $Q_M$")
    left.set_title(r"(a)  frequenza esatta", fontsize=10, loc="left")
    left.legend(frameon=False, fontsize=8.5, loc="lower right", handlelength=1.6)

    # --- (b) condizionamento ----------------------------------------------
    data = sensitivity()
    floor = data[0][1]
    delta = np.array([d for d, _ in data[1:]])
    error = np.array([e for _, e in data[1:]])

    right.axhline(floor, color="0.75", linewidth=1.0, linestyle=(0, (1, 2)), zorder=1)
    right.annotate("pavimento a frequenza esatta", xy=(delta[0], floor),
                   xytext=(2, 5), textcoords="offset points",
                   color="0.45", fontsize=8.5)

    span = np.array([delta.min() * 0.6, delta.max() * 1.6])
    right.plot(span, span, linestyle=(0, (5, 4)), color="0.6", linewidth=1.0, zorder=1)
    right.annotate("nessuna amplificazione", xy=(span[1], span[1]),
                   xytext=(-4, 6), textcoords="offset points", ha="right",
                   color="0.45", fontsize=8.5)

    right.plot(delta, error, marker="^", color=GREEN, linewidth=2.0, markersize=6.5,
               markeredgecolor="white", markeredgewidth=0.8, zorder=3,
               label=r"errore su $Q_M$")
    strong = error > 3.0 * floor
    factor = float(np.median(error[strong] / delta[strong]))
    right.annotate(rf"$\times\,{factor:.0f}$", xy=(delta[-2], error[-2]),
                   xytext=(-30, 4), textcoords="offset points",
                   color=GREEN, fontsize=10, fontweight="bold")

    right.set_xscale("log"); right.set_yscale("log")
    right.set_xlabel(r"errore relativo su $\omega$")
    right.set_ylabel(r"errore relativo su $Q_M$")
    right.set_title(r"(b)  frequenza perturbata", fontsize=10, loc="left")

    for axes in (left, right):
        axes.grid(True, which="major", linewidth=0.5, color="0.92", zorder=0)
        axes.set_axisbelow(True)
        for side in ("top", "right"):
            axes.spines[side].set_visible(False)

    output = Path(__file__).resolve().parent / "fig2_condizionamento"
    figure.savefig(output.with_suffix(".pdf"))
    figure.savefig(output.with_suffix(".png"), dpi=200)
    plt.close(figure)

    print("(a) convergenza a frequenza esatta:")
    for r in rows:
        print(f"    L={r['L']:6.1f}   sqrt(q): {r['order0']:.3e}   +u2: {r['order2']:.3e}")
    print(f"\n(b) condizionamento, L=50, pavimento {floor:.3e}:")
    for d, e in data[1:]:
        print(f"    delta={d:.0e}   errore={e:.3e}   rapporto={e/d:7.1f}")
    print(f"\namplificazione mediana (regime dominato dal segnale): {factor:.0f}")
    print(f"scritto: {output}.pdf e .png")


if __name__ == "__main__":
    main()
