#!/usr/bin/env python3
"""Plot the standard Schwarzschild Dirac partners and the 2026 v1 comparison."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tmp" / "pdfs"
OUT.mkdir(parents=True, exist_ok=True)
PLOT_PATH = OUT / "dirac_partner_potentials.png"


def main() -> None:
    r = np.linspace(2.001, 12.0, 1800)
    f = 1.0 - 2.0 / r
    root_f = np.sqrt(f)
    v_plus = f / r**2 + root_f * (3.0 - r) / r**3
    v_minus = f / r**2 - root_f * (3.0 - r) / r**3
    v_preprint = f * (1.0 / r**2 - 1.0 / (r**3 * root_f))

    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=180)
    ax.plot(r, v_plus, color="#246B8E", lw=2.2, label=r"$V_+$ standard")
    ax.plot(r, v_minus, color="#2A8C82", lw=1.9, label=r"$V_-$ standard")
    ax.plot(
        r,
        v_preprint,
        color="#C9982D",
        lw=1.7,
        ls="--",
        label="potenziale arXiv:2605.28887v1, Eq. (7.36)",
    )
    ax.axvline(3.0, color="#7D8790", lw=0.9, ls=":", label=r"$r=3M$")
    ax.set_xlim(2.0, 9.0)
    ax.set_ylim(-0.002, 0.052)
    ax.set_xlabel(r"$r/M$")
    ax.set_ylabel(r"$M^2 V(r)$")
    ax.set_title(r"Dirac massless, $M=1$, $|\kappa|=1$")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, fontsize=7.7, loc="upper right")
    fig.tight_layout()
    fig.savefig(PLOT_PATH, bbox_inches="tight")
    plt.close(fig)
    print(PLOT_PATH)


if __name__ == "__main__":
    main()
