#!/usr/bin/env python3
"""Tabella WKB1/WKB3 per piu' multipoli e overtones."""

from __future__ import annotations

import argparse

from schwarzschild_wkb import qnm_wkb


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spin", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--ell-min", type=int, default=None)
    parser.add_argument("--ell-max", type=int, default=6)
    parser.add_argument("--n-max", type=int, default=2)
    args = parser.parse_args()

    ell_min = args.spin if args.ell_min is None else args.ell_min
    print(" s ell  n |             WKB1 (M omega) |             WKB3 (M omega) | |Delta|")
    print("-" * 91)
    for ell in range(ell_min, args.ell_max + 1):
        for n in range(args.n_max + 1):
            w1 = qnm_wkb(ell, n, args.spin, 1).omega_M
            w3 = qnm_wkb(ell, n, args.spin, 3).omega_M
            delta = abs(w3 - w1)
            print(
                f" {args.spin:d} {ell:3d} {n:2d} | "
                f"{w1.real: .9f}{w1.imag:+.9f}i | "
                f"{w3.real: .9f}{w3.imag:+.9f}i | {delta:.3e}"
            )


if __name__ == "__main__":
    main()
