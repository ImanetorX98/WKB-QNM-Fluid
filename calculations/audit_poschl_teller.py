#!/usr/bin/env python3
"""Reproduce the spectral audit; no changes to the existing benchmark solvers.

The Jost/Newton baseline uses only L and a trial frequency. Exact eigenvalues
are used separately for validation. This is a standard spectral baseline, not
a new Madelung estimator. All results are printed; no files are overwritten.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from calculations.poschl_teller_madelung_benchmark import (
    madelung_indicator,
    wkb_frequency,
)


def exact_omega(scale: mp.mpf, n: int) -> mp.mpc:
    return mp.sqrt(scale**2 - mp.mpf(1)/4) - 1j * (mp.mpf(n) + mp.mpf(1)/2)


def trial_omega(scale: mp.mpf, n: int, order: int) -> mp.mpc:
    overtone = mp.mpf(n) + mp.mpf(1)/2
    if order == 1:
        return mp.sqrt(scale**2 - 2j * overtone * scale)
    if order == 3:
        return mp.sqrt(scale**2 - overtone**2 - mp.mpf(1)/4
                       - 2j * overtone * scale * (1 - 1/(8 * scale**2)))
    raise ValueError("order must be 1 or 3")


def parameters(scale: mp.mpf, omega: mp.mpc) -> tuple:
    spectral_lambda = mp.sqrt(scale**2 - mp.mpf(1)/4)
    a = mp.mpf(1)/2 + 1j * spectral_lambda - 1j * omega
    b = mp.mpf(1)/2 - 1j * spectral_lambda - 1j * omega
    return a, b, 1 - 1j * omega


def incoming_coefficient(scale: mp.mpf, omega: mp.mpc) -> mp.mpc:
    a, b, c = parameters(scale, omega)
    return mp.gamma(c) * mp.gamma(-1j * omega) * mp.rgamma(a) * mp.rgamma(b)


def newton_correction(scale: mp.mpf, omega: mp.mpc) -> mp.mpc:
    a, b, c = parameters(scale, omega)
    # Intended for off-root trials; digamma has poles exactly at a QNM.
    slope = -1j * (mp.digamma(c) + mp.digamma(-1j * omega)
                   - mp.digamma(a) - mp.digamma(b))
    return -1/slope


def check_identities() -> None:
    # Use exactly the fixed grid/weight of the original n=0 indicator.
    sigma = 1/np.sqrt(2.0)
    y = np.linspace(-2.25 * sigma, 2.25 * sigma, 4001)
    weight = np.exp(-0.5 * (y/sigma)**2)
    moment = np.trapezoid(weight * np.tanh(y)**2, y)/np.trapezoid(weight, y)
    max_relative = 0.0
    for scale in (1.5, 2, 4, 8, 16, 32, 64):
        numerical, _, nodal = madelung_indicator(scale, 0)
        eps2 = 1/scale**2
        closed = eps2 * (2 - moment)/(8 - eps2 * moment)
        relative = abs(numerical/closed - 1)
        assert not nodal and relative < 1e-10, (scale, relative)
        max_relative = max(max_relative, relative)
    print(f"Fixed-weight moment c = {moment:.15g}")
    print(f"Closed-form indicator: max relative discrepancy = {max_relative:.3e}")

    tolerance = mp.mpf(10)**(-mp.mp.dps//2)
    max_connection = mp.mpf(0)
    max_derivative = mp.mpf(0)
    for scale in map(mp.mpf, (4, 16)):
        for n in (0, 1, 2):
            for order in (1, 3):
                trial = trial_omega(scale, n, order)
                assert abs(complex(trial) - wkb_frequency(float(scale), n, order)) < 1e-12
                a, b, c = parameters(scale, trial)
                outgoing = mp.gamma(c)*mp.gamma(c-a-b)*mp.rgamma(c-a)*mp.rgamma(c-b)
                incoming = incoming_coefficient(scale, trial)
                # Independent Gauss connection check, including normalization.
                for z in map(mp.mpf, ("0.37", "0.8")):
                    direct = mp.hyp2f1(a, b, c, z)
                    connected = (outgoing * mp.hyp2f1(a, b, a+b-c+1, 1-z)
                                 + incoming * (1-z)**(c-a-b)
                                 * mp.hyp2f1(c-a, c-b, c-a-b+1, 1-z))
                    relative = abs(direct-connected)/max(1, abs(direct))
                    assert relative < tolerance, relative
                    max_connection = max(max_connection, relative)
                differentiated = mp.diff(lambda w: incoming_coefficient(scale, w), trial)
                correction = newton_correction(scale, trial)
                derivative_error = abs(correction/(-incoming/differentiated) - 1)
                assert derivative_error < tolerance, derivative_error
                max_derivative = max(max_derivative, derivative_error)
                root = exact_omega(scale, n)
                assert abs(incoming_coefficient(scale, root)) < tolerance
    print("Gauss connection: max scaled discrepancy =", mp.nstr(max_connection, 5))
    print("Jost derivative: max relative discrepancy =", mp.nstr(max_derivative, 5))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--uniform-tail", action="store_true",
                        help="also rerun the original uniform indicator for L=16..64")
    args = parser.parse_args()
    if args.dps < 40:
        parser.error("--dps must be at least 40")
    mp.mp.dps = args.dps
    print(f"mpmath {mp.__version__}, working precision {mp.mp.dps} digits")
    check_identities()
    print("\nn L order relative_error estimator/true_error error_after_one_Newton_step")
    for n in (0, 1, 2):
        for scale in map(mp.mpf, (4, 16, 64)):
            for order in (1, 3):
                trial = trial_omega(scale, n, order)
                exact = exact_omega(scale, n)
                correction = newton_correction(scale, trial)
                error = abs(trial-exact)/abs(exact)
                # Use the trial frequency, not the exact one, in the estimator.
                estimate = abs(correction)/abs(trial)
                corrected_error = abs(trial+correction-exact)/abs(exact)
                print(n, int(scale), order, *(mp.nstr(v, 12) for v in
                      (error, estimate/error, corrected_error)))
                # No assertion of unconditional improvement: n=2,L=4,WKB1 worsens.

    print("\nAsymptotic coefficients (L=256): L^2 E1 -> B/2; L^5 E3 -> N/128")
    scale = mp.mpf(256)
    for n in (0, 1, 2):
        overtone = mp.mpf(n) + mp.mpf(1)/2
        exact = exact_omega(scale, n)
        measured = [scale**power * abs(trial_omega(scale, n, order)/exact-1)
                    for order, power in ((1, 2), (3, 5))]
        expected = [(overtone**2 + mp.mpf(1)/4)/2, overtone/128]
        for actual, limit in zip(measured, expected):
            assert abs(actual/limit-1) < mp.mpf("0.001")
        print(n, "measured", *(mp.nstr(v, 12) for v in measured),
              "limits", *(mp.nstr(v, 12) for v in expected))
    if args.uniform_tail:
        from calculations.uniform_madelung_defect import benchmark_row, power_law_exponent
        print("\nUniform tail fits for L=16,24,32,48,64 (double precision ODE)")
        for n in (0, 2):
            rows = [benchmark_row(scale, n) for scale in (16, 24, 32, 48, 64)]
            print(n, {field: power_law_exponent(rows, field, minimum_scale=16)
                      for field in ("uniform_defect", "full_madelung",
                                    "wkb1_relative_error", "wkb3_relative_error")})
    print("All numerical audit assertions passed.")


if __name__ == "__main__":
    main()
