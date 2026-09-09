#!/usr/bin/env python3
"""Esiste un termine di ordine eps^1 nel settore bosonico di Kerr?

Nel settore bosonico statico non esiste: con la sostituzione di Langer
ell(ell+1)=L^2-1/4, e il -1/4 e' puro eps^2.  In Kerr l'autovalore sferoidale
A_{lm}(c), c=a*omega, sostituisce ell(ell+1).  Se la sua espansione eikonale

    A/L^2 = A0 + A1/L + A2/L^2 + ...     (mu=m/L, chat=c/L fissati)

ha A1 != 0, allora la rotazione produce un termine di ordine eps^1 nel
potenziale efficace, esattamente come la connessione di spin lo produce per
Dirac.  In quel caso la pendenza 1 del criterio non e' una firma fermionica.

Robustezza rispetto alle convenzioni: le forme alternative dell'equazione
angolare differiscono per termini a^2 omega^2 = L^2 chat^2 e 2 a m omega =
2 L^2 mu chat, entrambi puramente di ordine L^2.  Spostano A0, non A1.

Non rimovibilita': ridefinendo L -> L+delta si ha A1 -> A1 - 2 delta A0.
Una costante unica potrebbe assorbire A1 solo se A1/(2A0) fosse la stessa per
tutti i modi.  Lo script lo verifica.

Il caso QNM ha c complessa.  Per un overtone fissato
omega = L*Omega - i(n+1/2)lambda, quindi Im(chat) = -eps*a(n+1/2)lambda e'
essa stessa di ordine eps: il test a chat complessa costante isola A1, e la
scalatura fisica aggiunge un ulteriore contributo eps^1, non lo cancella.

Uso: python3.13 kerr_eikonal_order_test.py
"""

from __future__ import annotations

import numpy as np


def cosine_matrix(m: int, ell_max: int) -> tuple[np.ndarray, np.ndarray]:
    """Elementi <l'|cos theta|l> nella base delle armoniche sferiche."""
    ells = np.arange(abs(m), ell_max + 1)
    size = ells.size
    matrix = np.zeros((size, size))
    for index, ell in enumerate(ells):
        if index + 1 < size:
            element = np.sqrt(
                (ell + 1 - m) * (ell + 1 + m) / ((2 * ell + 1) * (2 * ell + 3))
            )
            matrix[index, index + 1] = element
            matrix[index + 1, index] = element
    return ells, matrix


def spheroidal_eigenvalue(
    ell: int, m: int, c: complex, pad: int = 60, steps: int = 24
) -> complex:
    """A_{lm}(c) per c complessa, con tracciamento robusto del modo.

    Il tracciamento per prossimita' di autovalore fallisce quando due
    autovalori si avvicinano: la continuazione salta modo, e la successione in
    ell diventa erratica.  Si procede quindi in due fasi.

    1. Parte reale.  Per c reale la matrice e' reale simmetrica, gli autovalori
       non si incrociano e l'ordinamento crescente e' conservato per
       continuita': il modo e' l'indice ell-|m| della lista ordinata.
    2. Parte immaginaria.  Si continua da Re(c) a c seguendo l'**autovettore**,
       non l'autovalore: la sovrapposizione e' insensibile agli avvicinamenti.
    """
    ells, cosine = cosine_matrix(m, ell + pad)
    cosine_squared = cosine @ cosine
    diagonal = np.diag(ells * (ells + 1.0))
    index = ell - abs(m)

    real_part = float(np.real(c))
    values, vectors = np.linalg.eigh(diagonal - (real_part**2) * cosine_squared)
    current_value = complex(values[index])
    current_vector = vectors[:, index].astype(complex)

    imaginary_part = float(np.imag(c))
    if imaginary_part == 0.0:
        return current_value

    complex_diagonal = diagonal.astype(complex)
    for step in range(1, steps + 1):
        c_step = real_part + 1j * imaginary_part * step / steps
        values, vectors = np.linalg.eig(complex_diagonal - (c_step**2) * cosine_squared)
        overlap = np.abs(current_vector.conj() @ vectors)
        best = int(np.argmax(overlap))
        current_value = complex(values[best])
        current_vector = vectors[:, best]
    return current_value


def eikonal_fit(
    mu: float,
    chat: complex,
    ells: tuple[int, ...] = (40, 60, 80, 120, 160, 240),
) -> dict[str, complex | float]:
    """Fit di A/L^2 = A0 + A1/L + A2/L^2 a mu e chat fissati."""
    scales = []
    values = []
    for ell in ells:
        large_l = ell + 0.5
        m = int(round(mu * large_l))
        eigenvalue = spheroidal_eigenvalue(ell, m, chat * large_l)
        scales.append(large_l)
        values.append(eigenvalue / large_l**2)
    scales = np.array(scales)
    values = np.array(values)
    design = np.vstack([np.ones_like(scales), 1.0 / scales, 1.0 / scales**2]).T
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    residual = float(np.max(np.abs(design @ coefficients - values)))
    a0, a1, a2 = coefficients
    return {
        "A0": complex(a0),
        "A1": complex(a1),
        "A2": complex(a2),
        "residual": residual,
        "langer_shift": complex(a1 / (2.0 * a0)),
    }


def main() -> None:
    print("=== Controllo: caso sferico, chat = 0 ===")
    control = eikonal_fit(0.5, 0.0 + 0.0j)
    print(f"  A0 = {control['A0']:.10f}")
    print(f"  A1 = {control['A1']:.3e}   (atteso 0: Langer)")
    print(f"  A2 = {control['A2']:.8f}   (atteso -0.25)")
    print()

    print("=== chat COMPLESSA: c = a*omega con omega di QNM ===")
    print("  mu    chat                 A1                     |A1|      delta_Langer")
    cases = [
        (0.5, 0.4 - 0.04j),
        (0.8, 0.6 - 0.06j),
        (0.3, 0.9 - 0.09j),
        (0.9, 0.75 - 0.05j),
        (0.6, 0.5 - 0.15j),
    ]
    shifts = []
    for mu, chat in cases:
        fit = eikonal_fit(mu, chat)
        shifts.append(fit["langer_shift"])
        print(
            f"  {mu:<5} {chat!s:<18}  {fit['A1'].real:+.6f}{fit['A1'].imag:+.6f}i   "
            f"{abs(fit['A1']):.6f}   {fit['langer_shift'].real:+.5f}{fit['langer_shift'].imag:+.5f}i"
            f"   (res {fit['residual']:.1e})"
        )

    print()
    print("=== Non rimovibilita' con una sostituzione di Langer ===")
    array = np.array(shifts)
    spread = float(np.max(np.abs(array - array.mean())))
    print(f"  delta richiesti: dispersione attorno alla media = {spread:.5f}")
    print(f"  rapporto dispersione/media = {spread / abs(array.mean()):.1%}")
    print("  Un solo delta costante potrebbe assorbire A1 solo se questa fosse ~0.")

    print()
    print("=== Scalatura fisica: Im(chat) e' essa stessa O(eps) ===")
    print("  A overtone fissato Im(chat) = -eps*a(n+1/2)lambda: si verifica")
    print("  che rimpicciolire Im(chat) NON annulla A1, che resta reale e finito.")
    print("   Im(chat)        A1")
    for imaginary in (0.0, -0.02, -0.05, -0.10):
        fit = eikonal_fit(0.5, 0.4 + 1j * imaginary)
        print(f"   {imaginary:+.2f}      {fit['A1'].real:+.6f}{fit['A1'].imag:+.6f}i")


if __name__ == "__main__":
    main()
