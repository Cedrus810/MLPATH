"""Mechanism-level path descriptors and the path-space distance (design 15).

V1 descriptor: the ordered sequence of per-frame bond sets that PRRS recorded
in observations.jsonl. Bond perception is authoritative on the PRRS side; PARS
only consumes it and never re-derives connectivity.

Local cost between two frames is the Jaccard distance between bond sets; the
sequence distance is the discrete Frechet distance, which respects the ORDER
of bond-change events: two paths visiting the same states in a different
order (proton transfer then cleavage vs. the reverse) are far apart even when
their frame sets coincide.
"""

from __future__ import annotations

Signature = frozenset[tuple[int, int]]


def signature_distance(a: Signature, b: Signature) -> float:
    union = a | b
    if not union:
        return 0.0
    return 1.0 - len(a & b) / len(union)


def frechet_path_distance(a: tuple[Signature, ...], b: tuple[Signature, ...]) -> float:
    """Discrete Frechet distance between two signature sequences."""
    if not a or not b:
        raise ValueError("path distance requires non-empty signature sequences")
    rows, cols = len(a), len(b)
    cost = [[0.0] * cols for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            local = signature_distance(a[i], b[j])
            if i == 0 and j == 0:
                value = local
            elif i == 0:
                value = max(cost[0][j - 1], local)
            elif j == 0:
                value = max(cost[i - 1][0], local)
            else:
                value = max(min(cost[i][j - 1], cost[i - 1][j], cost[i - 1][j - 1]), local)
            cost[i][j] = value
    return cost[-1][-1]


def recrossing_count(signatures: tuple[Signature, ...]) -> int:
    """Number of returns to the source basin after the first departure.

    Basin membership is approximated by exact bond-set equality with the first
    frame (PRRS bond perception already carries hysteresis). A clean
    A->B trajectory has zero recrossings; A->B->A->B has one.
    """
    if not signatures:
        raise ValueError("recrossing requires a non-empty signature sequence")
    source = signatures[0]
    outside = [signature != source for signature in signatures]
    return sum(1 for left, arrived in zip(outside, outside[1:]) if left and not arrived)
