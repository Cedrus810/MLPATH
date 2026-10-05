from __future__ import annotations

import pytest

from pars.paths.fingerprint import (
    frechet_path_distance,
    recrossing_count,
    signature_distance,
)

A = frozenset({(0, 1)})
AB = frozenset({(0, 1), (1, 2)})


def test_signature_distance_is_jaccard():
    assert signature_distance(A, A) == 0.0
    assert signature_distance(A, AB) == pytest.approx(1 / 2)
    assert signature_distance(A, frozenset({(1, 2)})) == 1.0


def test_frechet_identical_sequences_are_zero():
    seq = (A, A, AB, AB)
    assert frechet_path_distance(seq, seq) == 0.0


def test_frechet_is_symmetric():
    first = (A, A, AB, AB)
    second = (A, AB, AB, AB)
    assert frechet_path_distance(first, second) == pytest.approx(
        frechet_path_distance(second, first)
    )


def test_frechet_respects_event_order():
    # same bond sets visited in opposite order are NOT the same mechanism:
    # the distance is bounded below by the A<->AB frame cost, which is > 0
    forward = (A, A, AB, AB)
    reverse = (AB, AB, A, A)
    assert frechet_path_distance(forward, reverse) == pytest.approx(0.5)
    assert frechet_path_distance((A, AB), (AB, A)) == pytest.approx(0.5)
    assert frechet_path_distance(forward, forward) == 0.0


def test_frechet_empty_sequences_rejected():
    with pytest.raises(ValueError, match="non-empty"):
        frechet_path_distance((), (A,))


def test_recrossing_on_signature_sequences():
    assert recrossing_count((A, A, AB, AB)) == 0
    assert recrossing_count((A, A, AB, A, AB, A)) == 2
    assert recrossing_count((A,)) == 0
