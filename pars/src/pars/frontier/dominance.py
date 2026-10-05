"""Nominal and uncertainty-aware dominance (design sections 11-12).

Objectives are dictionaries with availability masks. Different masks are
incomparable: pair-specific intersections can create cyclic dominance and an
empty frontier. Certain dominance requires complete uncertainty intervals for
both candidates, or an explicit tolerance policy when neither has intervals.

Tolerance only ever makes a claim harder. Nominal dominance needs no
objective worse and one better by MORE than the tolerance. Over a common set
of available objectives that is a strict partial order (irreflexive,
antisymmetric, transitive) for every tolerance, so the frontier of a finite
set is never empty and every dominated candidate is dominated by a frontier
member. Certain dominance needs interval separation
by at least the tolerance on every objective.
"""

from __future__ import annotations

ObjectiveVector = dict[str, float | None]


def comparable_keys(a: ObjectiveVector, b: ObjectiveVector) -> list[str]:
    available_a = {k for k, value in a.items() if value is not None}
    available_b = {k for k, value in b.items() if value is not None}
    return sorted(available_a) if available_a == available_b else []


def dominates(a: ObjectiveVector, b: ObjectiveVector, tolerance: float = 0.0) -> bool:
    """a <= b everywhere and a < b - tolerance somewhere, over jointly
    available objectives with identical availability masks. tolerance=0 is strict Pareto dominance; a larger
    tolerance treats smaller improvements as indistinguishable, which can only
    widen the frontier."""
    keys = comparable_keys(a, b)
    if not keys:
        return False
    weakly = all(a[k] <= b[k] for k in keys)
    strictly = any(a[k] < b[k] - tolerance for k in keys)
    return weakly and strictly


def certainly_dominates(
    a: ObjectiveVector,
    b: ObjectiveVector,
    sigma_a: dict[str, float] | None = None,
    sigma_b: dict[str, float] | None = None,
    tolerance: float = 0.0,
) -> bool:
    """Conservative interval test (design section 12): on every comparable
    available objective the upper bound of a must sit at least ``tolerance``
    below the lower bound of b. Hard deletion requires this, not the nominal
    point estimates. A positive tolerance can explicitly license separation of
    point estimates only when neither candidate supplies uncertainty data;
    partial intervals never imply zero uncertainty on missing entries."""
    keys = comparable_keys(a, b)
    if not keys or not dominates(a, b, tolerance):
        return False
    if not sigma_a and not sigma_b:
        if tolerance <= 0:
            return False
        sigma_a = dict.fromkeys(keys, 0.0)
        sigma_b = dict.fromkeys(keys, 0.0)
    elif (
        sigma_a is None
        or sigma_b is None
        or any(sigma_a.get(k) is None or sigma_b.get(k) is None for k in keys)
    ):
        return False
    return all(a[k] + float(sigma_a[k]) + tolerance <= b[k] - float(sigma_b[k]) for k in keys)
