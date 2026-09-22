"""Conditional completeness certificates for basin membership (PLAN item 3, 2026-09-20).

`search.py` samples basins; it never claims a basin list is complete. This module is
the machinery for the one honest thing that CAN be said in the other direction:
given a covering radius delta, a measured gradient bound L, a subspace V and a ball
radius R, and a membership table that says every sample point quenches back to A,
record the four quantities and check the three conditions of the draft proposition
(decision 4, DECISIONS_PENDING 2026-09-20). It does NOT issue certificates: the
proposition is a draft awaiting ratification, so every text it renders is stamped
DRAFT, and that stamp only comes off by explicit user ratification.

What the code owns vs what the mathematics owns, stated once:

  code:   delta, L, V, R, the membership table, and whether the three recorded
          conditions hold on the measurements.
  maths:  whether (i)+(ii)+(iii) actually imply "no sublevel component below X
          outside A in B cap V". That implication is the thing being ratified.
          Nothing here silently upgrades a measured L (a sup over sampled pairs,
          hence a LOWER bound of the true Lipschitz constant) into a proven upper
          bound: callers may supply `lipschitz_upper` from an analytic Hessian norm,
          and L_used is the max of the two, declared in the record.

Cost is exponential in d, so V is built by the caller (d <= 3 enforced) and its
construction travels inside every record -- a certificate over an unjustified
subspace certifies nothing and must not look like one.
"""

import numpy as np

from .reliability import GuardedCalculator
from .runner import quench

MAX_SUBSPACE_DIM = 3


def cover(centre, basis, radius, delta):
    """A d-dimensional grid over the ball with covering radius at most delta.

    Axis-aligned grid with spacing h = 2*delta/sqrt(d): the farthest any point of a
    cell can be from its centre is h*sqrt(d)/2 = delta, so every ball point lies
    within delta of some grid point that is itself kept only if inside the ball.
    Deterministic -- no low-discrepancy sequence, nothing seeded.
    """
    basis = np.asarray(basis, dtype=float)
    if basis.ndim != 2 or basis.shape[1] != len(centre):
        raise ValueError("basis must be d columns over the flat coordinate space")
    d = basis.shape[0]
    if not 1 <= d <= MAX_SUBSPACE_DIM:
        raise ValueError(f"subspace dimension {d} exceeds the enforced cap {MAX_SUBSPACE_DIM}")
    if delta <= 0 or radius <= 0:
        raise ValueError("radius and delta must be positive")
    h = 2.0 * delta / np.sqrt(d)
    axes = [
        np.arange(np.ceil(-(radius + delta) / h), np.floor((radius + delta) / h) + 1) * h
        for _ in range(d)
    ]
    points = []
    for combo in np.array(np.meshgrid(*axes, indexing="ij")).reshape(d, -1).T:
        # Keep grid points a little OUTSIDE the ball too: a ball point near the
        # boundary has its nearest grid point outside, and filtering that away is
        # how a cover silently stops covering.
        if np.linalg.norm(combo) <= radius + delta + 1e-12:
            points.append(np.asarray(centre) + basis.T @ combo)
    if not points:
        raise ValueError("empty cover: delta or radius degenerate")
    return points


def membership(point, factory, config, basin_minimum, basin_rmsd_A):
    """Quench from `point`; the answer is True only if it lands on the basin minimum.

    Fail-closed by construction: a quench that does not converge, or a guard gate
    that rejects the geometry, is recorded as its own outcome and is NEVER True.
    "We could not even ask" is not evidence of membership.
    """
    atoms = basin_minimum.copy()
    atoms.set_positions(np.asarray(point).reshape(atoms.positions.shape))
    atoms.set_momenta(np.zeros(atoms.positions.shape))
    guard = GuardedCalculator(factory(), config)
    atoms.calc = guard
    try:
        converged, _ = quench(atoms, guard, config)
    except Exception as exc:  # noqa: BLE001 - gate rejections are data here, not failures
        return {"member": False, "outcome": f"gate_{type(exc).__name__}", "detail": str(exc)}
    if not converged:
        return {"member": False, "outcome": "quench_unconverged"}
    delta = atoms.get_positions() - basin_minimum.get_positions()
    rmsd = float(np.sqrt((delta**2).sum(axis=1).mean()))
    return {
        "member": rmsd <= basin_rmsd_A,
        "outcome": "basin" if rmsd <= basin_rmsd_A else "other_basin",
        "rmsd_A": rmsd,
    }


def gradient_bound(points, gradients, min_separation=1e-6):
    """Measured Lipschitz estimate: sup over pairs of ||dE difference|| / distance.

    This is a sup over SAMPLED pairs, so it is a lower bound of the true Lipschitz
    constant -- declared in the record, never silently promoted. Pairs closer than
    `min_separation` are skipped: at separations the covering argument never uses,
    floating-point noise in the gradient difference dominates the ratio and the
    bound reports noise, not curvature.
    """
    points = np.asarray(points, dtype=float)
    gradients = np.asarray(gradients, dtype=float)
    if len(points) < 2:
        raise ValueError("gradient bound needs at least two points")
    flat_p = points.reshape(len(points), -1)
    flat_g = gradients.reshape(len(gradients), -1)
    best = 0.0
    for i in range(len(flat_p)):
        for j in range(i + 1, len(flat_p)):
            distance = np.linalg.norm(flat_p[i] - flat_p[j])
            if distance < max(min_separation, 1e-12):
                continue
            best = max(best, np.linalg.norm(flat_g[i] - flat_g[j]) / distance)
    return float(best)


def audit_basin(
    minimum,
    factory,
    config,
    basis,
    radius,
    delta,
    energy_ceiling_X,
    gap_w,
    lipschitz_upper=None,
):
    """Assemble the certificate record and evaluate the draft proposition's checks.

    The three conditions, numbered as in the draft (decision 4):

      (i)   every sample of the annulus [gap_w, radius] quenches back to A;
      (ii)  L_used * delta < g_min, with g_min the smallest measured ||grad E||
            over the annulus samples -- the covering argument's engine: any point of
            the annulus lies within delta of a sample, so its gradient cannot vanish;
      (iii) every boundary-shell sample (distance in [radius - delta, radius]) has
            E > X -- a sublevel component that reaches the boundary could leave the
            ball without a critical point, and then nothing here rules it out.

    The four quantities delta / L / V / R are all recorded; any missing one is a
    ValueError rather than a quieter record (fail closed).
    """
    for name, value in (
        ("radius", radius),
        ("delta", delta),
        ("gap_w", gap_w),
        ("energy_ceiling_X", energy_ceiling_X),
    ):
        if value is None:
            raise ValueError(f"certificate quantity {name} is missing; refusing (fail closed)")
    basis = np.asarray(basis, dtype=float)
    d = basis.shape[0]
    points = cover(minimum.get_positions().ravel(), basis, radius, delta)

    table = []
    energies, gradients, distances = [], [], []
    for point in points:
        atoms = minimum.copy()
        atoms.set_positions(np.asarray(point).reshape(atoms.positions.shape))
        atoms.calc = factory()
        energies.append(float(atoms.get_potential_energy()))
        gradients.append(atoms.get_forces().ravel().copy())
        centre = minimum.get_positions().ravel()
        offset = np.asarray(point) - centre
        # distance measured inside V (the certificate never speaks outside its subspace)
        distances.append(float(np.linalg.norm(basis @ offset)))

    for point, energy, gradient, distance in zip(points, energies, gradients, distances):
        if distance < gap_w:
            table.append({"in": "core", "member": True, "distance": distance})
            continue
        outcome = membership(point, factory, config, minimum, config.basin_rmsd_A)
        table.append(
            {
                "in": "annulus" if distance < radius - 1e-12 else "boundary_shell",
                **outcome,
                "energy_eV": energy,
                "gradient_norm": float(np.linalg.norm(gradient)),
                "distance": distance,
            }
        )

    annulus = [row for row in table if row["in"] in ("annulus", "boundary_shell")]
    g_min = min((row["gradient_norm"] for row in annulus), default=None)
    l_measured = gradient_bound(
        [p for p, row in zip(points, table) if row["in"] != "core"],
        [g for g, row in zip(gradients, table) if row["in"] != "core"],
    )
    l_used = max(l_measured, lipschitz_upper or 0.0)

    conditions = {
        "(i) membership all basin": all(row["member"] for row in table),
        "(ii) L_used*delta < g_min": (g_min is not None and l_used * delta < g_min),
        "(iii) boundary shell above X": all(
            row["energy_eV"] > energy_ceiling_X
            for row in table
            if row["in"] == "boundary_shell"
        ),
    }
    return {
        "subspace": {
            "dimension": d,
            "basis": basis.tolist(),
            "construction": "caller-supplied; a certificate over an unjustified "
            "subspace certifies nothing",
            "cap": MAX_SUBSPACE_DIM,
        },
        "radius_R": radius,
        "delta": delta,
        "gap_w": gap_w,
        "energy_ceiling_X": energy_ceiling_X,
        "l_measured": l_measured,
        "l_used": l_used,
        "g_min": g_min,
        "samples": len(points),
        "membership_table": table,
        "conditions": conditions,
        "verdict": "conditions_met" if all(conditions.values()) else "refused",
    }


_DRAFT_STAMP = (
    "DRAFT -- proposition not ratified (decision 4, 2026-09-20). This text is a "
    "measurement record, not an issued certificate; nothing may cite it as one."
)


def certificate_text(record):
    """Render a record as text. Always DRAFT-stamped until the proposition is ratified."""
    lines = [
        _DRAFT_STAMP,
        "",
        f"verdict: {record['verdict']}",
        f"V: d={record['subspace']['dimension']} (cap {record['subspace']['cap']}), "
        f"construction: {record['subspace']['construction']}",
        f"R = {record['radius_R']} A, delta = {record['delta']} A, "
        f"X = {record['energy_ceiling_X']} eV, gap w = {record['gap_w']} A",
        f"L measured = {record['l_measured']:.6f}, L used = {record['l_used']:.6f} eV/A^2, "
        f"g_min = {record['g_min']} eV/A",
        f"samples: {record['samples']}",
    ]
    for name, holds in record["conditions"].items():
        lines.append(f"  {name}: {'met' if holds else 'NOT met'}")
    return "\n".join(lines)
