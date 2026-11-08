"""Bounded constrained-relaxation probe (DISCOVERY_LAYER_DESIGN_2026-10-13 §4, v1).

At each target value of the intent coordinate q(x) = sum_k c_k r_k the structure is
relaxed on the unbiased PES with q held exactly (a FixConstraint on the bond combination); every
other internal degree of freedom is free. The energy recorded is the PES energy --
there is no penalty term -- and ``reaction_force`` is dE*/dq read from the
unconstrained forces, which at a constrained minimum is the Lagrange multiplier.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase.constraints import FixConstraint
from ase.optimize import BFGS

from .ledger import BudgetExhausted


@dataclass(frozen=True)
class Intent:
    intent_id: str
    pairs: tuple[tuple[int, int], ...]
    coefs: tuple[float, ...]

    def value(self, positions) -> float:
        return float(
            sum(
                c * np.linalg.norm(positions[j] - positions[i])
                for (i, j), c in zip(self.pairs, self.coefs)
            )
        )

    def gradient(self, positions) -> np.ndarray:
        grad = np.zeros_like(positions, dtype=float)
        for (i, j), c in zip(self.pairs, self.coefs):
            v = positions[j] - positions[i]
            u = v / np.linalg.norm(v)
            grad[j] += c * u
            grad[i] -= c * u
        return grad


# BFGS trust radius per step. The default 0.2 A lets the first Hessian update overshoot
# into the overlap guard on a soft, nearly flat surface (Ar3 LJ test: 0.2 -> gate at step 6,
# 0.1 / 0.05 / 0.02 all converge). Not a result gate: it only shortens trial steps.
_MAX_STEP_A = 0.05


class FixBondCombo(FixConstraint):
    """Hold q(x) = sum_k c_k r_k at ``target`` exactly.

    ASE's FixInternals projects out six rigid-body modes unconditionally, which is
    singular for linear systems (a diatomic has only five rigid modes, and the axial
    rotation vector is zero), so it cannot hold a bond length on a diatomic. q is
    translation- and rotation-invariant, so projecting the forces against grad q and
    restoring q by Newton steps along grad q is the complete constraint.
    """

    def __init__(self, intent, target, tol=1e-12, max_iter=50):
        self.intent = intent
        self.target = float(target)
        self.tol = tol
        self.max_iter = max_iter

    def adjust_positions(self, atoms, new):
        for _ in range(self.max_iter):
            err = self.intent.value(new) - self.target
            if abs(err) < self.tol:
                break
            g = self.intent.gradient(new)
            new -= err * g / np.vdot(g, g)

    def adjust_forces(self, atoms, forces):
        g = self.intent.gradient(atoms.positions)
        forces -= np.vdot(forces, g) / np.vdot(g, g) * g

    def get_removed_dof(self, atoms):
        return 1

    def index_shuffle(self, atoms, ind):
        raise NotImplementedError

    def todict(self):
        return {"name": "FixBondCombo", "kwargs": {"target": self.target}}


@dataclass(frozen=True)
class ScanSettings:
    deltas_A: tuple[float, ...]
    fmax_eV_A: float
    max_steps_per_point: int
    max_rise_eV: float
    reuse_curvature: bool = True

    def __post_init__(self):
        if not self.deltas_A or self.deltas_A[0] != 0.0:
            raise ValueError("deltas_A must start at 0.0 (the intent's own start)")


@dataclass
class ScanPoint:
    delta_A: float
    q_target: float
    q_actual: float
    energy_eV: float
    reaction_force: float
    fmax_projected: float
    converged: bool
    steps: int
    edges: list
    positions: np.ndarray
    forces_pes: np.ndarray | None = None


@dataclass
class ScanResult:
    intent_id: str
    q0: float
    points: list
    status: str
    rejected_positions: np.ndarray | None = None


def _feasible(intent, q0, deltas, numbers, config):
    """Design §3.2: reject an unreachable constraint before any evaluation. Every distance
    is bounded below by the existing overlap guard, r_k >= min_distance_ratio *
    (R_i + R_j), so with all coefficients positive q cannot go below sum c_k floor_k
    (all negative: above it). Mixed signs leave q unbounded both ways."""
    from ase.data import covalent_radii

    floors = [
        config.min_distance_ratio * (covalent_radii[numbers[i]] + covalent_radii[numbers[j]])
        for i, j in intent.pairs
    ]
    bound = sum(c * f for c, f in zip(intent.coefs, floors))
    if all(c > 0 for c in intent.coefs):
        return all(q0 + d > bound for d in deltas)
    if all(c < 0 for c in intent.coefs):
        return all(q0 + d < bound for d in deltas)
    return True


def feasible_prefix(intent, q0, deltas, numbers, config):
    """The longest leading part of ``deltas`` inside the overlap-floor domain (design
    §3.2). ``deltas`` start at 0 and move away monotonically, so once a target is
    infeasible every later one is too. Fewer than two points left means no scan."""
    k = len(deltas)
    while k > 0 and not _feasible(intent, q0, deltas[:k], numbers, config):
        k -= 1
    return tuple(deltas[:k])


def reverse_plan(points, jump, settings):
    """Start positions and settings of the reverse scan that checks a jump for
    hysteresis: from the post-jump point back to the pre-jump target."""
    before, after = points[jump["before"]], points[jump["after"]]
    back = ScanSettings(
        deltas_A=(0.0, before.delta_A - after.delta_A),
        fmax_eV_A=settings.fmax_eV_A,
        max_steps_per_point=settings.max_steps_per_point,
        max_rise_eV=settings.max_rise_eV,
        reuse_curvature=settings.reuse_curvature,
    )
    return after.positions.copy(), back


def constrained_scan(
    atoms,
    intent,
    settings,
    factory,
    config,
    ledger,
    task,
    attempt=0,
    reserve=None,
    artifact=None,
):
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.state import encode

    work = atoms.copy()
    work.set_constraint()
    q0 = intent.value(work.positions)
    points = []
    if not _feasible(intent, q0, settings.deltas_A, work.numbers, config):
        return ScanResult(intent.intent_id, q0, points, "infeasible")
    work.calc = GuardedCalculator(factory(), config)
    status = "completed"
    e_start = None
    rejected = None
    hessian = None
    try:
        with ledger.segment(task=task, stage="probe", attempt=attempt, reserve=reserve) as seg:
            seg.artifact = artifact
            try:
                for delta in settings.deltas_A:
                    target = q0 + delta
                    work.set_constraint(FixBondCombo(intent, target))
                    # BFGS only applies the constraint on a step; a point with zero projected
                    # force would otherwise never be moved onto the target.
                    work.set_positions(work.positions)
                    opt = BFGS(work, logfile=None, maxstep=_MAX_STEP_A)
                    if settings.reuse_curvature and hessian is not None:
                        # Only transfer curvature. A new optimizer has no secant
                        # history, so the target move never enters a BFGS update.
                        # Reset the constraint-normal block to ASE's own alpha;
                        # retain curvature on the NEW tangent space.
                        normal = intent.gradient(work.positions).ravel()
                        normal /= np.linalg.norm(normal)
                        nn = np.outer(normal, normal)
                        tangent = np.eye(len(normal)) - nn
                        opt.H0 = tangent @ hessian @ tangent + opt.alpha * nn
                    converged = bool(
                        opt.run(fmax=settings.fmax_eV_A, steps=settings.max_steps_per_point)
                    )
                    # ASE versions expose H as either an attribute or a property
                    # of an initialized state. Zero-step points may have neither.
                    learned = getattr(opt, "H", None)
                    hessian = learned.copy() if converged and learned is not None else None
                    energy = float(work.get_potential_energy())
                    f_pes = work.get_forces(apply_constraint=False)
                    g = intent.gradient(work.positions)
                    point = ScanPoint(
                        delta_A=float(delta),
                        q_target=float(target),
                        q_actual=intent.value(work.positions),
                        energy_eV=energy,
                        reaction_force=-float(np.vdot(f_pes, g) / np.vdot(g, g)),
                        fmax_projected=float(np.max(np.linalg.norm(work.get_forces(), axis=1))),
                        converged=converged,
                        steps=int(opt.nsteps),
                        edges=sorted(map(list, encode(work, config.bond_scale).edges)),
                        positions=work.positions.copy(),
                        forces_pes=f_pes.copy(),
                    )
                    points.append(point)
                    e_start = energy if e_start is None else e_start
                    if energy - e_start > settings.max_rise_eV:
                        status = "censored_energy"
                        break
            except BudgetExhausted:
                status = "budget_exhausted"
            except GateRejected as exc:
                status = f"gate_rejected:{exc.code}"
                rejected = exc.atoms.positions.copy()
            seg.outcome = status
    except BudgetExhausted:  # the reservation itself was refused: no segment opened
        status = "budget_exhausted"
    finally:
        work.set_constraint()
    return ScanResult(intent.intent_id, q0, points, status, rejected)


def read_scan(points, noise_eV):
    """Release triggers from a finished scan (§4.3). Pure: no evaluations."""
    triggers = []
    energies = [p.energy_eV for p in points]
    for k in range(1, len(points) - 1):
        if not all(p.converged for p in points[: k + 2]):
            break
        # F3 (TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24 §4): the drop must come from the
        # contiguous converged run right after the peak; a later unconverged point is kept
        # in the scan but is not evidence of a descent.
        end = k + 1
        while end + 1 < len(points) and points[end + 1].converged:
            end += 1
        rises = energies[k] - energies[0] > noise_eV
        falls = energies[k] - min(energies[k + 1 : end + 1]) > noise_eV
        flips = points[k - 1].reaction_force * points[k + 1].reaction_force < 0
        if (
            rises
            and falls
            and flips
            and energies[k] >= energies[k - 1]
            and energies[k] >= energies[k + 1]
        ):
            triggers.append({"kind": "ridge_suspected", "peak": k, "release": k + 1})
            break
    for k in range(1, len(points)):
        a, b = points[k - 1], points[k]
        if (
            a.converged
            and b.converged
            and a.edges != b.edges
            and b.energy_eV < a.energy_eV - noise_eV
        ):
            triggers.append({"kind": "jump", "before": k - 1, "after": k})
    return triggers


def bracket_jump(points, jump, reverse_point):
    """A jump is bracketed when, scanned back to ``before``, the structure stays on the
    graph it jumped to: two branches coexist at that target (hysteresis)."""
    before, after = points[jump["before"]], points[jump["after"]]
    if reverse_point.converged and reverse_point.edges == after.edges != before.edges:
        return {
            "kind": "branch_bracketed",
            "reactant_side": jump["before"],
            "product_side_positions": reverse_point.positions,
        }
    return None


def reverse_check(
    atoms_template,
    points,
    jump,
    intent,
    settings,
    factory,
    config,
    ledger,
    task,
    attempt=0,
    reserve=None,
    artifact=None,
):
    positions, back = reverse_plan(points, jump, settings)
    start = atoms_template.copy()
    start.positions = positions
    result = constrained_scan(
        start,
        intent,
        back,
        factory,
        config,
        ledger,
        task,
        attempt=attempt,
        reserve=reserve,
        artifact=artifact,
    )
    if result.status != "completed" or len(result.points) < 2:
        return None
    return bracket_jump(points, jump, result.points[-1])


def release(
    positions,
    template,
    source,
    factory,
    config,
    ledger,
    task,
    seed,
    attempt=0,
    reserve=None,
    artifact=None,
):
    """Drop every external constraint, quench on the unbiased PES, confirm the minimum,
    classify against the source. Static: no momenta, no claim of a thermal crossing.

    ``reserve`` is held back for each of the two sequential segments (quench,
    minimum_check) in turn. It guarantees availability, it does not cap spending: forwards
    draw it down first, then continue against the global remaining budget. The first
    segment's unused part is returned before the second opens."""
    from prrs.io import snapshot
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.runner import confirm_minimum, quench
    from prrs.state import classify

    release_positions = np.asarray(positions, float).copy()
    constraints_removed = [repr(c) for c in template.constraints]
    atoms = template.copy()
    atoms.set_constraint()
    atoms.positions = release_positions.copy()
    atoms.set_momenta(np.zeros((len(atoms), 3)))
    out = {
        "status": None,
        "label": None,
        "fingerprint": None,
        "endpoint": None,
        "release_positions": release_positions,
        "release_forces_pes": None,
        "seed": seed,
        "constraints_removed": constraints_removed,
    }
    try:
        guard = GuardedCalculator(factory(), config)
        atoms.calc = guard
        seg_args = {"task": task, "attempt": attempt, "reserve": reserve}
        with ledger.segment(stage="quench", **seg_args) as seg:
            seg.artifact = artifact
            try:
                out["release_forces_pes"] = atoms.get_forces().copy()
                converged, _ = quench(atoms, guard, config)
            except BudgetExhausted:
                seg.outcome = "budget_exhausted"
                raise
            except GateRejected as exc:
                seg.outcome = f"gate_rejected:{exc.code}"
                raise
            seg.outcome = "converged" if converged else "quench_failed"
            if converged:
                endpoint = snapshot(atoms, guard)
        if not converged:
            out["status"] = "quench_failed"
            return out
        with ledger.segment(stage="minimum_check", **seg_args) as seg:
            seg.artifact = artifact
            try:
                confirmed, info = confirm_minimum(endpoint, factory, config, seed)
            except BudgetExhausted:
                seg.outcome = "budget_exhausted"
                raise
            except GateRejected as exc:
                seg.outcome = f"gate_rejected:{exc.code}"
                raise
            reason = str((info or {}).get("reason", ""))
            gated = not confirmed and reason.startswith("gate_")
            if gated:
                seg.outcome = f"gate_rejected:{reason[len('gate_') :]}"
            else:
                seg.outcome = "minimum" if confirmed else "not_minimum"
        if gated:
            out["status"] = seg.outcome
            return out
        if not confirmed:
            out["status"] = "not_minimum"
            return out
        label, fingerprint = classify(source, endpoint, config)
        out.update(
            status="c1_candidate", label=label, fingerprint=fingerprint, endpoint=endpoint
        )
    except BudgetExhausted:
        out["status"] = "budget_exhausted"
    except GateRejected as exc:
        out["status"] = f"gate_rejected:{exc.code}"
    return out
