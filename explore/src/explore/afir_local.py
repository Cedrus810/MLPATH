"""MIN, TS and connection from an AFIR path with E/F only (design
docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md §4; plan docs/PLAN_AFIR_MC_2026-11-14.md).

This is the project's own part. GRRM turns AFIR-path minima / maxima into MINs / TSs with
RFO and exact Hessians and connects with a Hessian-started IRC; here:

- MIN: unbiased E/F minimization of a path minimum, confirmed by perturbed re-descents
  that return to the same basin (no Hessian).
- TS: the AFIR path is a fixed-force path, which jumps the part of the rho_AB branch where
  d2E*/drho2 < 0, so its top is only near a TS -- and far from it when gamma is much
  larger than the barrier (Chem. Rec. 2016 §2.2; on the TwoSaddle bench the path top sits
  0.66 A in rho past the saddle). Path frames are relaxed at their own fixed rho_AB
  (displacement control, unbiased E/F) into ScanPoints, starting at the top and walking
  along the path towards the side the top's reaction force points to, until the reaction
  force dE/drho changes sign; that pair brackets a maximum on the branch and the HF E/F
  kernel (`continue_branches_ef` -> `refine_ts_ef`) locates it.
- Branch jump (design §4.2; user 2026-11-15): when the kernel on the push coordinate ends
  in ``unresolved:hysteresis``, each jump's two sides (this branch's last structure and
  the other branch at the same q) define the distance-space displacement coordinate q_d
  (`afir.DisplacementCoordinate`); points between them are relaxed at fixed q_d, the
  first + -> - pair of dE/dq_d is handed to the same kernel (`rebracket_jump`). The
  concerted motion of all atoms across the jump thus enters the locating coordinate.
- Connection: E/F descent from both sides along the search's own direction
  (`branch_starts` + `trace_sides_ef`, as `d2._hc_run`).
"""

from __future__ import annotations

import numpy as np
from ase.optimize import BFGS

from .ledger import BudgetExhausted
from .probe import ScanPoint
from .ts_branch import BranchSettings, continue_branches_ef
from .ts_connect import ConnectSettings, branch_starts, trace_sides_ef
from .ts_ef import EFSettings, _align, _Constraint, _relax, _Stop

TS_MIN_STEP_A = 0.05  # smallest rho change between two relaxed path frames
TS_MAX_POINTS = 24  # relaxed frames per bracket search


def path_minima(frames, noise_eV):
    """Indices of energy minima along the path (each confirmed by a rise > noise on
    both sides) and the last frame."""
    e = [f["energy_eV"] for f in frames]
    out = []
    for k in range(1, len(e) - 1):
        if e[k] <= e[k - 1] and e[k] <= e[k + 1]:
            left, right = max(e[:k]), max(e[k + 1 :])
            if left - e[k] > noise_eV and right - e[k] > noise_eV:
                out.append(k)
    if frames and (len(frames) - 1) not in out:
        out.append(len(frames) - 1)
    return out


def relax_min(template, positions, factory, config, ledger, task, *, rng, attempt=0,
              n_perturb=3, perturb_A=0.03, max_steps=1000):  # fmt: skip
    """Unbiased E/F minimum from ``positions`` and its perturbation check."""
    from prrs.state import aligned_rmsd, encode

    def descend(x):
        a = template.copy()
        a.set_constraint()
        a.positions = np.asarray(x, float)
        a.calc = factory()
        ok = BFGS(a, logfile=None, maxstep=0.1).run(
            fmax=config.quench_fmax_eV_A, steps=max_steps
        )
        return a, bool(ok)

    out = {"status": "unresolved", "positions": None}
    try:
        with ledger.segment(task=task, stage="afir_min", attempt=attempt) as seg:
            a, ok = descend(positions)
            if not ok:
                out["status"] = "descent_not_converged"
            else:
                e = float(a.get_potential_energy())
                g = encode(a, config.bond_scale, active=config.active_atoms).edges
                trials = []
                for _ in range(n_perturb):
                    b, okb = descend(
                        a.positions + rng.normal(0.0, perturb_A, a.positions.shape)
                    )
                    same = (
                        okb
                        and encode(b, config.bond_scale, active=config.active_atoms).edges == g
                        and abs(float(b.get_potential_energy()) - e) <= config.basin_energy_eV
                        and aligned_rmsd(a, b, config.active_atoms) <= config.basin_rmsd_A
                    )
                    trials.append(bool(same))
                out = {
                    "status": "confirmed" if all(trials) else "unstable",
                    "positions": a.positions.copy(),
                    "energy_eV": e,
                    "edges": sorted([int(i), int(j)] for i, j in g),
                    "trials": trials,
                }
            seg.outcome = out["status"]
    except BudgetExhausted:
        out["status"] = "budget_exhausted"
    return out


def _scan_point(template, coordinate, frame, factory, config, cfg, rho0):
    """Unbiased relaxation at the frame's own rho_AB; a ScanPoint for the HF kernel."""
    from prrs.state import encode

    a = template.copy()
    a.set_constraint()
    a.positions = np.asarray(frame["positions"], float)
    q = coordinate.value(a.positions)
    a.set_constraint(_Constraint(coordinate, q, tol=1e-8))
    a.set_positions(a.positions)
    a.calc = factory()
    opt = BFGS(a, logfile=None, maxstep=0.05)
    try:
        converged = bool(opt.run(fmax=cfg.transverse_fmax_eV_A, steps=cfg.max_relax_steps))
    except _Stop:
        return None
    f = a.get_forces(apply_constraint=False)
    g = coordinate.gradient(a.positions)
    rf = -float(np.vdot(f, g) / np.vdot(g, g))
    proj = f + rf * g
    a.set_constraint()
    return ScanPoint(
        delta_A=q - rho0,
        q_target=q,
        q_actual=coordinate.value(a.positions),
        energy_eV=float(a.get_potential_energy()),
        reaction_force=rf,
        fmax_projected=float(np.linalg.norm(proj, axis=1).max()),
        converged=converged,
        steps=int(opt.nsteps),
        edges=sorted([int(i), int(j)] for i, j in encode(a, config.bond_scale, active=config.active_atoms).edges),
        positions=a.positions.copy(),
        forces_pes=np.asarray(f, float),
    )  # fmt: skip


def ts_from_path(template, coordinate, frames, top, factory, config, ledger, task, *,
                 ef_settings, branch_settings, ceiling_eV, attempt=0):  # fmt: skip
    """Bracket and locate a TS from path frame ``top``; returns (status, points, branch).

    A maximum of the rho branch is where dE/drho changes from + to - with increasing rho.
    A fixed-force path reaches its energy top only after it has passed the TS region (it
    jumps the unstable part), so frames are relaxed from the top back towards the path
    start, each at least TS_MIN_STEP_A in rho from the last one, until two consecutive
    relaxed frames show that sign change; failing that, forward from the top. The reaction
    force at the top alone does not say where the TS is: on the product-side wall it
    points away from the TS (TwoSaddle bench, 2026-11-14)."""
    cfg = ef_settings
    rho0 = frames[0]["rho_A"]
    points, pair = [], None

    def scan(k):
        p = _scan_point(template, coordinate, frames[k], factory, config, cfg, rho0)
        if p is not None:
            points.append(p)
        return p

    def brackets(a, b):
        lo, hi = sorted((a, b), key=lambda p: p.q_target)
        return (
            a.converged
            and b.converged
            and hi.q_target - lo.q_target > 1e-6
            and lo.reaction_force > 0 > hi.reaction_force
        )

    try:
        with ledger.segment(task=task, stage="afir_ts_bracket", attempt=attempt) as seg:
            first = scan(top)
            for step in (-1, +1):
                if first is None:
                    break
                last, k = first, top
                while len(points) < TS_MAX_POINTS:
                    k += step
                    if not 0 <= k < len(frames):
                        break
                    if abs(frames[k]["rho_A"] - last.q_target) < TS_MIN_STEP_A:
                        continue
                    p = scan(k)
                    if p is None:
                        break
                    if brackets(last, p):
                        pair = tuple(sorted((last, p), key=lambda x: x.q_target))
                        break
                    last = p
                if pair is not None:
                    break
            seg.outcome = f"{len(points)} points"
    except BudgetExhausted:
        return "budget_exhausted", points, None
    if pair is None:
        return "no_bracket", points, None
    branch = continue_branches_ef(
        template.copy(),
        coordinate,
        pair,
        factory,
        config,
        ledger,
        task,
        branch_settings=branch_settings,
        settings=cfg,
        ceiling_eV=ceiling_eV,
        attempt=attempt,
    )
    return branch.status, points, branch


def link(
    template, branch, factory, config, ledger, task, *, ef_settings, sigma, seed, attempt=0
):
    """Both sides of a located candidate (as d2._hc_run): E/F descent and endpoint checks."""
    found = branch.refine
    frames = list(found["history"])
    for side in branch.sides:
        if side.get("status") == "maximum_bracketed":
            frames += side["points"]
    starts = branch_starts(
        frames,
        found["metrics"]["q_actual"],
        sigma,
        transverse_fmax=ef_settings.transverse_fmax_eV_A,
        descent_fmax=config.quench_fmax_eV_A,
    )
    return trace_sides_ef(
        template.copy(),
        found["positions"],
        found["energy_eV"],
        starts,
        found["direction_cartesian"],
        factory,
        config,
        ledger,
        task,
        settings=ConnectSettings(reaction_force_sigma=sigma, seed=seed),
        attempt=attempt,
    )


REBRACKET_POINTS = 5  # q_d values from the pre-jump to the post-jump structure


def rebracket_jump(template, branch, factory, config, ledger, task, *, ef_settings,
                   branch_settings, ceiling_eV=None, attempt=0, n_points=REBRACKET_POINTS):  # fmt: skip
    """Re-bracket an ``unresolved:hysteresis`` kernel result on the displacement coordinate.

    For each side that ended in a jump: q_d from its last accepted structure to the
    other branch at the same q (the opposite side's accepted point nearest in q; the
    kernel's post-jump frame if the opposite side has none), in distance space
    (`afir.DisplacementCoordinate`); points started on the straight line between the two
    (the second rigidly aligned onto the first); ``n_points`` evenly spaced q_d in [0, |d|], started on the
    straight line, each relaxed at fixed q_d (E/F, kernel settings); the first neighbouring
    converged pair with dE/dq_d from + to - goes to ``continue_branches_ef`` on q_d.
    Returns a dict with ``status`` (the kernel status, or ``branch_unresolved`` when no
    jump yields a bracket), per-jump scans, and ``branch`` (the BranchResult or None)."""
    from types import SimpleNamespace

    from prrs.reliability import GateRejected

    from .afir import DisplacementCoordinate

    cfg = ef_settings
    out = {"status": "branch_unresolved", "jumps": [], "branch": None}
    sides = list(branch.sides)
    for k, side in enumerate(sides):
        if side.get("termination") != "jump":
            continue
        last = side["points"][-1]
        before = np.asarray(last["positions"], float)
        # the other branch at the same q: the opposite side's accepted point nearest in q
        # (in a hysteresis window both branches exist there); the kernel's post-jump frame
        # only when the opposite side has none (near a fold the "jump" frame can still be
        # on this branch: fold bench 2026-11-15, d = 0.004 A)
        other = [p for j, sd in enumerate(sides) if j != k for p in sd.get("points", [])]
        if other:
            near = min(other, key=lambda p: abs(p["q_actual"] - last["q_actual"]))
            after, after_from = np.asarray(near["positions"], float), "other_branch"
        else:
            after, after_from = np.asarray(side["jump"]["positions_after"], float), "jump_frame"
        rec = {"origin": side.get("origin"), "status": "no_bracket", "scan": [],
               "after_from": after_from}  # fmt: skip
        out["jumps"].append(rec)
        try:
            coord = DisplacementCoordinate.build(before, after)
        except ValueError:
            rec["status"] = "no_displacement"
            continue
        rec["length_A"] = coord.length_A
        after_aligned = _align(before, after)
        pts = []
        with ledger.segment(task=task, stage="rebracket_qd", attempt=attempt) as seg:
            for s in np.linspace(0.0, 1.0, n_points):
                q = s * coord.length_A
                start = before + s * (after_aligned - before)  # q restored by the constraint
                try:
                    r = _relax(template, coord, start, q, "qd", factory, config, cfg, [],
                               ceiling_eV)  # fmt: skip
                except _Stop as exc:
                    rec["scan"].append({"q": q, "status": str(exc)})
                    continue
                except GateRejected as exc:  # a point, not the search (as the kernel)
                    rec["scan"].append({"q": q, "status": f"gate_rejected:{exc.code}"})
                    continue
                rec["scan"].append({"q": q, "status": "converged",
                                    "energy_eV": r["energy_eV"],
                                    "reaction_force": r["reaction_force"]})  # fmt: skip
                pts.append(SimpleNamespace(q_target=r["q_actual"], positions=r["positions"],
                                           reaction_force=r["reaction_force"]))  # fmt: skip
            seg.outcome = f"{len(pts)} of {n_points} converged"
        pair = next(
            ((a, b) for a, b in zip(pts, pts[1:])
             if a.reaction_force > cfg.reaction_force_tol
             and b.reaction_force < -cfg.reaction_force_tol),
            None,
        )  # fmt: skip
        if pair is None:
            continue
        rec["status"] = "bracketed"
        rec["bracket"] = [pair[0].q_target, pair[1].q_target]
        br = continue_branches_ef(template.copy(), coord, list(pair), factory, config, ledger,
                                  task, branch_settings=branch_settings, settings=cfg,
                                  ceiling_eV=ceiling_eV, attempt=attempt)  # fmt: skip
        rec["kernel_status"] = br.status
        out["status"], out["branch"] = br.status, br
        if br.status == "converged_candidate":
            break
    return out


__all__ = [
    "BranchSettings",
    "EFSettings",
    "link",
    "path_minima",
    "rebracket_jump",
    "relax_min",
    "ts_from_path",
]
