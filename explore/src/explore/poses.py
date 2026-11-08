"""L0 for discovery v1 (DISCOVERY_LAYER_DESIGN_2026-10-13 §2): source conformers on a
torsion grid, and rigid placements of two fragments. Pure geometry except
`qualify_conformers`, which spends model forwards through the ledger (shared cost P)."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, replace

import numpy as np

from .streams import stream

RULE_VERSION = "L0-v1"


@dataclass(frozen=True)
class Fragment:
    key: str
    atoms: tuple[int, ...]
    ordinal: int = 0


@dataclass
class Conformer:
    conformer_id: str
    fragment_key: str
    atoms: tuple[int, ...]
    positions: np.ndarray
    grid: tuple[int, ...]
    status: str


def _kabsch_onto(moving, target):
    """Rigidly superpose `moving` onto `target` (same atom order); exact for n <= 2."""
    a = moving - moving.mean(axis=0)
    b = target - target.mean(axis=0)
    if len(a) < 3:
        return (
            a + target.mean(axis=0) if len(a) == 1 else _align_line(a, b) + target.mean(axis=0)
        )
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt)) or 1.0
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return a @ rot + target.mean(axis=0)


def _align_line(a, b):
    # a diatomic: rotate its axis onto the target axis; the axial angle is undefined
    ua, ub = a[1] - a[0], b[1] - b[0]
    ua, ub = ua / np.linalg.norm(ua), ub / np.linalg.norm(ub)
    from scipy.spatial.transform import Rotation

    rot, _ = Rotation.align_vectors([ub], [ua])
    return rot.apply(a)


def fragments(atoms, config):
    from prrs.chemistry import _formula, _fragments, canonical_labels
    from prrs.state import encode

    graph = encode(atoms, config.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    out = []
    for group in _fragments(graph.edges, graph.index):
        colours = ",".join(sorted(str(labels[i]) for i in group))
        out.append(
            Fragment(f"{_formula(atoms.numbers, group)}|{colours}", tuple(sorted(group)))
        )
    # identical fragments tie on key: break by centroid, which relabelling cannot move
    ordered = sorted(
        out,
        key=lambda f: (f.key, tuple(np.round(atoms.positions[list(f.atoms)].mean(axis=0), 8))),
    )
    return [replace(f, ordinal=i) for i, f in enumerate(ordered)]


def _fragment_atoms(atoms, fragment):
    return atoms[list(fragment.atoms)]


def _ordered_torsions(sub, config):
    """Every genuine torsion of the fragment (identity orbits: two equivalent methoxy
    groups are two coordinates), in an order relabelling cannot change."""
    from prrs.chemistry import canonical_labels
    from prrs.perturbations import rotatable_torsions
    from prrs.state import encode

    graph = encode(sub, config.bond_scale)
    labels = canonical_labels(sub.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels, {int(i): int(i) for i in graph.index})

    def key(t):
        j, k = t.indices[1], t.indices[2]
        mid = (sub.positions[j] + sub.positions[k]) / 2
        return (tuple(sorted((str(labels[j]), str(labels[k])))), tuple(np.round(mid, 8)))

    return sorted(genuine, key=key), rotors, graph


def _clashes(positions, numbers, config):
    from ase.data import covalent_radii

    r = covalent_radii[numbers]
    d = np.linalg.norm(positions[:, None] - positions[None], axis=-1)
    np.fill_diagonal(d, np.inf)
    return bool(np.any(d < config.min_distance_ratio * (r[:, None] + r[None])))


def conformer_grid(atoms, fragment, config, steps_per_turn, cap, seed):
    """Torsion-grid conformers of one fragment. Grid-index sign may differ under relabelling;
    the set of geometries is what is invariant."""
    from prrs.perturbations import step_torsion
    from prrs.state import encode

    sub = _fragment_atoms(atoms, fragment)
    torsions, rotors, graph = _ordered_torsions(sub, config)
    axes = [range(max(1, steps_per_turn // t.symmetry_order)) for t in torsions]
    grid = list(itertools.product(*axes))
    origin = tuple(0 for _ in torsions)
    if len(grid) > cap:
        rest = [g for g in grid if g != origin]
        pick = stream(seed, "conformer_cap", f"f{fragment.ordinal}|{fragment.key}").permutation(
            len(rest)
        )[: cap - 1]
        chosen = [origin] + [rest[i] for i in sorted(pick)]
    else:
        chosen = grid
    confs, counts = [], {}
    for g in chosen:
        work = sub.copy()
        for t, idx in zip(torsions, g):
            n = max(1, steps_per_turn // t.symmetry_order)
            step_torsion(work, t, idx * (2 * math.pi / t.symmetry_order) / n)
        status = "generated"
        if _clashes(work.positions, work.numbers, config):
            status = "overlap"
        elif encode(work, config.bond_scale).edges != graph.edges:
            status = "graph_changed"
        counts[status] = counts.get(status, 0) + 1
        confs.append(
            Conformer(
                conformer_id=f"f{fragment.ordinal}|{fragment.key}:t{'.'.join(map(str, g))}:{RULE_VERSION}",
                fragment_key=fragment.key,
                atoms=fragment.atoms,
                positions=_kabsch_onto(work.positions, sub.positions),
                grid=g,
                status=status,
            )
        )
    report = {
        "genuine_torsions": len(torsions),
        "rotors_not_enumerated": len(rotors),
        "rings_not_enumerated": True,
        "grid_size": len(grid),
        "censored": len(grid) - len(chosen),
        "status_counts": counts,
    }
    return confs, report


def _jsonable(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [_jsonable(v) for v in x]
    if x is None or isinstance(x, (bool, int, float, str)):
        return x
    return str(x)


def qualify_conformers(
    atoms, conformers, factory, config, ledger, seed, fragment_states, attempt=0
):
    """Quench each generated conformer, confirm the minimum, and merge strict duplicates
    with the PRRS basin rule (same bond graph, `basin_energy_eV`, `basin_rmsd_A` minimised
    over the graph automorphisms). Every forward is shared preparation cost P (stage
    "prepare"); failures keep their status and are never dropped.

    `fragment_states` maps `Conformer.fragment_key` to `{"charge": int, "multiplicity":
    int}`; it is stamped on each fragment as `info["charge"]` / `info["spin"]` (the keys
    PRRS `_stamped_with_charge` uses). Nothing is inferred: a missing entry is the status
    "electronic_state_missing" and costs no forward.

    `quench_failed` / `not_minimum` rows keep the geometry reached and the diagnostics of
    `quench` / `confirm_minimum` under `row["diagnostics"]`. An unexpected exception for one
    conformer is the status `error:<ExceptionName>`; the loop goes on. `attempt` is passed to
    the ledger segments."""
    from prrs.chemistry import automorphisms, symmetric_rmsd
    from prrs.io import snapshot
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.runner import confirm_minimum, quench
    from prrs.state import encode

    from .ledger import BudgetExhausted

    kept, rows = [], []
    for conf in conformers:
        row = {
            "conformer_id": conf.conformer_id,
            "status": None,
            "positions": None,
            "energy_eV": None,
            "duplicate_of": None,
            "diagnostics": None,
        }
        rows.append(row)
        if conf.status != "generated":
            row["status"] = f"skipped:{conf.status}"
            continue
        state = fragment_states.get(conf.fragment_key)
        if state is None:
            row["status"] = "electronic_state_missing"
            continue
        sub = atoms[list(conf.atoms)].copy()
        sub.set_constraint()
        sub.info["charge"] = int(state["charge"])
        sub.info["spin"] = int(state["multiplicity"])
        sub.positions = conf.positions.copy()
        diag = {"quench": None, "confirm": None}
        try:
            guard = GuardedCalculator(factory(), config)
            sub.calc = guard
            with ledger.segment(
                task=conf.conformer_id, stage="prepare", attempt=attempt
            ) as seg:
                converged, qinfo = quench(sub, guard, config)
                diag["quench"] = _jsonable(qinfo)
                seg.outcome = "converged" if converged else "quench_failed"
                endpoint = snapshot(sub, guard) if converged else None
            if endpoint is None:
                row.update(
                    status="quench_failed",
                    positions=sub.positions.copy(),
                    diagnostics=diag,
                )
                continue
            with ledger.segment(
                task=conf.conformer_id, stage="prepare", attempt=attempt
            ) as seg:
                confirmed, info = confirm_minimum(endpoint, factory, config, seed)
                diag["confirm"] = _jsonable(info)
                reason = str(info.get("reason", ""))
                gated = not confirmed and reason.startswith("gate_")
                seg.outcome = (
                    "minimum" if confirmed else ("gate_rejected" if gated else "not_minimum")
                )
            if gated:
                row["status"] = f"gate_rejected:{reason[len('gate_') :]}"
                continue
            if not confirmed:
                row.update(
                    status="not_minimum",
                    positions=endpoint.positions.copy(),
                    diagnostics=diag,
                )
                continue
        except BudgetExhausted:
            row["status"] = "budget_exhausted"
            continue
        except GateRejected as exc:
            row["status"] = f"gate_rejected:{exc.code}"
            continue
        except Exception as exc:
            row["status"] = f"error:{type(exc).__name__}"
            continue
        energy = float(endpoint.get_potential_energy())
        row.update(positions=endpoint.positions.copy(), energy_eV=energy)
        edges = encode(endpoint, config.bond_scale).edges
        group, _ = automorphisms(endpoint, config.bond_scale, None, config.automorphism_limit)
        for other_id, other, other_edges in kept:
            if (
                edges == other_edges
                and abs(energy - other.get_potential_energy()) <= config.basin_energy_eV
                and symmetric_rmsd(endpoint, other, group) <= config.basin_rmsd_A
            ):
                row.update(status="duplicate", duplicate_of=other_id)
                break
        else:
            row["status"] = "qualified"
            kept.append((conf.conformer_id, endpoint, edges))
    return rows


def conformers_from_qualification(conformers, rows):
    """Conformers for the rows with status "qualified" (matched by conformer_id), carrying the
    relaxed positions; poses are built from these (ruling R9)."""
    by_id = {r["conformer_id"]: r for r in rows}
    out = []
    for c in conformers:
        row = by_id.get(c.conformer_id)
        if row is None or row["status"] != "qualified":
            continue
        out.append(
            Conformer(
                conformer_id=c.conformer_id,
                fragment_key=c.fragment_key,
                atoms=c.atoms,
                positions=np.asarray(row["positions"], float).copy(),
                grid=c.grid,
                status="qualified",
            )
        )
    return out


@dataclass
class Pose:
    pose_id: str
    conformers: tuple[str, ...]
    sample: int | None
    shell: float | None
    direction: np.ndarray | None
    quaternion: np.ndarray | None
    positions: np.ndarray
    min_contact_ratio: float | None
    status: str
    t_A: float | None = None


def _radii(numbers):
    from ase.data import covalent_radii

    return covalent_radii[np.asarray(numbers)]


def contact_ratio(pos_a, num_a, pos_b, num_b):
    d = np.linalg.norm(pos_a[:, None] - pos_b[None], axis=-1)
    return float(np.min(d / (_radii(num_a)[:, None] + _radii(num_b)[None])))


_MARCH_CAP, _REFINE_CAP, _BISECT_CAP = 10_000, 1_000, 200


def place(pos_a, num_a, pos_b, num_b, direction, quaternion, shell):
    """Put B's centroid on the ray from A's centroid along `direction`, B rotated about
    its centroid by `quaternion`, at the OUTERMOST distance where the closest-contact
    ratio equals `shell`. March inward with the Lipschitz-safe step
    (g - shell) * min(r_i + r_j), which can never jump over a crossing, then bisect."""
    from scipy.spatial.transform import Rotation

    direction = _unit(np.asarray(direction, float))
    centre_a = pos_a.mean(axis=0)
    body = Rotation.from_quat(quaternion).apply(pos_b - pos_b.mean(axis=0))
    rmin = float(np.min(_radii(num_a)[:, None] + _radii(num_b)[None]))
    reach = np.max(np.linalg.norm(pos_a - centre_a, axis=1)) + np.max(
        np.linalg.norm(body, axis=1)
    )

    def g(t):
        return contact_ratio(pos_a, num_a, body + centre_a + t * direction, num_b)

    rmax = float(np.max(_radii(num_a)[:, None] + _radii(num_b)[None]))
    # every pair distance >= t - reach, so g(outer) >= (shell * rmax + 1) / rmax > shell
    outer = reach + shell * rmax + 1.0

    def result(t, status):
        t = t if np.isfinite(t) else outer
        return body + centre_a + t * direction, g(t), status, float(t)

    t = outer
    for _ in range(_MARCH_CAP):
        gt = g(t)
        if gt <= shell:
            break
        step = (gt - shell) * rmin
        if step < 1e-9:
            break
        t -= step
        if not np.isfinite(t):
            return result(outer, "placement_failed")
    else:
        return result(t, "placement_failed")
    lo, hi = t, outer
    for _ in range(_REFINE_CAP):  # the march stops at or just above the crossing
        if g(lo) <= shell:
            break
        lo -= 1e-9
    else:
        return result(lo, "placement_failed")
    for _ in range(_BISECT_CAP):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if g(mid) <= shell else (lo, mid)
        if hi - lo < 1e-10:
            break
    else:
        return result(hi, "placement_failed")
    return result(hi, "placed")


def _unit(v):
    return v / np.linalg.norm(v)


def sample_poses(atoms, frags, conf_a, conf_b, shells, n_samples, seed):
    """Uniform approach directions (normalised 3-D Gaussian) and uniform orientations
    (normalised 4-D Gaussian = Haar on SO(3)), one keyed stream per sample index, so
    the k-th sample is the same however many are drawn. Each sample is placed at every
    shell."""
    if len(frags) != 2:
        raise ValueError(
            f"v1 two-fragment poses need exactly 2 fragments, got {len(frags)}; "
            "use single_fragment_poses for 1"
        )
    fa, fb = frags
    if conf_a.atoms != fa.atoms or conf_b.atoms != fb.atoms:
        raise ValueError(
            "conformer atoms do not match the fragments (order: frags[0], frags[1])"
        )
    bad = next(
        (c.status for c in (conf_a, conf_b) if c.status not in ("generated", "qualified")), None
    )
    num_a, num_b = atoms.numbers[list(fa.atoms)], atoms.numbers[list(fb.atoms)]
    pair_key = f"{conf_a.conformer_id}+{conf_b.conformer_id}"
    poses = []
    for k in range(n_samples):
        rng = stream(seed, "pose", pair_key, k)
        direction = _unit(rng.normal(size=3))
        quaternion = _unit(rng.normal(size=4))
        for shell in shells:
            new_b, ratio, status, t_a = place(
                conf_a.positions, num_a, conf_b.positions, num_b, direction, quaternion, shell
            )
            if bad is not None:
                status = f"conformer_{bad}"
            positions = atoms.positions.copy()
            positions[list(fa.atoms)] = conf_a.positions
            positions[list(fb.atoms)] = new_b
            poses.append(
                Pose(
                    pose_id=f"{pair_key}:k{k}:f{shell:g}:{RULE_VERSION}",
                    conformers=(conf_a.conformer_id, conf_b.conformer_id),
                    sample=k,
                    shell=float(shell),
                    direction=direction,
                    quaternion=quaternion,
                    positions=positions,
                    min_contact_ratio=ratio,
                    status=status,
                    t_A=t_a,
                )
            )
    return poses


def single_fragment_poses(atoms, conformers):
    poses = []
    for c in conformers:
        positions = atoms.positions.copy()
        positions[list(c.atoms)] = c.positions
        poses.append(
            Pose(
                c.conformer_id,
                (c.conformer_id,),
                None,
                None,
                None,
                None,
                positions,
                None,
                c.status,
            )
        )
    return poses


def write_l0(path, poses, electronic_state, generation, fragment_atoms):
    import json
    from pathlib import Path

    with Path(path).open("w", encoding="utf-8") as fh:
        for p in poses:
            fh.write(
                json.dumps(
                    {
                        "source_microstate": None,  # set when the pose is qualified as a basin
                        "conformer": list(p.conformers),
                        "pose": p.pose_id,
                        "atoms_map": {"fragments": [list(f) for f in fragment_atoms]},
                        "electronic_state": electronic_state,
                        "generation_rule": {**generation, "shell": p.shell, "sample": p.sample},
                        "seed": generation.get("seed"),
                        "preparation": "raw_placement",
                        "qualification": p.status,
                        "cost": {"forwards": 0},
                        "direction": None if p.direction is None else p.direction.tolist(),
                        "quaternion": None if p.quaternion is None else p.quaternion.tolist(),
                        "min_contact_ratio": p.min_contact_ratio,
                        "t_A": p.t_A,
                        "positions": p.positions.tolist(),
                    }
                )
                + "\n"
            )
