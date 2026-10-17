"""prepare / shoot / report for one maturation batch (plan M2-M4).

prepare (M2)  PARS decides WHICH paths (its refinement order, one per reaction
              tube first); the candidate geometries come from the source PRRS
              trajectory itself, free segment only. PARS's representative files
              are resampled chains -- interpolated, rigid-aligned -- and a
              candidate taken from them would pass interpolation off as MD.
              ``CandidateRule.proposal_ids`` keeps only the listed candidates;
              with ``from_batch`` they are copied from an earlier batch, so a
              confirm or sidecheck batch shoots what explore diagnosed (M5).
shoot (M3)    One ``prrs.runner.shoot(shots=1)`` per task: thermal momenta from
              the task's seed, no probe delivered, first-entry termination on the
              region's two graphs. Judging reads the source registry
              (``Registry.identify``) and never writes to it -- the same reading
              ``search.py`` applies to its own committor shots.
report (M4)   Rebuilt from commit markers every time, so a re-report never
              counts a shot twice.

Statistics vocabulary (plan M3): only A and B enter p_B_given_AB. Everything
else -- unabsorbed, third basin, gate, failure -- keeps its own label, and a
candidate whose excluded share exceeds ``other_fraction_max`` is unresolved.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import time
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from ase.io import read as ase_read
from ase.io import write as ase_write

from .io import (
    MANIFEST_NAME,
    PROPOSALS_NAME,
    OutputConflict,
    OutputStore,
    read_json_strict,
    read_jsonl,
    sha256_file,
    verify_file_hashes,
    write_json_atomic,
)
from .schema import (
    MicrostateRef,
    PathProposal,
    SchemaError,
    ShotRecord,
    ShootingEstimate,
    TransitionDiagnosis,
    sha256_object,
    utc_now,
    wilson_interval,
)

PREPARE_NAME = "prepare.json"
REPORT_NAME = "report.json"
HANDOFF_NAME = "refinement_handoff.json"
FREE_PHASES = ("free_start", "free")


# ---- the source run ---------------------------------------------------------


def prrs_implementation_sha256() -> str:
    """The fingerprint ``search.py`` writes into a run manifest: name + bytes of
    every ``src/prrs/*.py``, sorted."""
    import prrs

    digest = hashlib.sha256()
    for path in sorted(Path(prrs.__file__).parent.glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def source_config(run_dir):
    from prrs.config import SearchConfig

    manifest = read_json_strict(Path(run_dir) / "config.json")
    return SearchConfig.from_dict(manifest["config"])


def _network(run_dir) -> dict:
    return read_json_strict(Path(run_dir) / "network.json")


def _edges_of(atoms, config) -> list:
    from prrs.state import encode

    return sorted(
        (int(i), int(j))
        for i, j in encode(atoms, config.bond_scale, active=config.active_atoms).edges
    )


def microstate_ref(run_dir, microstate_id) -> MicrostateRef:
    """A region side from the source run's own published microstate."""
    run_dir = Path(run_dir)
    config = source_config(run_dir)
    for node in _network(run_dir)["chemical_nodes"]:
        for state in node["microstates"]:
            if state["id"] == microstate_id:
                path = run_dir / state["structure"]
                return MicrostateRef(
                    node_id=node["id"],
                    microstate_id=microstate_id,
                    chemical_key=node["chemical_key"],
                    structure_path=state["structure"],
                    structure_sha256=sha256_file(path),
                    edges=_edges_of(ase_read(str(path)), config),
                )
    raise SchemaError(f"microstate {microstate_id!r} is not in {run_dir}/network.json")


def verify_region(protocol) -> None:
    """Each side's structure file is the one hashed, and encodes to its graph.

    First-entry termination compares the shot's graph with these edges exactly,
    so edges that do not come from the structure under the run's own bond rule
    would make a side unreachable by construction.
    """
    run_dir = Path(protocol.source.run_dir)
    config = source_config(run_dir)
    for side, ref in sorted(protocol.regions.sides.items()):
        path = run_dir / ref.structure_path
        if not path.is_file():
            raise OutputConflict(f"region side {side}: {path} is missing")
        if sha256_file(path) != ref.structure_sha256:
            raise OutputConflict(f"region side {side}: {path} changed since the protocol")
        if [list(e) for e in _edges_of(ase_read(str(path)), config)] != [
            list(e) for e in ref.edges
        ]:
            raise OutputConflict(
                f"region side {side}: the structure does not encode to the recorded "
                "edges under the source run's bond rule"
            )


# ---- M2: candidate extraction -----------------------------------------------


def _kabsch_onto(moving, reference):
    a = moving - moving.mean(axis=0)
    b = reference - reference.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    return a @ rotation + reference.mean(axis=0)


def _arc_lengths(frames, align) -> np.ndarray:
    positions = [f.positions for f in frames]
    if align:
        positions = [_kabsch_onto(p, positions[0]) for p in positions]
    steps = [float(np.linalg.norm(b - a)) for a, b in zip(positions, positions[1:])]
    return np.concatenate([[0.0], np.cumsum(steps)])


def _region_pair(protocol) -> set:
    return {protocol.regions.sides["A"].node_id, protocol.regions.sides["B"].node_id}


def _choose_paths(channel, rule) -> tuple:
    """Up to max_paths in PARS refinement order, distinct tubes first (plan M2)."""
    tube_of = {m: t["tube_id"] for t in channel["tubes"] for m in t["member_ids"]}
    order = list(channel["refinement_order"])
    chosen, tubes_used = [], set()
    for cid in order:
        if len(chosen) < rule.max_paths_per_channel and tube_of.get(cid) not in tubes_used:
            chosen.append(cid)
            tubes_used.add(tube_of.get(cid))
    for cid in order:
        if len(chosen) < rule.max_paths_per_channel and cid not in chosen:
            chosen.append(cid)
    not_sampled = sorted(
        t["tube_id"] for t in channel["tubes"] if t["tube_id"] not in tubes_used
    )
    return chosen, not_sampled


def extract_proposals(protocol, out_dir, write=True) -> tuple:
    """(proposals, ledger). Geometry files are written under ``out_dir`` when
    ``write``; the ledger names every channel, tube and position not proposed."""
    source = protocol.source
    run_dir = Path(source.run_dir)
    if source.pars_selection_dir is None:
        raise SchemaError("prepare needs SourceRef.pars_selection_dir (a PARS select output)")
    selection_dir = Path(source.pars_selection_dir)
    selection_path = selection_dir / "pars_selection.json"
    if (
        source.pars_selection_sha256 is not None
        and sha256_file(selection_path) != source.pars_selection_sha256
    ):
        raise OutputConflict(f"{selection_path} changed since the protocol was written")
    selection = read_json_strict(selection_path)
    candidates = {}
    for line in (selection_dir / "candidates.jsonl").read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            candidates[record["candidate_id"]] = record
    if Path(selection["run_dir"]).resolve() != run_dir.resolve():
        raise OutputConflict(
            f"the PARS selection was made on {selection['run_dir']}, not {run_dir}"
        )
    rule = protocol.candidates
    wanted = set(rule.proposal_ids) if rule.proposal_ids is not None else None
    reference = ase_read(str(run_dir / protocol.regions.sides["A"].structure_path))
    pair = _region_pair(protocol)
    graphs = protocol.regions.graph_sha256s
    proposals, ledger = [], {"channels": {}, "meaning": _LEDGER_MEANING}
    for channel_id, channel in sorted(selection["channels"].items()):
        members = [candidates[m] for m in sorted({*channel["refinement_order"]})]
        ends = {
            frozenset(
                (m["provenance"]["reactant_state_id"], m["provenance"]["product_state_id"])
            )
            for m in members
        }
        if ends != {frozenset(pair)}:
            ledger["channels"][channel_id] = {"status": "not_this_region"}
            continue
        chosen, not_sampled = _choose_paths(channel, rule)
        entry = {
            "status": "sampled",
            "paths": [],
            "not_sampled_budget": not_sampled,
        }
        for cid in chosen:
            provenance = candidates[cid]["provenance"]
            attempt = provenance["attempt_id"]
            trajectory = run_dir / "trials" / attempt / "response.traj"
            frames = ase_read(str(trajectory), index=":")
            free = [i for i, f in enumerate(frames) if f.info.get("phase") in FREE_PHASES]
            path_entry = {"candidate_id": cid, "attempt_id": attempt, "positions": []}
            entry["paths"].append(path_entry)
            if len(free) < 3 or free != list(range(free[0], free[-1] + 1)):
                path_entry["skipped"] = "free segment shorter than three contiguous frames"
                continue
            if any(not np.array_equal(frames[i].numbers, reference.numbers) for i in free):
                raise OutputConflict(
                    f"{trajectory}: atoms differ from region side A; all atoms, every "
                    "hydrogen included, must be kept (plan M2)"
                )
            segment = [frames[i] for i in free]
            arc = _arc_lengths(segment, align=not any(reference.pbc))
            if arc[-1] <= 0.0:
                path_entry["skipped"] = "free segment does not move"
                continue
            parent_sha = sha256_file(trajectory)
            seen = set()
            for k, fraction in enumerate(rule.arc_fractions):
                target = fraction * arc[-1]
                j = int(np.searchsorted(arc, target, side="right")) - 1
                j = min(max(j, 0), len(segment) - 2)
                w = (target - arc[j]) / (arc[j + 1] - arc[j]) if arc[j + 1] > arc[j] else 0.0
                interpolate = rule.allow_interpolation and 1e-6 < w < 1 - 1e-6
                if interpolate:
                    positions = (1 - w) * segment[j].positions + w * segment[j + 1].positions
                    key = ("interp", j, round(w, 6))
                else:
                    near = j if w < 0.5 else j + 1
                    positions = segment[near].positions.copy()
                    key = ("frame", near)
                if key in seen:
                    path_entry["positions"].append(
                        {"arc_fraction": fraction, "skipped": "same frame as an earlier one"}
                    )
                    continue
                seen.add(key)
                proposal_id = f"mp-{channel_id}-{attempt}-k{k}"
                if wanted is not None and proposal_id not in wanted:
                    path_entry["positions"].append(
                        {
                            "arc_fraction": fraction,
                            "skipped": "not in CandidateRule.proposal_ids",
                        }
                    )
                    continue
                geometry_file = f"geometries/{proposal_id}.extxyz"
                atoms = reference.copy()
                atoms.positions = positions
                atoms.info = {"proposal_id": proposal_id, "arc_fraction": fraction}
                geometry_sha = None
                if write:
                    path = Path(out_dir) / geometry_file
                    path.parent.mkdir(parents=True, exist_ok=True)
                    ase_write(str(path), atoms)
                    geometry_sha = sha256_file(path)
                common = dict(
                    proposal_id=proposal_id,
                    batch_role=protocol.batch_role,
                    channel_id=channel_id,
                    event_key=provenance.get("event_key"),
                    tube_id=candidates[cid].get("tube_id"),
                    source_candidate_id=cid,
                    source_run_dir=str(run_dir),
                    source_attempt_id=attempt,
                    parent_frame_start=free[0],
                    parent_frame_stop=free[-1] + 1,
                    parent_path_sha256=parent_sha,
                    arc_fraction=fraction,
                    geometry_file=geometry_file,
                    geometry_sha256=geometry_sha or sha256_object(positions.tolist()),
                    n_atoms=len(atoms),
                    region_graph_sha256s=dict(graphs),
                )
                if interpolate:
                    proposal = PathProposal(
                        **common,
                        proposal_kind="path_interpolation",
                        interpolation={
                            "method": "linear",
                            "weights": [1 - w, w],
                            "frames": [free[j], free[j + 1]],
                        },
                    )
                else:
                    proposal = PathProposal(
                        **common, proposal_kind="observed_frame", frame_index=free[key[1]]
                    )
                proposals.append(proposal)
                path_entry["positions"].append(
                    {"arc_fraction": fraction, "proposal_id": proposal_id}
                )
        ledger["channels"][channel_id] = entry
    if wanted is not None:
        produced = {p.proposal_id for p in proposals}
        missing = [i for i in rule.proposal_ids if i not in produced]
        if missing:
            raise OutputConflict(
                f"CandidateRule.proposal_ids: {missing} are not proposed from this "
                "selection under this rule (prepare --dry-run lists the ids it proposes)"
            )
    return proposals, ledger


def proposals_from_batch(protocol, out_dir, write=True) -> tuple:
    """(proposals, ledger) for listed candidates an earlier batch already holds.

    The geometry is that batch's file, byte for byte: a confirm or sidecheck
    batch shoots what the explore batch diagnosed (plan M5), and a
    re-extraction would silently differ once the selection or the rule changed.
    """
    rule = protocol.candidates
    batch = Path(rule.from_batch)
    manifest_path = batch / MANIFEST_NAME
    if not manifest_path.is_file():
        raise OutputConflict(f"CandidateRule.from_batch: no {MANIFEST_NAME} in {batch}")
    manifest = read_json_strict(manifest_path)
    if manifest.get("kind") != "mlpath-maturation":
        raise OutputConflict(f"CandidateRule.from_batch: {batch} is not a maturation batch")
    recorded = {d["proposal_id"]: d for d in read_jsonl(batch / PROPOSALS_NAME)}
    missing = [i for i in rule.proposal_ids if i not in recorded]
    if missing:
        raise OutputConflict(f"CandidateRule.proposal_ids: {missing} are not in {batch}")
    origin = {
        "dir": str(batch),
        "protocol_sha256": manifest["protocol_sha256"],
        "batch_role": manifest["batch_role"],
    }
    graphs = protocol.regions.graph_sha256s
    proposals = []
    for proposal_id in rule.proposal_ids:
        data = recorded[proposal_id]
        if data["region_graph_sha256s"] != graphs:
            raise OutputConflict(
                f"{proposal_id} was proposed against another A/B region than this "
                "protocol's; its shots would be judged against other sides"
            )
        geometry = batch / data["geometry_file"]
        verify_file_hashes({geometry: data["geometry_sha256"]})
        if write:
            target = Path(out_dir) / data["geometry_file"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(geometry, target)
        proposals.append(
            PathProposal.from_dict(
                {**data, "batch_role": protocol.batch_role, "origin_batch": origin}
            )
        )
    ledger = {
        "from_batch": origin,
        "proposals": list(rule.proposal_ids),
        "meaning": "candidates listed by the protocol, copied from an earlier batch",
    }
    return proposals, ledger


def _candidates(protocol, out_dir, write):
    if protocol.candidates.from_batch is not None:
        return proposals_from_batch(protocol, out_dir, write)
    return extract_proposals(protocol, out_dir, write)


_LEDGER_MEANING = (
    "not_this_region: the channel does not connect the protocol's A and B nodes; "
    "not_sampled_budget: tubes beyond max_paths_per_channel, not absent or unimportant; "
    "a skipped position names why it was not proposed"
)


def prepare(protocol, out_dir, dry_run=False) -> dict:
    verify_region(protocol)
    fingerprint = prrs_implementation_sha256()
    # a listed id that cannot be produced refuses here, before anything is written
    proposals, ledger = _candidates(protocol, out_dir, write=False)
    if dry_run:
        shots = protocol.shooting.shots_for_role()
        return {
            "dry_run": True,
            "proposals": [p.proposal_id for p in proposals],
            "tasks": len(proposals) * shots,
            "ledger": ledger,
        }
    store = OutputStore.initialize(out_dir, protocol, prrs_implementation_sha256=fingerprint)
    proposals, ledger = _candidates(protocol, out_dir, write=True)
    for proposal in proposals:
        store.add_proposal(proposal)
    tasks = store.build_task_table(protocol, proposals)
    summary = {
        "prepared_utc": utc_now(),
        "prrs_implementation_sha256": fingerprint,
        "proposals": [p.proposal_id for p in proposals],
        "tasks": len(tasks),
        "ledger": ledger,
    }
    write_json_atomic(Path(out_dir) / PREPARE_NAME, summary)
    return summary


# ---- M3: one shot -----------------------------------------------------------


class _Counted:
    """Counts calls into a calculator's ``calculate`` -- the real model
    evaluations. ASE's own property cache sits in front of ``calculate`` and is
    left exactly as it was, so counting does not change what gets computed."""

    def __init__(self, calculator):
        self.calculator = calculator
        self.calls = 0
        inner = calculator.calculate

        def calculate(*args, **kwargs):
            self.calls += 1
            return inner(*args, **kwargs)

        calculator.calculate = calculate


@dataclass
class ShotContext:
    protocol: object
    store: OutputStore
    config: object
    registry: object
    counted: _Counted
    fingerprint: str
    first_entry: dict
    side_states: dict
    probes: dict


def shot_context(protocol, store) -> ShotContext:
    from prrs.calculators import load_factory
    from prrs.network import Registry

    run_dir = Path(protocol.source.run_dir)
    base = source_config(run_dir)
    dynamics = protocol.dynamics
    config = replace(
        base,
        temperature_K=dynamics.temperature_K,
        timestep_fs=dynamics.timestep_fs,
        response_steps=dynamics.free_step_cap,
    )
    sides = protocol.regions.sides
    attempts = {}
    for line in (run_dir / "attempts.jsonl").read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            attempts[record["attempt_id"]] = record
    calculator = load_factory(protocol.source.model_factory, protocol.source.model_kwargs)()
    return ShotContext(
        protocol=protocol,
        store=store,
        config=config,
        registry=Registry.from_network(_network(run_dir), run_dir, base),
        counted=_Counted(calculator),
        fingerprint=prrs_implementation_sha256(),
        first_entry={
            "graphs": {side: {tuple(e) for e in ref.edges} for side, ref in sides.items()},
            "k": dynamics.first_entry_k,
        },
        side_states={side: ref.microstate_id for side, ref in sides.items()},
        probes=attempts,
    )


def _probe_for(context, proposal, seed):
    """The source attempt's probe, carried as a label only: shots run with
    ``deliver=False`` (plan M3: no second perturbation)."""
    from prrs.perturbations import Probe

    data = context.probes[proposal.source_attempt_id]["probe"]
    return Probe(
        family=data["family"],
        indices=tuple(data["indices"]),
        sign=data["sign"],
        amplitude=data["amplitude"],
        seed=seed,
    )


def judge(record, landed, side_states, cap) -> dict:
    """ShotRecord fields from one trial record (plan M3 统计与分类).

    A/B is the first-entry side. A disagreement with the quench is counted on
    the record and never re-judges the shot. Without an entry the shot is
    ``other`` when its quench converged somewhere that is neither side (a third
    basin, identified or not), else ``unabsorbed``.
    """
    status = record["status"]
    failure = record.get("failure") or {}
    entry = record.get("first_entry") or {}
    fields = {
        "landed_microstate": landed,
        "max_nve_drift_eV_atom": record.get("max_nve_drift_eV_atom"),
        "quench_converged": record.get("quench_converged"),
    }
    if status == "rejected":
        return {**fields, "outcome": "gate_rejected", "outcome_reason": failure.get("code")}
    if status == "failed":
        reason = f"{failure.get('code')}: {failure.get('type')}: {failure.get('message')}"
        return {**fields, "outcome": "failed", "outcome_reason": reason}
    side = entry.get("side")
    if side in ("A", "B"):
        return {
            **fields,
            "outcome": side,
            "entry_side": side,
            "entry_step": entry["step"],
            "entry_time_fs": entry.get("time_fs"),
            "entry_quench_consistent": landed == side_states[side],
            "cost_free_steps": entry["step"],
        }
    fields["cost_free_steps"] = cap
    if record.get("quench_converged") and landed not in side_states.values():
        return {
            **fields,
            "outcome": "other",
            "outcome_reason": f"landed:{landed or 'unidentified'}",
        }
    return {
        **fields,
        "outcome": "unabsorbed",
        "outcome_reason": f"no {entry.get('k')}-step entry within {cap} free steps",
    }


def run_shot(context, task, proposal) -> ShotRecord:
    from prrs.runner import shoot

    store = context.store
    shot_dir = store.output_dir / "shots" / task["proposal_id"] / task["shot_id"]
    atoms = ase_read(str(store.output_dir / proposal.geometry_file))
    atoms.info = {}
    atoms.calc = context.counted.calculator
    atoms.get_potential_energy()
    calls_before = context.counted.calls
    started = time.monotonic()
    (outcome,) = shoot(
        atoms,
        _probe_for(context, proposal, task["seed"]),
        lambda: context.counted.calculator,
        context.config,
        shot_dir,
        1,
        task["seed"],
        prefix="shot",
        first_entry=context.first_entry,
    )
    wall = time.monotonic() - started
    record = outcome.record
    write_json_atomic(shot_dir / "trials" / "shot0000" / "result.json", record)
    landed = None
    if outcome.endpoint is not None:
        _, state = context.registry.identify(outcome.endpoint)
        landed = state.id if state is not None else None
    endpoint = shot_dir / "trials" / "shot0000" / "endpoint.extxyz"
    protocol = context.protocol
    return ShotRecord(
        proposal_id=task["proposal_id"],
        shot_id=task["shot_id"],
        batch_role=task["batch_role"],
        seed=task["seed"],
        protocol_sha256=protocol.sha256,
        implementation_sha256=context.fingerprint,
        trial_dir=f"shots/{task['proposal_id']}/{task['shot_id']}/trials/shot0000",
        cost_force_evaluations=context.counted.calls - calls_before,
        cost_wall_s=wall,
        endpoint_sha256=sha256_file(endpoint) if endpoint.exists() else None,
        model_checkpoint_sha256=protocol.source.model_checkpoint_sha256,
        **judge(record, landed, context.side_states, protocol.dynamics.free_step_cap),
    )


# ---- M4: stopping, budget, the shoot loop --------------------------------------


def _excluded(records) -> int:
    return sum(1 for r in records if r.outcome not in ("A", "B"))


def stop_reason(protocol, records) -> str | None:
    """Why a candidate takes no further shots, or None. A pure function of its
    committed records, so a resumed run re-derives the same decision. Explore
    batches only: confirm and sidecheck batches are fixed-sample estimates.

    - excluded_share_already_failed: more shots outside A/B than
      other_fraction_max allows over the full sample -- the candidate is
      unresolved whatever the remaining shots do.
    - interim_side (explore, interim_look set): the first ``interim_look`` shots'
      Wilson interval lies wholly outside the diagnostic band.
    """
    rule = protocol.shooting
    if protocol.batch_role != "explore":
        return None
    full = rule.shots_for_role()
    if _excluded(records) > rule.other_fraction_max * full:
        return "excluded_share_already_failed"
    if rule.interim_look is not None:
        first = {f"s{i:04d}" for i in range(rule.interim_look)}
        look = [r for r in records if r.shot_id in first]
        if len(look) == rule.interim_look:
            b = sum(1 for r in look if r.outcome == "B")
            ab = sum(1 for r in look if r.outcome in ("A", "B"))
            interval = wilson_interval(b, ab)
            low, high = rule.diagnostic_band
            if interval and (interval[1] < low or interval[0] > high):
                return "interim_side"
    return None


def budget_state(protocol, records) -> dict:
    """Usage against every limit, and whether one more shot fits.

    A shot's cost is bounded before it runs only for free steps (the cap); for
    force evaluations and wall time the largest per-shot cost seen so far is the
    reservation, so the last shot cannot overrun by more than the spread of
    per-shot costs.
    """
    budget = protocol.budget
    cap = protocol.dynamics.free_step_cap
    used = {
        "shots": len(records),
        # a shot that died before reporting its steps is charged the cap
        "free_steps": sum(
            cap if r.cost_free_steps is None else r.cost_free_steps for r in records
        ),
        "force_evaluations": sum(r.cost_force_evaluations or 0 for r in records),
        "wall_minutes": sum(r.cost_wall_s or 0.0 for r in records) / 60.0,
    }
    worst_calls = max((r.cost_force_evaluations or 0 for r in records), default=0)
    worst_minutes = max((r.cost_wall_s or 0.0 for r in records), default=0.0) / 60.0
    blocked = []
    if budget.max_shots is not None and used["shots"] + 1 > budget.max_shots:
        blocked.append("max_shots")
    if budget.max_free_steps is not None and used["free_steps"] + cap > budget.max_free_steps:
        blocked.append("max_free_steps")
    if (
        budget.max_force_evaluations is not None
        and used["force_evaluations"] + worst_calls > budget.max_force_evaluations
    ):
        blocked.append("max_force_evaluations")
    if (
        budget.max_wall_minutes is not None
        and used["wall_minutes"] + worst_minutes > budget.max_wall_minutes
    ):
        blocked.append("max_wall_minutes")
    return {"used": used, "limits": budget.to_dict(), "blocked_by": blocked}


def open_batch(out_dir, protocol=None):
    from .schema import MaturationProtocol

    stored = MaturationProtocol.from_dict(read_json_strict(Path(out_dir) / "protocol.json"))
    protocol = protocol or stored
    store = OutputStore.resume(
        out_dir, protocol, prrs_implementation_sha256=prrs_implementation_sha256()
    )
    verify_region(protocol)
    return protocol, store


def shoot_batch(out_dir, max_new_shots=None, dry_run=False, log=print) -> dict:
    protocol, store = open_batch(out_dir)
    proposals = {p.proposal_id: p for p in store.proposals()}
    pending = store.pending_tasks()
    if dry_run:
        return {"dry_run": True, "pending": len(pending), "budget": _budget(protocol, store)}
    context = None
    fired, skipped_stop = 0, 0
    state = "complete"
    for task in pending:
        records = store.committed_records()
        mine = [r for r in records if r.proposal_id == task["proposal_id"]]
        if stop_reason(protocol, mine) is not None:
            skipped_stop += 1
            continue
        if budget_state(protocol, records)["blocked_by"]:
            state = "budget_exhausted"
            break
        if max_new_shots is not None and fired >= max_new_shots:
            state = "paused"
            break
        store.isolate_interrupted(task["proposal_id"], task["shot_id"])
        context = context or shot_context(protocol, store)
        record = run_shot(context, task, proposals[task["proposal_id"]])
        store.commit(record)
        fired += 1
        log(f"{task['task_id']}: {record.outcome}")
    return {
        "state": state,
        "fired": fired,
        "skipped_by_stop_rule": skipped_stop,
        "budget": _budget(protocol, store),
    }


def _budget(protocol, store):
    return budget_state(protocol, store.committed_records())


# ---- M4: report ---------------------------------------------------------------


def diagnose(protocol, proposal, records) -> tuple:
    """(estimate, diagnosis) for one candidate's committed shots (plan M3)."""
    rule = protocol.shooting
    estimate = ShootingEstimate.from_records(records, band=rule.diagnostic_band)
    excluded_share = _excluded(records) / len(records)
    drifts = [r.max_nve_drift_eV_atom for r in records if r.max_nve_drift_eV_atom is not None]
    checks = {
        "entry_quench_mismatches": estimate.entry_quench_mismatches,
        "entry_quench_within_limit": estimate.entry_quench_mismatches
        <= rule.entry_quench_mismatch_max,
        "excluded_share": excluded_share,
        "excluded_share_within_limit": excluded_share <= rule.other_fraction_max,
        "max_nve_drift_eV_atom": max(drifts) if drifts else None,
        "drift_within_tolerance": (
            None
            if rule.drift_tolerance_eV_atom is None
            else bool(drifts) and max(drifts) <= rule.drift_tolerance_eV_atom
        ),
        "shots": len(records),
        "shots_planned": rule.shots_for_role(),
        "stop_reason": stop_reason(protocol, records),
    }
    complete = len(records) == rule.shots_for_role() or checks["stop_reason"] is not None
    if protocol.batch_role == "zero_measure":
        label, basis = "unresolved", "zero_measure_batch"
    elif not complete:
        label, basis = "unresolved", "shots_incomplete"
    elif not checks["entry_quench_within_limit"]:
        label, basis = "unresolved", "entry_quench_mismatch"
    elif not checks["excluded_share_within_limit"]:
        label, basis = "unresolved", "excluded_share_above_limit"
    elif checks["drift_within_tolerance"] is not True:
        label, basis = "unresolved", "drift_not_within_tolerance"
    else:
        label, basis = estimate.band_class, "wilson_interval_vs_band"
    use = "no_action"
    if label == "transition_region_compatible":
        use = (
            "independent_confirmation" if protocol.batch_role == "explore" else "neb_refinement"
        )
    diagnosis = TransitionDiagnosis(
        proposal_id=proposal.proposal_id,
        label=label,
        basis=basis,
        checks=checks,
        recommended_use=use,
        parent_path_sha256=proposal.parent_path_sha256,
    )
    return estimate, diagnosis


def report(out_dir) -> dict:
    protocol, store = open_batch(out_dir)
    records = store.committed_records()
    by_proposal = {}
    for record in records:
        by_proposal.setdefault(record.proposal_id, []).append(record)
    rows, diagnostics, handoff = [], [], []
    for proposal in sorted(store.proposals(), key=lambda p: p.proposal_id):
        mine = by_proposal.get(proposal.proposal_id, [])
        if not mine:
            rows.append({"proposal_id": proposal.proposal_id, "status": "not_shot"})
            continue
        estimate, diagnosis = diagnose(protocol, proposal, mine)
        diagnostics.append(diagnosis.to_dict())
        rows.append(
            {
                "proposal_id": proposal.proposal_id,
                "channel_id": proposal.channel_id,
                "tube_id": proposal.tube_id,
                "arc_fraction": proposal.arc_fraction,
                "proposal_kind": proposal.proposal_kind,
                "estimate": estimate.to_dict(),
                "diagnosis": diagnosis.to_dict(),
            }
        )
        if diagnosis.recommended_use != "no_action":
            handoff.append(
                {
                    "proposal_id": proposal.proposal_id,
                    "recommended_use": diagnosis.recommended_use,
                    "geometry_file": proposal.geometry_file,
                    "geometry_sha256": proposal.geometry_sha256,
                    "proposal_kind": proposal.proposal_kind,
                    "parent": {
                        "run_dir": proposal.source_run_dir,
                        "attempt_id": proposal.source_attempt_id,
                        "candidate_id": proposal.source_candidate_id,
                        "trajectory_sha256": proposal.parent_path_sha256,
                        "free_frames": [
                            proposal.parent_frame_start,
                            proposal.parent_frame_stop,
                        ],
                    },
                    "wilson_95": estimate.wilson_95,
                    "origin_batch": proposal.origin_batch,
                }
            )
    prepare_path = Path(out_dir) / PREPARE_NAME
    summary = {
        "reported_utc": utc_now(),
        "batch_role": protocol.batch_role,
        "protocol_sha256": protocol.sha256,
        "pending_tasks": len(store.pending_tasks()),
        "budget": budget_state(protocol, records),
        "candidates": rows,
        "ledger": read_json_strict(prepare_path)["ledger"] if prepare_path.exists() else None,
        "meaning": (
            "p_B_given_AB is conditional on reaching A or B; transition_region_compatible "
            "means the Wilson interval overlaps the band, not that q is in it and not a "
            "saddle; saddle evidence is the refinement handoff's job (plan E6)"
        ),
    }
    write_json_atomic(Path(out_dir) / REPORT_NAME, summary)
    write_json_atomic(
        Path(out_dir) / HANDOFF_NAME,
        {
            "candidates": handoff,
            "meaning": "geometries to refine or confirm; none of them is a verified saddle",
        },
    )
    _write_jsonl_atomic(Path(out_dir) / "diagnostics.jsonl", diagnostics)
    return summary


def _write_jsonl_atomic(path, rows) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        "".join(json.dumps(r, sort_keys=True, allow_nan=False) + "\n" for r in rows),
        encoding="utf-8",
    )
    tmp.replace(path)


__all__ = [
    "extract_proposals",
    "proposals_from_batch",
    "prepare",
    "shoot_batch",
    "report",
    "judge",
    "stop_reason",
    "budget_state",
    "diagnose",
    "microstate_ref",
    "verify_region",
]
