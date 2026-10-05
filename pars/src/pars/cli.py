"""Command line entry points: ``pars select`` and ``pars refine``.

``select`` consumes a PRRS run directory and writes the channel-conditioned
Pareto selection into an output directory. ``refine`` takes that selection and
refines the representative paths through a backend. The two stages communicate
only through files, so they can run at different times and the selection
record stays auditable.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read as ase_read
from ase.io import write as ase_write

from . import __version__, electronic
from .adapters.prrs_run import PrrsRun
from .config import ParsConfig
from .frontier.archive import ParetoArchive
from .paths.candidate import build_candidates
from .paths.ingest import channel_members_not_ingested, ingest_run, skipped_attempts
from .paths.resample import align_frames, resample_positions
from .refinement.ase_neb import AseNebBackend
from .provenance import (
    calculator_identity,
    condition_groups,
    implementation_sha256,
    numerical_versions,
)
from .schema import PathCandidate, Provenance, RefinedPath


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="pars",
        description=(
            "PARS: Pareto-Adaptive Reaction Path Search. Selects representative "
            "reaction tubes from PRRS candidate paths before expensive refinement."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    select = sub.add_parser(
        "select", help="select representative paths from a PRRS run directory"
    )
    select.add_argument(
        "--run-dir",
        required=True,
        help="directory containing network.json (a <seed>/<arm> dir)",
    )
    select.add_argument("--out", required=True, help="output directory (must not exist)")
    select.add_argument("--config", default=None, help="ParsConfig JSON file")
    select.add_argument(
        "--channels", default=None, help="comma-separated channel ids to process (default: all)"
    )

    refine = sub.add_parser(
        "refine", help="refine representative paths from a selection directory"
    )
    refine.add_argument(
        "--selection", required=True, help="output directory of a previous 'pars select'"
    )
    refine.add_argument(
        "--calculator-factory",
        required=True,
        help="module:attribute resolving to a callable returning an "
        "ASE calculator, e.g. a thin wrapper around a MACE factory",
    )
    refine.add_argument("--only", default=None, help="refine a single candidate id")
    refine.add_argument(
        "--calculator-metadata",
        default=None,
        help="JSON metadata naming model files and calculator settings; otherwise use the factory module's calculator_metadata() hook",
    )
    refine.add_argument(
        "--redo",
        action="store_true",
        help="refine again and replace existing results, even ones recorded "
        "under the same conditions (default: resume, skipping those)",
    )
    refine.add_argument("--images", type=int, default=None)
    refine.add_argument("--fmax-neb", type=float, default=None)
    refine.add_argument(
        "--uncertainty-refusal",
        type=float,
        default=None,
        help="eV/A; exceeding it on any image reports UNRESOLVED_MLP",
    )

    args = parser.parse_args(argv)
    if args.command == "select":
        return _select(args)
    return _refine(args)


# -- select -------------------------------------------------------------------


def _select(args) -> int:
    out = Path(args.out)
    if out.exists():
        raise FileExistsError(f"Will not overwrite {out}")
    config = ParsConfig.from_dict(_load_json(args.config)) if args.config else ParsConfig()
    run = PrrsRun(args.run_dir)
    traces = ingest_run(run, config)
    if args.channels:
        # filter before extraction, so a trace outside the requested channels
        # (even one no channel can place) cannot fail the selection
        wanted = {c.strip() for c in args.channels.split(",") if c.strip()}
        traces = [t for t in traces if t.provenance.channel_id in wanted]
    candidates = build_candidates(traces, config)
    if not candidates:
        raise ValueError(
            "no candidates after ingestion; skipped attempts by status/classification: "
            f"{skipped_attempts(run, config)} (only completed trials with "
            f"classification in {list(config.include_classifications)} become paths)"
        )

    grouped: dict[str, list] = {}
    for candidate in candidates:
        grouped.setdefault(candidate.channel_id, []).append(candidate)
    archives = {
        channel_id: ParetoArchive(channel_id, members, config)
        for channel_id, members in sorted(grouped.items())
    }

    selection = {
        "pars_version": __version__,
        "config": config.to_dict(),
        "run_dir": str(Path(args.run_dir)),
        "units": run.units,
        "counts": {
            "attempts_total": len(run.attempts()),
            "traces": len(traces),
            "candidates": len(candidates),
            "channels": len(archives),
            "skipped_attempts": skipped_attempts(run, config),
            "channel_members_not_ingested": channel_members_not_ingested(run, config),
        },
        "channels": {
            channel_id: archive.selection_dict() for channel_id, archive in archives.items()
        },
    }
    _write_json(out / "pars_selection.json", selection)
    _write_jsonl(
        out / "candidates.jsonl",
        [c.to_summary() for c in sorted(candidates, key=lambda c: c.candidate_id)],
    )
    representatives_dir = out / "representatives"
    for archive in archives.values():
        for candidate_id in archive.representatives:
            _write_representative(
                representatives_dir / f"rep-{candidate_id}.extxyz",
                archive.candidates[candidate_id],
                config,
            )

    total_reps = sum(len(a.representatives) for a in archives.values())
    print(
        f"pars select: {len(candidates)} candidates in {len(archives)} channels; "
        f"{total_reps} representative paths written to {out}"
    )
    for channel_id, archive in sorted(archives.items()):
        print(
            f"  {channel_id}: {archive.selection_dict()['n_candidates']} candidates, "
            f"{len(archive.frontier())} on the frontier, "
            f"{len(archive.tubes)} tubes, "
            f"{len(archive.representatives)} representatives"
        )
    return 0


def _write_representative(path: Path, candidate: PathCandidate, config: ParsConfig) -> None:
    frames = candidate.positions
    if config.refinement_remove_rigid_motion:
        frames = align_frames(frames)
    chain = resample_positions(frames, config.refinement_images)
    frames = []
    for index, positions in enumerate(chain):
        atoms = Atoms(numbers=candidate.numbers, positions=positions)
        atoms.info = {
            "purpose": "pars_refinement_seed",
            "candidate_id": candidate.candidate_id,
            "channel_id": candidate.channel_id,
            "image": index,
            "pars_provenance": json.dumps(candidate.provenance.to_dict()),
        }
        electronic.stamp(atoms.info, _seed_state(candidate.provenance, candidate.candidate_id))
        frames.append(atoms)
    path.parent.mkdir(parents=True, exist_ok=True)
    ase_write(str(path), frames)


# -- refine -------------------------------------------------------------------


def _refine(args) -> int:
    selection_dir = Path(args.selection)
    selection = _load_json(selection_dir / "pars_selection.json")
    config = ParsConfig.from_dict(selection["config"])
    overrides = {}
    if args.images is not None:
        overrides["refinement_images"] = args.images
    if args.fmax_neb is not None:
        overrides["refinement_fmax_neb"] = args.fmax_neb
    if args.uncertainty_refusal is not None:
        overrides["refinement_uncertainty_refusal"] = args.uncertainty_refusal
    if overrides:
        config = dataclasses.replace(config, **overrides)

    representative_ids: list[str] = []
    for channel in selection["channels"].values():
        representative_ids.extend(channel.get("refinement_order", channel["representatives"]))
    if args.only:
        if args.only not in representative_ids:
            raise ValueError(f"{args.only} is not a representative in this selection")
        representative_ids = [args.only]

    refined_dir = selection_dir / "refined"
    replacement_names = (
        {f"refined-{cid}.json" for cid in representative_ids} if args.redo else set()
    )
    for path in sorted(refined_dir.glob("refined-*.json")):
        recorded = _load_json(path)
        replacing = path.name in replacement_names
        if not recorded.get("conditions") and not replacing:
            raise ValueError(
                f"{path.name}: no recorded conditions; redo all legacy records before aggregation"
            )
        if not replacing:
            _verify_refined_artifact(path, recorded)
    factory = _resolve_factory(args.calculator_factory)
    calculator = factory()
    metadata_path = getattr(args, "calculator_metadata", None)
    identity = calculator_identity(
        factory, calculator, _load_json(metadata_path) if metadata_path else None
    )
    seeds = {
        cid: selection_dir / "representatives" / f"rep-{cid}.extxyz"
        for cid in representative_ids
    }
    conditions = {
        cid: _refine_conditions(seeds[cid], args.calculator_factory, config, identity)
        for cid in representative_ids
    }
    # Selected records must match to resume, or refuse before force evaluations.
    # Different cohorts created by
    # --only/--redo are explicitly separated in refined_groups.json. The per-candidate
    # JSON is the commit marker and binds the geometry content: an extxyz
    # without a matching committed hash is an interrupted or damaged result.
    todo, skipped, conflicts = [], [], []
    for cid in representative_ids:
        record_path = refined_dir / f"refined-{cid}.json"
        if args.redo or not record_path.exists():
            todo.append(cid)
            continue
        recorded = _load_json(record_path).get("conditions")
        if recorded is None:
            conflicts.append(f"{cid}: no recorded conditions")
        elif recorded != conditions[cid]:
            differing = _differing_conditions(recorded, conditions[cid])
            conflicts.append(f"{cid}: differs in {differing}")
        else:
            skipped.append(cid)
    if conflicts:
        raise ValueError(
            "existing refined results were produced under other conditions; "
            "pass --redo to replace them:\n  " + "\n  ".join(conflicts)
        )

    if todo:
        backend = AseNebBackend(calculator)
        refined_dir.mkdir(parents=True, exist_ok=True)
    for candidate_id in todo:
        candidate = _candidate_from_seed(seeds[candidate_id])
        refined = backend.refine(candidate, config)
        artifact = _write_refined(refined_dir / f"refined-{candidate_id}.extxyz", refined)
        # one record per candidate as soon as it exists: a crash later in the
        # loop keeps every finished refinement for the next resume
        summary = refined.to_summary()
        summary["conditions"] = conditions[candidate_id]
        summary["geometry_artifact"] = artifact
        _write_json(refined_dir / f"refined-{candidate_id}.json", summary)
        print(f"  {candidate_id}: {refined.status.value} {refined.message}".rstrip())
    results = []
    for path in sorted(refined_dir.glob("refined-*.json")):
        recorded = _load_json(path)
        _verify_refined_artifact(path, recorded)
        results.append(recorded)
    if results:
        groups = condition_groups(results)
        _write_json(refined_dir / "refined_groups.json", groups)
        _write_json(refined_dir / "refined_paths.json", results)
        if not groups["homogeneous"]:
            print(
                f"  mixed refinement conditions: {len(groups['groups'])} cohorts; compare within refined_groups.json groups"
            )
    print(
        f"pars refine: {len(todo)} refined, {len(skipped)} already done under the "
        f"same conditions; {len(results)} recorded under {refined_dir}"
    )
    return 0


def _refine_conditions(
    seed: Path, calculator_factory: str, config: ParsConfig, identity=None
) -> dict:
    """Freeze implementation, numerical environment, model and candidate input."""
    return {
        "pars_version": __version__,
        "seed_sha256": hashlib.sha256(seed.read_bytes()).hexdigest(),
        "calculator_factory": calculator_factory,
        "calculator_identity": identity,
        "pars_implementation_sha256": implementation_sha256(),
        "numerical_versions": numerical_versions(),
        "refinement": {
            name: value
            for name, value in config.to_dict().items()
            if name.startswith("refinement_")
        },
    }


def _differing_conditions(recorded: dict, current: dict) -> list[str]:
    names = []
    for key in sorted(set(recorded) | set(current)):
        if key == "refinement":
            old, new = recorded.get(key) or {}, current.get(key) or {}
            names.extend(k for k in sorted(set(old) | set(new)) if old.get(k) != new.get(k))
        elif recorded.get(key) != current.get(key):
            names.append(key)
    return names


def _seed_state(provenance: Provenance, where: str) -> dict:
    if provenance.electronic_state is None:
        raise ValueError(
            f"{where}: no electronic state in provenance; a refinement seed must name "
            "the charge/spin of the PES it was found on (re-run pars select)"
        )
    return electronic.checked_state(provenance.electronic_state, where)


def _candidate_from_seed(path: Path) -> PathCandidate:
    frames = ase_read(str(path), index=":")
    info = frames[0].info
    provenance = Provenance.from_dict(json.loads(info["pars_provenance"]))
    state = _seed_state(provenance, path.name)
    for index, frame in enumerate(frames):
        if electronic.info_state(frame.info, f"{path.name} image {index}") != state:
            raise ValueError(
                f"{path.name}: image {index} electronic state "
                f"{electronic.info_state(frame.info, path.name)} differs from provenance {state}"
            )
    energies = []
    for frame in frames:
        energy = frame.get_potential_energy() if frame.calc is not None else 0.0
        energies.append(float(energy))
    return PathCandidate(
        candidate_id=info["candidate_id"],
        channel_id=info["channel_id"],
        trace_id=provenance.attempt_id,
        numbers=frames[0].get_atomic_numbers(),
        positions=np.asarray([f.positions for f in frames], dtype=float),
        energies=np.asarray(energies, dtype=float),
        uncertainty=None,
        signatures=None,
        provenance=dataclasses.replace(
            provenance, extra={**provenance.extra, "seed_path": str(path)}
        ),
    )


def _write_refined(path: Path, refined: RefinedPath) -> dict:
    frames = []
    if refined.positions is not None and refined.numbers is not None:
        for index, positions in enumerate(refined.positions):
            atoms = Atoms(numbers=refined.numbers, positions=positions)
            atoms.info = {
                "purpose": "pars_refined_path",
                "candidate_id": refined.candidate_id,
                "status": refined.status.value,
                "image": index,
            }
            if isinstance(refined.provenance.get("electronic_state"), dict):
                electronic.stamp(atoms.info, refined.provenance["electronic_state"])
            if refined.energies is not None:
                atoms.calc = SinglePointCalculator(atoms, energy=float(refined.energies[index]))
            frames.append(atoms)
    # Preserve the old artifact until the new one has been fully written. If
    # interrupted after replace but before the JSON commit, its old hash will
    # fail verification; a partial/new geometry can never be resumed silently.
    handle = tempfile.NamedTemporaryFile(dir=path.parent, suffix=".extxyz", delete=False)
    temporary = Path(handle.name)
    handle.close()
    try:
        if frames:
            ase_write(str(temporary), frames, format="extxyz")
        with temporary.open("rb") as stream:
            os.fsync(stream.fileno())
        digest = hashlib.sha256(temporary.read_bytes()).hexdigest()
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return {"filename": path.name, "sha256": digest, "frames": len(frames)}


def _verify_refined_artifact(record_path: Path, recorded: dict) -> None:
    """Verify the committed geometry before resuming or aggregating any cohort."""
    candidate_id = recorded.get("candidate_id")
    artifact = recorded.get("geometry_artifact")
    expected = f"refined-{candidate_id}.extxyz"
    if record_path.name != f"refined-{candidate_id}.json":
        raise ValueError(f"{record_path.name}: candidate identity mismatch; pass --redo")
    if not isinstance(artifact, dict) or artifact.get("filename") != expected:
        raise ValueError(f"{record_path.name}: no valid geometry artifact binding; pass --redo")
    path = record_path.parent / expected
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != artifact.get(
        "sha256"
    ):
        raise ValueError(
            f"{record_path.name}: missing or changed geometry artifact; pass --redo"
        )
    count = artifact.get("frames")
    if type(count) is not int or count < 0:
        raise ValueError(f"{record_path.name}: invalid geometry frame count; pass --redo")
    if count == 0:
        if path.stat().st_size or recorded.get("status") != "failed":
            raise ValueError(
                f"{record_path.name}: invalid empty geometry artifact; pass --redo"
            )
        return
    try:
        frames = ase_read(str(path), index=":")
    except Exception as error:
        raise ValueError(
            f"{record_path.name}: unreadable geometry artifact; pass --redo"
        ) from error
    if len(frames) != count or any(
        f.info.get("candidate_id") != candidate_id
        or f.info.get("status") != recorded.get("status")
        or f.info.get("image") != index
        for index, f in enumerate(frames)
    ):
        raise ValueError(
            f"{record_path.name}: geometry identity/frame count mismatch; pass --redo"
        )


def _resolve_factory(spec: str):
    module_name, _, attribute = spec.partition(":")
    if not module_name or not attribute:
        raise ValueError(f"--calculator-factory must be 'module:attribute', got {spec!r}")
    module = importlib.import_module(module_name)
    factory = getattr(module, attribute)
    if not callable(factory):
        raise ValueError(f"{spec} did not resolve to a callable")
    return factory


# -- small helpers ------------------------------------------------------------


def _load_json(path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def _write_json(path, payload) -> None:
    from .iox import write_json_atomic

    write_json_atomic(path, payload)


def _write_jsonl(path, rows) -> None:
    from .iox import write_jsonl_atomic

    write_jsonl_atomic(path, rows)
