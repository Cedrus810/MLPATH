"""Durable output for one maturation batch (plan sections 3.2 and M1/M4).

The layout under the output directory is the plan's:

    manifest.json             schema, sources, protocol/model hashes
    protocol.json             the frozen MaturationProtocol, canonical JSON
    proposals.jsonl           one record per candidate
    geometries/<id>.extxyz    candidate geometries (all-atom)
    tasks.jsonl               the deterministic proposal/shot/seed allocation
    shots/<pid>/<sid>/        per-shot trial output
    commits/<pid>/<sid>.json  the atomic completion marker (plan 3.2)
    diagnostics.jsonl         per-proposal diagnosis
    interrupted/              half-written shot dirs, quarantined, never deleted

Completion is defined by the commit marker, nothing else: a shot directory
without a marker is unfinished work, gets quarantined on the next execution,
and the task is redone. Aggregates are rebuilt from markers only, so "done"
is never inferred from files that merely look finished. Resume verifies the
protocol hash, the source file hashes and the running implementation -- a
changed model, region definition or checkpoint cannot continue in this
directory (M1 acceptance).
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path

from . import __version__
from .schema import (
    SCHEMA_VERSION,
    SchemaError,
    ShotRecord,
    ShootingEstimate,
    canonical_json,
    sha256_object,
    utc_now,
)

MANIFEST_NAME = "manifest.json"
PROTOCOL_NAME = "protocol.json"
PROPOSALS_NAME = "proposals.jsonl"
TASKS_NAME = "tasks.jsonl"
DIAGNOSTICS_NAME = "diagnostics.jsonl"


class OutputConflict(ValueError):
    """A record or environment conflicts with what this directory already holds."""


def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _dumps_json(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)


def read_json_strict(path):
    """Parse JSON refusing NaN/Infinity literals anywhere (they are not values)."""
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig")
    try:
        return json.loads(text, parse_constant=_reject_constant)
    except ValueError as exc:
        raise ValueError(f"{path}: {exc}") from None


def _reject_constant(name):
    raise ValueError(f"JSON constant {name!r} is not a value this layer accepts")


def write_json_atomic(path, obj) -> None:
    """Write canonical JSON via tmp + fsync + os.replace: readers never see half a file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        handle.write(_dumps_json(obj))
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def append_jsonl(path, obj) -> None:
    """Append one line and fsync: a killed process leaves a whole line or none."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    line = canonical_json(obj)
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def read_jsonl(path):
    path = Path(path)
    if not path.exists():
        return []
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            records.append(json.loads(line, parse_constant=_reject_constant))
    return records


def verify_file_hashes(pairs) -> None:
    """Every (path -> sha256) pair must hold; missing or differing is a refusal.

    This is the E4 gate's substance: a checkpoint (or any input) whose content
    changed under the same name cannot silently continue a batch.
    """
    for path, expected in sorted(pairs.items()):
        path = Path(path)
        if not path.exists():
            raise OutputConflict(f"source file missing: {path}")
        actual = sha256_file(path)
        if actual != expected:
            raise OutputConflict(
                f"source file content changed under the same name: {path} "
                f"(recorded {expected[:12]}..., now {actual[:12]}...)"
            )


class OutputStore:
    """Owns one output directory; every write goes through it."""

    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)

    # ---- lifecycle ----------------------------------------------------------

    @classmethod
    def initialize(cls, output_dir, protocol, prrs_implementation_sha256=None) -> "OutputStore":
        """Create a fresh batch directory; refuse anything that already holds data."""
        output_dir = Path(output_dir)
        if output_dir.exists() and any(output_dir.iterdir()):
            raise OutputConflict(
                f"refusing to initialise into non-empty {output_dir}: a maturation "
                "batch owns its directory from birth (resume instead: OutputStore.resume)"
            )
        _verify_protocol_sources(protocol)
        for sub in ("geometries", "shots", "commits", "interrupted"):
            (output_dir / sub).mkdir(parents=True, exist_ok=True)
        write_json_atomic(output_dir / PROTOCOL_NAME, protocol.to_dict())
        write_json_atomic(
            output_dir / MANIFEST_NAME,
            {
                "kind": "mlpath-maturation",
                "schema_version": SCHEMA_VERSION,
                "mature_version": __version__,
                "batch_role": protocol.batch_role,
                "protocol_sha256": protocol.sha256,
                "prrs_implementation_sha256": prrs_implementation_sha256,
                "created_utc": utc_now(),
                "source": protocol.source.to_dict(),
            },
        )
        return cls(output_dir)

    @classmethod
    def resume(cls, output_dir, protocol, prrs_implementation_sha256=None) -> "OutputStore":
        """Reopen an existing batch; refuse on any identity mismatch (plan M4).

        Refusals, each with the reason: no manifest (not a maturation
        directory), a different protocol hash (a changed model, region
        definition, budget or seed table -- plan M1: "换模型/区域定义不能续在
        原目录"), a stored protocol that differs from the supplied object, any
        source or geometry file whose content changed, or a different prrs
        implementation than the batch started under (a mixed batch -- this
        project already retired one of those).
        """
        output_dir = Path(output_dir)
        manifest_path = output_dir / MANIFEST_NAME
        if not manifest_path.exists():
            raise OutputConflict(f"resume refused: no {MANIFEST_NAME} in {output_dir}")
        manifest = read_json_strict(manifest_path)
        if manifest.get("kind") != "mlpath-maturation":
            raise OutputConflict(f"resume refused: {output_dir} is not a maturation directory")
        if manifest.get("schema_version") != SCHEMA_VERSION:
            raise OutputConflict(
                f"resume refused: manifest schema {manifest.get('schema_version')!r} "
                f"!= {SCHEMA_VERSION}"
            )
        if manifest.get("protocol_sha256") != protocol.sha256:
            raise OutputConflict(
                "resume refused: the protocol differs from the one this directory was "
                "started with (model, region definition, dynamics, rules or budget); "
                "start a new directory instead"
            )
        stored = read_json_strict(output_dir / PROTOCOL_NAME)
        if stored != protocol.to_dict():
            raise OutputConflict(
                "resume refused: stored protocol.json differs from the supplied "
                "protocol despite a matching hash -- refusing on principle"
            )
        _verify_protocol_sources(protocol)
        _verify_geometry_hashes(output_dir, protocol)
        recorded = manifest.get("prrs_implementation_sha256")
        if (
            prrs_implementation_sha256 is not None
            and recorded is not None
            and recorded != prrs_implementation_sha256
        ):
            raise OutputConflict(
                "resume refused: the prrs implementation fingerprint changed since "
                "this batch started "
                f"({str(recorded)[:12]}... -> {str(prrs_implementation_sha256)[:12]}...); "
                "continuing would mix implementations"
            )
        return cls(output_dir)

    # ---- proposals ----------------------------------------------------------

    def add_proposal(self, proposal) -> None:
        """Record a candidate; an identical re-record is a no-op, a clash refuses."""
        record = proposal.to_dict()
        existing = {
            data["proposal_id"]: data for data in read_jsonl(self.output_dir / PROPOSALS_NAME)
        }
        prior = existing.get(proposal.proposal_id)
        if prior is not None:
            if canonical_json(prior) == canonical_json(record):
                return
            raise OutputConflict(
                f"proposal id clash: {proposal.proposal_id} is already recorded with "
                "different content; stable ids must be content-derived, not reused"
            )
        geometry = self.output_dir / proposal.geometry_file
        if not geometry.exists():
            raise OutputConflict(f"proposal geometry missing: {geometry}")
        if sha256_file(geometry) != proposal.geometry_sha256:
            raise OutputConflict(
                f"proposal geometry hash mismatch: {geometry} does not hash to "
                f"{proposal.geometry_sha256[:12]}..."
            )
        append_jsonl(self.output_dir / PROPOSALS_NAME, record)

    def proposals(self) -> list:
        from .schema import PathProposal

        return [PathProposal.from_dict(d) for d in read_jsonl(self.output_dir / PROPOSALS_NAME)]

    # ---- tasks and completion ----------------------------------------------

    def build_task_table(self, protocol, proposals) -> list:
        """The deterministic proposal/shot/seed allocation; written once.

        Seeds are ``seed_base + index`` over the tasks sorted by
        (proposal_id, shot index) -- independent of execution order and of
        Python's randomised hash (plan M4: "不依赖执行顺序或 Python 的随机
        hash"). An existing identical table is a no-op; a differing one
        refuses.
        """
        shots = protocol.shooting.shots_for_role()
        if shots is None:
            raise SchemaError(
                f"batch role {protocol.shooting.batch_role!r} does not define a shot "
                "count for the task table"
            )
        tasks = []
        for proposal in sorted(proposals, key=lambda p: p.proposal_id):
            for shot_index in range(shots):
                shot_id = f"s{shot_index:04d}"
                tasks.append(
                    {
                        "task_id": f"{proposal.proposal_id}#{shot_id}",
                        "proposal_id": proposal.proposal_id,
                        "shot_id": shot_id,
                        "seed": protocol.shooting.seed_base + len(tasks),
                        "batch_role": protocol.shooting.batch_role,
                    }
                )
        table_path = self.output_dir / TASKS_NAME
        if table_path.exists():
            stored = read_jsonl(table_path)
            if [canonical_json(t) for t in stored] != [canonical_json(t) for t in tasks]:
                raise OutputConflict(
                    "task table conflict: the stored allocation differs from a "
                    "rebuild; proposals or protocol changed under the same directory"
                )
            return tasks
        for task in tasks:
            append_jsonl(table_path, task)
        return tasks

    def tasks(self) -> list:
        return read_jsonl(self.output_dir / TASKS_NAME)

    def pending_tasks(self) -> list:
        done = self.done_task_ids()
        return [t for t in self.tasks() if t["task_id"] not in done]

    def done_task_ids(self) -> set:
        """Completion comes from commit markers only (plan M4)."""
        done = set()
        commits = self.output_dir / "commits"
        if not commits.exists():
            return done
        for marker_path in sorted(commits.rglob("*.json")):
            marker = read_json_strict(marker_path)
            done.add(marker["task_id"])
        return done

    # ---- committing results -------------------------------------------------

    def isolate_interrupted(self, proposal_id, shot_id) -> str | None:
        """Quarantine a shot directory that has no commit marker; never delete."""
        shot_dir = self.output_dir / "shots" / proposal_id / shot_id
        marker = self.output_dir / "commits" / proposal_id / f"{shot_id}.json"
        if not shot_dir.exists() or marker.exists():
            return None
        destination = (
            self.output_dir
            / "interrupted"
            / proposal_id
            / (f"{shot_id}.{utc_now().replace(':', '').replace('+', 'p')}")
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(shot_dir), str(destination))
        return str(destination)

    def commit(self, record: ShotRecord) -> None:
        """Write the shot's result, then its marker -- in that order (plan M4).

        A crash between the two leaves an unmarked shot directory: unfinished
        work that the next execution quarantines and redoes. A marker with the
        same content hash is an idempotent re-commit; different content for
        the same task refuses.
        """
        self._check_batch_conditions(record)
        record_dict = record.to_dict()
        record_sha = sha256_object(record_dict)
        result_path = (
            self.output_dir / "shots" / record.proposal_id / record.shot_id / ("result.json")
        )
        marker_path = (
            self.output_dir / "commits" / record.proposal_id / (f"{record.shot_id}.json")
        )
        if marker_path.exists():
            marker = read_json_strict(marker_path)
            if marker.get("record_sha256") == record_sha:
                return
            raise OutputConflict(
                f"commit conflict: {marker_path} already records task "
                f"{record.proposal_id}/{record.shot_id} with different content"
            )
        write_json_atomic(result_path, record_dict)
        write_json_atomic(
            marker_path,
            {
                "task_id": f"{record.proposal_id}#{record.shot_id}",
                "proposal_id": record.proposal_id,
                "shot_id": record.shot_id,
                "record_sha256": record_sha,
                "record": record_dict,
                "committed_utc": utc_now(),
            },
        )

    def _check_batch_conditions(self, record: ShotRecord) -> None:
        """A shot counts only for the task it was allocated under this batch.

        Estimates group by proposal id alone, so a record from another
        protocol, seed, role or implementation would otherwise occupy a task
        and mix into one estimate with this batch's shots.
        """
        task_id = f"{record.proposal_id}#{record.shot_id}"
        task = next((t for t in self.tasks() if t["task_id"] == task_id), None)
        if task is None:
            raise OutputConflict(f"commit refused: {task_id} is not in this batch's task table")
        manifest = read_json_strict(self.output_dir / MANIFEST_NAME)
        expected = {
            "seed": task["seed"],
            "batch_role": task["batch_role"],
            "protocol_sha256": manifest["protocol_sha256"],
        }
        if manifest.get("prrs_implementation_sha256") is not None:
            expected["implementation_sha256"] = manifest["prrs_implementation_sha256"]
        differing = sorted(k for k, v in expected.items() if getattr(record, k) != v)
        if differing:
            raise OutputConflict(
                f"commit refused: {task_id} ran under other conditions than this batch "
                f"allocated ({', '.join(differing)})"
            )

    def committed_records(self) -> list:
        records = []
        for marker_path in sorted((self.output_dir / "commits").rglob("*.json")):
            marker = read_json_strict(marker_path)
            records.append(ShotRecord.from_dict(marker["record"]))
        return records

    # ---- aggregation --------------------------------------------------------

    def estimates(self, band=(0.4, 0.6)) -> dict:
        """One ShootingEstimate per proposal with commits, rebuilt from markers."""
        by_proposal = {}
        for record in self.committed_records():
            by_proposal.setdefault(record.proposal_id, []).append(record)
        return {
            proposal_id: ShootingEstimate.from_records(records, band=band)
            for proposal_id, records in sorted(by_proposal.items())
        }

    def add_diagnosis(self, diagnosis) -> None:
        append_jsonl(self.output_dir / DIAGNOSTICS_NAME, diagnosis.to_dict())

    def diagnoses(self) -> list:
        return read_jsonl(self.output_dir / DIAGNOSTICS_NAME)


# ---- verification helpers ---------------------------------------------------


def _verify_protocol_sources(protocol) -> None:
    source = protocol.source
    root = Path(source.run_dir)
    if not root.is_dir():
        raise OutputConflict(f"source run directory missing: {root}")
    verify_file_hashes({root / rel: sha for rel, sha in sorted(source.files.items())})
    if source.model_checkpoint_path and source.model_checkpoint_sha256:
        verify_file_hashes({Path(source.model_checkpoint_path): source.model_checkpoint_sha256})


def _verify_geometry_hashes(output_dir, protocol) -> None:
    """Recheck every recorded proposal geometry against its hash (E4)."""
    proposals_path = output_dir / PROPOSALS_NAME
    if not proposals_path.exists():
        return
    verify_file_hashes(
        {
            output_dir / data["geometry_file"]: data["geometry_sha256"]
            for data in read_jsonl(proposals_path)
        }
    )
