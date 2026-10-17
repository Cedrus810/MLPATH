"""Geometry/metadata-only seed probes; no PES, MD, NEB or historical run reads."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "pars/src")]

import numpy as np  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.io import read, write  # noqa: E402
from pars.cli import _write_representative  # noqa: E402
from pars.config import ParsConfig  # noqa: E402
from pars.schema import PathCandidate, Provenance  # noqa: E402
from prrs.search import _dividing_frame  # noqa: E402


def write_rows(directory, phases, events):
    directory.mkdir()
    frames = [
        Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [0, 1 + i * 0.1, 0]])
        for i in range(len(phases))
    ]
    write(str(directory / "response.traj"), frames)
    rows = [
        {"frame": i, "phase": phase, "events": event}
        for i, (phase, event) in enumerate(zip(phases, events))
    ]
    (directory / "observations.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in rows)
    )
    write(str(directory / "raw_endpoint.extxyz"), frames[-1])
    return frames


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    output = parser.parse_args().out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = [
        "src/prrs/search.py",
        "src/prrs/runner.py",
        "pars/src/pars/cli.py",
        "pars/src/pars/paths/resample.py",
        "benchmarks/project_audit/seed_stage_probes.py",
    ]
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    empty = {"formed": [], "broken": []}
    forced_dir = output / "first-forced-event"
    phases = ("source", "perturbed", "forced", "free_start", "free")
    events = (
        empty,
        empty,
        {"formed": [[0, 1]], "broken": []},
        empty,
        {"formed": [[1, 2]], "broken": [[0, 1]]},
    )
    write_rows(forced_dir, phases, events)
    _, description, _, tangent_source = _dividing_frame(forced_dir)
    assert description == "topology_switch:forced:2"

    no_event_dir = output / "forced-no-event"
    phases_without_event = ("source", "perturbed", "forced", "forced", "free_start", "free")
    frames = write_rows(
        no_event_dir, phases_without_event, (empty,) * len(phases_without_event)
    )
    _, fallback_description, tangent, fallback_source = _dividing_frame(no_event_dir)
    last_forced_step_norm = float(np.linalg.norm(frames[3].positions - frames[2].positions))
    assert tangent is None and last_forced_step_norm > 0

    # An artificial delivery gap of 1 A followed by a 0.05 A free segment.
    positions = np.repeat(np.array([[[0.0, 0, 0], [0, 2.0, 0], [0, 0, 1.0]]]), 3, axis=0)
    positions[:, 0, 0] = [-1.0, 0.0, 0.05]
    provenance = Provenance(
        run_dir="synthetic-not-a-run",
        attempt_id="gap",
        channel_id="synthetic",
        event_key=None,
        reactant_state_id=None,
        product_state_id=None,
        frame_start=0,
        frame_stop=4,
        electronic_state={"charge": 0, "spin": 1},
        extra={"kept_frames": [0, 2, 3]},
    )
    path = PathCandidate(
        candidate_id="delivery-gap",
        channel_id="synthetic",
        trace_id="gap",
        numbers=np.ones(3, dtype=int),
        positions=positions,
        energies=np.zeros(3),
        uncertainty=None,
        signatures=(frozenset(),) * 3,
        provenance=provenance,
    )
    seed_path = output / "representative.extxyz"
    _write_representative(seed_path, path, ParsConfig(refinement_remove_rigid_motion=False))
    seed = read(str(seed_path), index=":")
    x = [float(atoms.positions[0, 0]) for atoms in seed]
    interior_in_free_segment = sum(value >= 0.0 for value in x[1:-1])
    assert interior_in_free_segment == 0
    report = {
        "kind": "geometry-and-metadata-seed-audit",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": hashes,
        "limits": "Synthetic helper inputs only; event labels are supplied, not physically measured. No PES/MD/NEB or real saddle outcomes.",
        "first_switch": {
            "description": description,
            "tangent_source": tangent_source,
            "later_free_switch_frame": 4,
        },
        "forced_phase_fallback": {
            "description": fallback_description,
            "tangent_missing_without_proposal_fallback": tangent is None,
            "tangent_source": fallback_source,
            "actual_last_forced_step_norm_A": last_forced_step_norm,
        },
        "delivery_gap_resampling": {
            "source_to_free_start_A": 1.0,
            "free_segment_A": 0.05,
            "images": len(seed),
            "interior_images_in_free_segment": interior_in_free_segment,
            "moving_atom_x_A": x,
            "rigid_removal": False,
        },
    }
    assert hashes == {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
