"""Discipline mirrors pars/tests/conftest.py: tests never read the
repository's ``runs/`` and never fabricate run artifacts. The builders here
write tiny synthetic source trees and geometries into pytest's tmp_path;
nothing imports prrs at the schema/io layer, so these suites run on any
interpreter with numpy/ase-free deps. The shooting tests that land with M3
will import prrs (allowed one-way for this module, resolved from the parent
checkout via pyproject ``pythonpath``) and use small analytic potentials,
never a neural checkpoint.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from mature.schema import MaturationProtocol, MicrostateRef, PathProposal

HEX = "0123456789abcdef"


def fake_sha(seed: str) -> str:
    return hashlib.sha256(seed.encode()).hexdigest()


def make_source_run(root: Path) -> dict:
    """A minimal fake source run: two files whose hashes pin the protocol."""
    (root / "run").mkdir(parents=True, exist_ok=True)
    files = {}
    for name, content in (
        ("network.json", '{"chemical_nodes": []}\n'),
        ("config.json", "{}\n"),
    ):
        (root / "run" / name).write_text(content, encoding="utf-8")
        files[name] = hashlib.sha256(content.encode()).hexdigest()
    return {"root": root / "run", "files": files}


def make_sides() -> dict:
    """Two sides with deliberately different bond graphs (P1-shaped: one bond moved)."""
    side_a = MicrostateRef(
        node_id="c0000",
        microstate_id="c0000/m0000",
        chemical_key=fake_sha("shared-chemical-key"),
        structure_path="chemical/c0000_m0000.extxyz",
        structure_sha256=fake_sha("structure-A"),
        edges=[[0, 1], [1, 2]],
        mapped_broken=[[1, 2]],
        mapped_formed=[[0, 2]],
    )
    side_b = MicrostateRef(
        node_id="c0000",
        microstate_id="c0000/m0001",
        chemical_key=fake_sha("shared-chemical-key"),
        structure_path="chemical/c0000_m0001.extxyz",
        structure_sha256=fake_sha("structure-B"),
        edges=[[0, 1], [0, 2]],
        mapped_broken=[[1, 2]],
        mapped_formed=[[0, 2]],
    )
    return {"A": side_a, "B": side_b}


def make_protocol(
    source_root: Path, source_files: dict, batch_role: str = "explore", **overrides
):
    kwargs = dict(
        batch_role=batch_role,
        title="synthetic batch",
        source={
            "run_dir": str(source_root),
            "implementation_sha256": fake_sha("source-run-fingerprint"),
            "model_factory": "mature.tests:analytic_factory",
            "model_kwargs": {"device": "cpu"},
            "files": dict(source_files),
        },
        regions={"sides": {side: ref.to_dict() for side, ref in make_sides().items()}},
        dynamics={
            "temperature_K": 300.0,
            "timestep_fs": 0.25,
            "first_entry_k": 160,
            "free_step_cap": 1280,
        },
        candidates={},
        shooting={
            "batch_role": batch_role,
            "seed_base": 100000,
            **(
                {"shots_per_candidate": 16} if batch_role in ("zero_measure", "explore") else {}
            ),
            **(
                {"confirm_shots": 100, "selection_rule": "band-first, then distance to 0.5"}
                if batch_role == "confirm"
                else {}
            ),
            **(
                {"sidecheck_shots": 50, "selection_rule": "frozen seed draw"}
                if batch_role == "sidecheck"
                else {}
            ),
            **({"drift_tolerance_eV_atom": 0.01} if batch_role != "zero_measure" else {}),
        },
        budget={"max_shots": 160},
    )
    kwargs.update(overrides)
    return MaturationProtocol.from_dict(kwargs)


def write_geometry(store_dir: Path, name: str, n_atoms: int = 4) -> tuple:
    """A tiny extxyz-shaped geometry file; content is arbitrary for io tests."""
    path = store_dir / "geometries" / f"{name}.extxyz"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps({"n_atoms": n_atoms, "name": name}) + "\n"
    path.write_text(payload, encoding="utf-8")
    return path, hashlib.sha256(payload.encode()).hexdigest()


def make_proposal(
    store_dir: Path, proposal_id: str, batch_role: str = "explore"
) -> PathProposal:
    from mature.schema import sha256_object

    _, geometry_sha = write_geometry(store_dir, proposal_id)
    sides = make_sides()
    return PathProposal(
        proposal_id=proposal_id,
        batch_role=batch_role,
        channel_id="rc0000",
        event_key=fake_sha("event"),
        tube_id="tube0000",
        source_candidate_id="pc-rc0000-t000038",
        source_run_dir="runs/synthetic",
        source_attempt_id="t000038",
        parent_frame_start=0,
        parent_frame_stop=16,
        proposal_kind="observed_frame",
        arc_fraction=2 / 6,
        frame_index=5,
        geometry_file=f"geometries/{proposal_id}.extxyz",
        geometry_sha256=geometry_sha,
        n_atoms=4,
        region_graph_sha256s={
            side: sha256_object({"edges": [list(pair) for pair in ref.edges]})
            for side, ref in sides.items()
        },
    )
