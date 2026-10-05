"""Regressions from the 2026-10-05 re-review (docs/REREVIEW_2026-10-05.md).

Electronic-state handoff, uncertainty of force-difference Hessian samples and
the order of authoritative-record checks versus classification filtering.
"""

from __future__ import annotations

import dataclasses
import json
import textwrap

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read as ase_read
from ase.io import write as ase_write

from conftest import build_run
from pars.adapters.prrs_run import PrrsRun
from pars.cli import main
from pars.config import ParsConfig
from pars.paths.ingest import channel_members_not_ingested, ingest_run, skipped_attempts
from pars.refinement.ase_neb import AseNebBackend
from pars.schema import RefinementStatus
from test_review_regressions import CONFIG, QuarticWell, quartic_candidate

# -- electronic state ---------------------------------------------------------

CHARGED_FACTORY = textwrap.dedent(
    """
    import numpy as np
    from ase.calculators.calculator import Calculator, all_changes

    SEEN = set()


    class ChargedWell(Calculator):
        \"\"\"Quartic double well whose depth scales with the structure's charge.\"\"\"

        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            SEEN.add((self.atoms.info.get("charge"), self.atoms.info.get("spin")))
            scale = 1 + abs(self.atoms.info.get("charge", 0))
            x, y, z = self.atoms.positions[0]
            self.results = {
                "energy": scale * (0.25 * x * x * (x - 2) ** 2 + 0.5 * (y * y + z * z)),
                "forces": scale * np.array([[-x * (x - 2) * (x - 1), -y, -z]]),
            }


    def make():
        return ChargedWell()
    """
)


def _charged_run(root, charge=-1, spin=3):
    """One reactive single-atom trial along the double-well axis."""
    run = build_run(root)
    network = json.loads((run / "network.json").read_text())
    network["reaction_channels"] = network["reaction_channels"][:1]
    network["reactions"] = network["reactions"][:1]
    (run / "network.json").write_text(json.dumps(network))
    ledger = [json.loads(line) for line in (run / "attempts.jsonl").read_text().splitlines()]
    (run / "attempts.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in ledger if row["attempt_id"] == "t000001")
    )
    manifest = json.loads((run / "config.json").read_text())
    manifest["config"] = {
        "charge_sensitive": True,
        "total_charge": charge,
        "multiplicity": spin,
    }
    (run / "config.json").write_text(json.dumps(manifest))
    trial = run / "trials/t000001"
    phases = ["source", "perturbed", "free_start", "free", "free", "free", "quench", "quenched"]
    xs = [0, 0.1, 0.2, 0.7, 1.2, 1.8, 2, 2]
    frames, rows = [], []
    for index, (phase, x) in enumerate(zip(phases, xs)):
        atoms = Atoms(
            "H", positions=[[x, 0, 0]], info={"phase": phase, "charge": charge, "spin": spin}
        )
        energy = 0.5 * x * x * (x - 2) ** 2
        atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=np.zeros((1, 3)))
        frames.append(atoms)
        rows.append({"frame": index, "phase": phase, "potential_eV": energy, "edges": []})
    ase_write(str(trial / "response.traj"), frames)
    ase_write(str(trial / "endpoint.extxyz"), frames[-1])
    (trial / "observations.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    return run


def _rewrite_frames(path, edit):
    frames = ase_read(str(path), index=":")
    for index, atoms in enumerate(frames):
        edit(index, atoms.info)
    ase_write(str(path), frames)


def _set_manifest_state(run, charge, spin):
    manifest = json.loads((run / "config.json").read_text())
    manifest["config"] = {
        "charge_sensitive": True,
        "total_charge": charge,
        "multiplicity": spin,
    }
    (run / "config.json").write_text(json.dumps(manifest))


@pytest.fixture
def charged_factory(tmp_path, monkeypatch):
    (tmp_path / "pars_charged_factory.py").write_text(CHARGED_FACTORY)
    monkeypatch.syspath_prepend(str(tmp_path))
    import importlib
    import sys

    sys.modules.pop("pars_charged_factory", None)
    return importlib.import_module("pars_charged_factory")


def _select_charged(run, out):
    config = out.parent / f"{out.name}-config.json"
    config.write_text(
        json.dumps({"refinement_remove_rigid_motion": False, "refinement_fmax_neb": 0.001})
    )
    assert (
        main(["select", "--run-dir", str(run), "--out", str(out), "--config", str(config)]) == 0
    )


def test_charged_state_reaches_every_backend_evaluation_and_output(tmp_path, charged_factory):
    run = _charged_run(tmp_path / "src")
    selection = tmp_path / "sel"
    _select_charged(run, selection)
    seed = selection / "representatives/rep-pc-rc0000-t000001.extxyz"
    for frame in ase_read(str(seed), index=":"):
        assert (frame.info["charge"], frame.info["spin"]) == (-1, 3)
    candidates = [json.loads(line) for line in (selection / "candidates.jsonl").open()]
    assert candidates[0]["provenance"]["electronic_state"] == {"charge": -1, "spin": 3}

    argv = ["refine", "--selection", str(selection)]
    assert main([*argv, "--calculator-factory", "pars_charged_factory:make"]) == 0
    record = json.loads((selection / "refined/refined-pc-rc0000-t000001.json").read_text())
    assert charged_factory.SEEN == {(-1, 3)}
    assert record["status"] == "ci_neb_converged"
    assert record["barrier_eV"] == pytest.approx(0.5, abs=1e-5)
    assert record["provenance"]["electronic_state"] == {"charge": -1, "spin": 3}
    refined = ase_read(str(selection / "refined/refined-pc-rc0000-t000001.extxyz"), index=":")
    assert {(f.info["charge"], f.info["spin"]) for f in refined} == {(-1, 3)}


def test_frames_contradicting_the_manifest_state_are_rejected(tmp_path):
    run = _charged_run(tmp_path / "src")
    _rewrite_frames(
        run / "trials/t000001/response.traj",
        lambda index, info: info.update(charge=0) if index == 4 else None,
    )
    with pytest.raises(ValueError, match="electronic state"):
        PrrsRun(run).load_trace("t000001")


def test_endpoint_contradicting_the_manifest_state_is_rejected(tmp_path):
    run = _charged_run(tmp_path / "src")
    _rewrite_frames(run / "trials/t000001/endpoint.extxyz", lambda _, info: info.update(spin=1))
    with pytest.raises(ValueError, match="electronic state"):
        PrrsRun(run).load_trace("t000001")


def test_unstamped_frames_cannot_certify_a_nondefault_declared_state(tmp_path):
    run = _charged_run(tmp_path / "src")
    for name in ("response.traj", "endpoint.extxyz"):
        _rewrite_frames(
            run / "trials/t000001" / name,
            lambda _, info: [info.pop(k) for k in ("charge", "spin")],
        )
    with pytest.raises(ValueError, match="electronic state"):
        PrrsRun(run).load_trace("t000001")


def test_partially_stamped_trajectory_is_rejected(tmp_path):
    run = build_run(tmp_path)
    _rewrite_frames(
        run / "trials/t000001/response.traj",
        lambda index, info: info.pop("spin") if index == 2 else None,
    )
    with pytest.raises(ValueError, match="electronic state"):
        PrrsRun(run).load_trace("t000001")


def test_unstamped_frames_take_the_default_declared_state(tmp_path):
    run = build_run(tmp_path)
    for trial in ("t000001", "t000002"):
        for name in ("response.traj", "endpoint.extxyz"):
            _rewrite_frames(
                run / "trials" / trial / name,
                lambda _, info: [info.pop(k) for k in ("charge", "spin")],
            )
    traces = ingest_run(run, ParsConfig())
    assert {t.provenance.electronic_state["charge"] for t in traces} == {0}
    assert {t.provenance.electronic_state["spin"] for t in traces} == {1}


@pytest.mark.parametrize(
    "config", [None, {"total_charge": -1.0, "multiplicity": 1}, {"multiplicity": 0}]
)
def test_missing_or_invalid_manifest_state_is_rejected(tmp_path, config):
    run = build_run(tmp_path)
    manifest = json.loads((run / "config.json").read_text())
    if config is None:
        manifest.pop("config")
    else:
        manifest["config"] = {"charge_sensitive": True, **config}
    (run / "config.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="electronic state|total_charge|multiplicity"):
        PrrsRun(run)


@pytest.mark.parametrize(
    "edit",
    [
        lambda index, info: info.update(charge=0) if index == 3 else None,
        lambda index, info: info.pop("spin"),
        lambda index, info: info.update(
            pars_provenance=_without_state(info["pars_provenance"])
        ),
    ],
)
def test_seed_state_is_verified_before_refinement(tmp_path, charged_factory, edit):
    run = _charged_run(tmp_path / "src")
    selection = tmp_path / "sel"
    _select_charged(run, selection)
    _rewrite_frames(selection / "representatives/rep-pc-rc0000-t000001.extxyz", edit)
    argv = ["refine", "--selection", str(selection)]
    with pytest.raises(ValueError, match="electronic state"):
        main([*argv, "--calculator-factory", "pars_charged_factory:make"])
    assert charged_factory.SEEN == set()


def _without_state(payload):
    provenance = json.loads(payload)
    provenance.pop("electronic_state")
    return json.dumps(provenance)


def test_backend_without_a_declared_state_says_so():
    result = AseNebBackend(QuarticWell()).refine(quartic_candidate(), CONFIG)
    assert result.provenance["electronic_state"] == "unspecified"


# -- Hessian sampling uncertainty -----------------------------------------------


class TransverseUncertainty(QuarticWell):
    """The band and descents stay on the x axis; only transverse Hessian
    displacements (y or z off zero) see ``transverse``, optionally only near x0."""

    def __init__(self, transverse, near=None, axis=0.005):
        super().__init__(axis)
        self.transverse = transverse
        self.near = near
        self.transverse_samples = 0

    def calculate(self, *args, **kwargs):
        super().calculate(*args, **kwargs)
        x, y, z = self.atoms.positions[0]
        if (abs(y) > 1e-8 or abs(z) > 1e-8) and (
            self.near is None or abs(x - self.near) < 0.05
        ):
            self.transverse_samples += 1
            if self.transverse is None:
                self.results.pop("force_uncertainty_eV_A")
            else:
                self.results["force_uncertainty_eV_A"] = self.transverse


GATED = dataclasses.replace(
    CONFIG, refinement_fmax_neb=0.001, refinement_uncertainty_refusal=0.01
)


@pytest.mark.parametrize(
    ("transverse", "status"),
    [
        (float("nan"), RefinementStatus.FAILED),
        (-1.0, RefinementStatus.FAILED),
        (1000.0, RefinementStatus.UNRESOLVED_MLP),
        (None, RefinementStatus.UNRESOLVED_MLP),
    ],
)
def test_endpoint_hessian_samples_pass_the_uncertainty_gate(transverse, status):
    calculator = TransverseUncertainty(transverse)
    candidate = quartic_candidate()
    result = AseNebBackend(calculator).refine(candidate, GATED)
    assert calculator.transverse_samples > 0
    assert result.status is status
    assert result.barrier_eV is None
    assert "Hessian" in result.message
    check = result.provenance["endpoint_checks"]["reactant"]
    assert check["status"] in {"failed", "unresolved_mlp"}
    assert check["uncertainty_observations"] >= 1
    json.dumps(result.to_summary(), allow_nan=False)


def test_saddle_hessian_samples_pass_the_uncertainty_gate():
    calculator = TransverseUncertainty(1000.0, near=1.0)
    result = AseNebBackend(calculator).refine(quartic_candidate(), GATED)
    assert calculator.transverse_samples > 0
    assert result.status is RefinementStatus.UNRESOLVED_MLP
    assert result.barrier_eV is None
    saddle = result.provenance["saddle_check"]
    assert saddle["status"] == "unresolved_mlp"
    assert saddle["uncertainty_max_eV_A"] == pytest.approx(1000.0)
    assert "negative_modes" not in saddle
    assert result.uncertainty_max_eV_A == pytest.approx(1000.0)
    top = saddle["image_index"]
    assert result.positions[top, 0, 1] == 0.0  # displacements restored


def test_hessian_sample_uncertainty_enters_the_reported_maximum():
    for config in (GATED, dataclasses.replace(GATED, refinement_uncertainty_refusal=None)):
        result = AseNebBackend(TransverseUncertainty(0.007)).refine(quartic_candidate(), config)
        assert result.status is RefinementStatus.CI_NEB_CONVERGED
        assert result.barrier_eV == pytest.approx(0.25, abs=1e-5)
        assert result.uncertainty_max_eV_A == pytest.approx(0.007)
        branch = result.provenance["connectivity_check"]["branches"][0]
        assert branch["minimum_check"]["uncertainty_max_eV_A"] == pytest.approx(0.007)


# -- authoritative records before filtering --------------------------------------


def _edit_ledger(run, attempt_id, **fields):
    path = run / "attempts.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    for row in rows:
        if row["attempt_id"] == attempt_id:
            row.update(fields)
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))


@pytest.mark.parametrize(
    ("attempt_id", "fields"),
    [
        ("t000001", {"classification": "elastic"}),  # ledger hides a reactive result
        ("t000000", {"classification": "reactive"}),  # ledger promotes an elastic result
        ("t000001", {"status": "unconfirmed"}),
    ],
)
def test_ledger_result_conflicts_are_caught_before_filtering(tmp_path, attempt_id, fields):
    run = build_run(tmp_path)
    _edit_ledger(run, attempt_id, **fields)
    prrs = PrrsRun(run)
    for view in (ingest_run, skipped_attempts, channel_members_not_ingested):
        with pytest.raises(ValueError, match="ledger mismatch"):
            view(prrs, ParsConfig())


def test_consistent_elastic_trials_are_still_skipped(run_dir):
    run = PrrsRun(run_dir)
    assert [t.trace_id for t in ingest_run(run, ParsConfig())] == ["t000001", "t000002"]
    assert skipped_attempts(run, ParsConfig()) == {"completed/elastic": 1}
