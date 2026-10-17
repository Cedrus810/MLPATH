"""End-to-end CLI glue on the synthetic run: select -> files -> refine."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest
from ase.io import read as ase_read
from ase.io import write as ase_write

from pars.cli import main


FACTORY_MODULE = textwrap.dedent(
    """
    import numpy as np
    from ase.calculators.calculator import Calculator, all_changes


    class PerAtomDoubleWell(Calculator):
        \"\"\"Each atom carries independent x/y wells at 0 and 1.5, with
        x relative to its synthetic source position and harmonic z confinement.\"\"\"

        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            r = self.atoms.positions
            xy = r[:, :2].copy()
            xy[:, 0] -= 1.1 * np.arange(len(r))
            energy = float(np.sum(0.25 * xy**2 * (xy - 1.5) ** 2))
            energy += float(0.5 * np.sum(r[:, 2] ** 2))
            forces = np.zeros_like(r)
            forces[:, :2] = -xy * (xy - 1.5) * (xy - 0.75)
            forces[:, 2] = -r[:, 2]
            self.results["energy"] = energy
            self.results["forces"] = forces


    def make_calculator():
        return PerAtomDoubleWell()
    """
)


def _select(run_dir, out, *extra):
    # the synthetic run and PerAtomDoubleWell are external fields (frames move by
    # pure translation), so rigid-motion removal must be off for them
    config = Path(out).parent / f"{Path(out).name}-config.json"
    config.write_text(json.dumps({"refinement_remove_rigid_motion": False}))
    argv = ["select", "--run-dir", str(run_dir), "--out", str(out), "--config", str(config)]
    return main([*argv, *extra])


@pytest.fixture
def selection_dir(run_dir, tmp_path) -> Path:
    out = tmp_path / "selection"
    exit_code = _select(run_dir, out)
    assert exit_code == 0
    return out


def test_select_writes_the_full_record(selection_dir):
    selection = json.loads((selection_dir / "pars_selection.json").read_text(encoding="utf-8"))
    assert selection["counts"]["candidates"] == 2
    assert set(selection["channels"]) == {"rc0000", "rc0001"}
    for channel in selection["channels"].values():
        assert channel["representatives"]
        assert channel["tiers"]["active_frontier"]
        assert channel["normalization"]["method"] == "robust"
        # one candidate per channel: center == its own raw excursion, and the
        # zero-IQR fallback leaves the scale at the configured physical 1.0
        rep = channel["representatives"][0]
        raw = channel["objectives_raw"][rep]["energy_excursion"]
        assert channel["normalization"]["center"]["energy_excursion"] == raw
        assert channel["normalization"]["scale"]["energy_excursion"] == 1.0
        assert len(channel["audit"]) <= 5
    assert selection["config"]["margin_gamma"] == 0.15

    summaries = [
        json.loads(line)
        for line in (selection_dir / "candidates.jsonl").read_text().splitlines()
    ]
    assert {s["candidate_id"] for s in summaries} == {
        "pc-rc0000-t000001",
        "pc-rc0001-t000002",
    }
    for summary in summaries:
        assert summary["provenance"]["attempt_id"]
        assert summary["objectives_raw"]["energy_excursion"] is not None


def test_representative_seeds_are_resampled_chains(selection_dir):
    config = json.loads((selection_dir / "pars_selection.json").read_text(encoding="utf-8"))[
        "config"
    ]
    seeds = sorted((selection_dir / "representatives").glob("rep-*.extxyz"))
    assert seeds
    for seed in seeds:
        frames = ase_read(str(seed), index=":")
        assert len(frames) == config["refinement_images"]
        assert frames[0].info["purpose"] == "pars_refinement_seed"
        provenance = json.loads(frames[0].info["pars_provenance"])
        assert provenance["attempt_id"].startswith("t0000")


def test_select_refuses_to_overwrite(run_dir, selection_dir):
    with pytest.raises(FileExistsError):
        _select(run_dir, selection_dir)


def test_refine_runs_end_to_end(run_dir, selection_dir, tmp_path, monkeypatch):
    factory_path = tmp_path / "pars_test_factory.py"
    factory_path.write_text(FACTORY_MODULE, encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))

    exit_code = main(
        [
            "refine",
            "--selection",
            str(selection_dir),
            "--calculator-factory",
            "pars_test_factory:make_calculator",
        ]
    )
    assert exit_code == 0
    refined_dir = selection_dir / "refined"
    summary = json.loads((refined_dir / "refined_paths.json").read_text())
    assert len(summary) == 2
    for record in summary:
        # Translating three independent wells crosses three simultaneous
        # unstable modes. The CLI must persist this refusal, rather than call
        # the synchronized stationary band an index-one reaction saddle.
        assert record["status"] == "failed"
        assert record["barrier_eV"] is None
        assert record["provenance"]["saddle_check"]["negative_modes"] == 3
        chain = ase_read(
            str(refined_dir / f"refined-{record['candidate_id']}.extxyz"), index=":"
        )
        assert len(chain) >= 3
        assert chain[0].info["purpose"] == "pars_refined_path"


def test_refine_unknown_candidate_rejected(selection_dir, tmp_path, monkeypatch):
    factory_path = tmp_path / "pars_test_factory.py"
    factory_path.write_text(FACTORY_MODULE, encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))
    with pytest.raises(ValueError, match="not a representative"):
        main(
            [
                "refine",
                "--selection",
                str(selection_dir),
                "--calculator-factory",
                "pars_test_factory:make_calculator",
                "--only",
                "pc-does-not-exist",
            ]
        )


def _rewrite_jsonl(path: Path, edit) -> None:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    path.write_text("".join(json.dumps(edit(row)) + "\n" for row in rows))


def _rewrite_trial_records(run_dir: Path, edit) -> None:
    """Edit the ledger and each authoritative result.json consistently."""
    _rewrite_jsonl(run_dir / "attempts.jsonl", edit)
    for path in sorted((run_dir / "trials").glob("*/result.json")):
        path.write_text(json.dumps(edit(json.loads(path.read_text()))))


def test_select_records_skipped_attempts(selection_dir):
    selection = json.loads((selection_dir / "pars_selection.json").read_text())
    assert selection["counts"]["skipped_attempts"] == {"completed/elastic": 1}
    assert selection["counts"]["channel_members_not_ingested"] == {}


def test_select_names_channel_members_it_did_not_ingest(run_dir, tmp_path):
    # PRRS counts an unconfirmed crossing as a channel observation; PARS does not
    # take it as an A->B path, and the record has to say so
    def unconfirm(row):
        if row["attempt_id"] == "t000000":
            row.update(status="unconfirmed", classification="reactive")
        return row

    _rewrite_trial_records(run_dir, unconfirm)
    network = json.loads((run_dir / "network.json").read_text())
    network["reaction_channels"][0]["observed_by"].append("t000000")
    (run_dir / "network.json").write_text(json.dumps(network))
    out = tmp_path / "sel"
    assert _select(run_dir, out) == 0
    counts = json.loads((out / "pars_selection.json").read_text())["counts"]
    assert counts["skipped_attempts"] == {"unconfirmed/reactive": 1}
    assert counts["channel_members_not_ingested"] == {"rc0000": ["t000000"]}


def test_empty_selection_says_what_was_skipped(run_dir, tmp_path):
    # every reactive trial stopped on a saddle (the p2_seed* runs): the error
    # must say so rather than only "no candidates"
    def ridge(row):
        if row["classification"] == "reactive":
            row["status"] = "ts_candidate"
        return row

    _rewrite_trial_records(run_dir, ridge)
    with pytest.raises(ValueError, match="ts_candidate/reactive"):
        _select(run_dir, tmp_path / "sel")


def test_channel_filter_applies_before_unplaceable_traces_fail(run_dir, tmp_path):
    network = json.loads((run_dir / "network.json").read_text())
    network["reaction_channels"][1]["observed_by"] = []  # t000002 loses its channel
    (run_dir / "network.json").write_text(json.dumps(network))
    with pytest.raises(ValueError, match="no reaction channel"):
        _select(run_dir, tmp_path / "all")
    out = tmp_path / "only"
    assert _select(run_dir, out, "--channels", "rc0000") == 0
    selection = json.loads((out / "pars_selection.json").read_text())
    assert set(selection["channels"]) == {"rc0000"}


def _refine(selection_dir, *extra):
    argv = ["refine", "--selection", str(selection_dir)]
    argv += ["--calculator-factory", "pars_test_factory:make_calculator", *extra]
    return main(argv)


@pytest.fixture
def factory_on_path(tmp_path, monkeypatch):
    (tmp_path / "pars_test_factory.py").write_text(FACTORY_MODULE, encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))


A, B = "pc-rc0000-t000001", "pc-rc0001-t000002"


def _record_path(selection_dir, candidate_id):
    return selection_dir / "refined" / f"refined-{candidate_id}.json"


def test_refine_resumes_finished_candidates(selection_dir, factory_on_path):
    assert _refine(selection_dir, "--only", A) == 0
    # a marker no refinement writes: it survives only if A is not refined again
    record = json.loads(_record_path(selection_dir, A).read_text())
    record["marker"] = "kept"
    _record_path(selection_dir, A).write_text(json.dumps(record))

    assert _refine(selection_dir) == 0
    assert json.loads(_record_path(selection_dir, A).read_text())["marker"] == "kept"
    summary = json.loads((selection_dir / "refined" / "refined_paths.json").read_text())
    assert [r["candidate_id"] for r in summary] == [A, B]


def test_refine_refuses_results_from_other_conditions(selection_dir, factory_on_path):
    assert _refine(selection_dir, "--only", A) == 0
    with pytest.raises(ValueError, match="refinement_fmax_neb") as error:
        _refine(selection_dir, "--fmax-neb", "0.02")
    assert A in str(error.value)
    assert not _record_path(selection_dir, B).exists()  # nothing ran


def test_refine_refuses_when_the_seed_changed(selection_dir, factory_on_path):
    assert _refine(selection_dir, "--only", A) == 0
    seed = selection_dir / "representatives" / f"rep-{A}.extxyz"
    seed.write_text(seed.read_text() + "\n")
    with pytest.raises(ValueError, match="seed_sha256"):
        _refine(selection_dir, "--only", A)


def _rewrite_seed(selection_dir, candidate_id, edit) -> None:
    seed = selection_dir / "representatives" / f"rep-{candidate_id}.extxyz"
    frames = ase_read(str(seed), index=":")
    ase_write(str(seed), edit(frames))


def _swap_middle_species(frames):
    frames[len(frames) // 2].numbers[0] = 2  # one He image in an otherwise unchanged chain
    return frames


def _reverse_images(frames):
    return frames[::-1]


@pytest.mark.parametrize(
    ("edit", "match"),
    [(_swap_middle_species, "species differ"), (_reverse_images, "labelled")],
)
def test_refine_refuses_inconsistent_seed_images(selection_dir, factory_on_path, edit, match):
    """Audit 2026-10-11 #1: the chain was rebuilt from image 0's species, so a
    seed whose images disagree was refined as atoms the file does not record."""
    _rewrite_seed(selection_dir, A, edit)
    with pytest.raises(ValueError, match=match):
        _refine(selection_dir, "--only", A)
    assert not _record_path(selection_dir, A).exists()


def test_refine_refuses_a_seed_written_for_another_candidate(selection_dir, factory_on_path):
    other = selection_dir / "representatives" / f"rep-{B}.extxyz"
    (selection_dir / "representatives" / f"rep-{A}.extxyz").write_text(other.read_text())
    with pytest.raises(ValueError, match="seed is for candidate"):
        _refine(selection_dir, "--only", A)


def test_refine_refuses_records_without_conditions(selection_dir, factory_on_path):
    assert _refine(selection_dir, "--only", A) == 0
    record = json.loads(_record_path(selection_dir, A).read_text())
    del record["conditions"]
    _record_path(selection_dir, A).write_text(json.dumps(record))
    with pytest.raises(ValueError, match="no recorded conditions"):
        _refine(selection_dir, "--only", A)


def test_refine_redo_replaces_results(selection_dir, factory_on_path):
    assert _refine(selection_dir, "--only", A) == 0
    assert _refine(selection_dir, "--only", A, "--fmax-neb", "0.02", "--redo") == 0
    record = json.loads(_record_path(selection_dir, A).read_text())
    assert record["conditions"]["refinement"]["refinement_fmax_neb"] == 0.02


@pytest.mark.parametrize(
    "corruption", ["missing", "truncated", "legacy", "wrong_count", "wrong_id"]
)
def test_refine_refuses_unbound_or_damaged_geometry(selection_dir, factory_on_path, corruption):
    assert _refine(selection_dir, "--only", A) == 0
    marker = _record_path(selection_dir, A)
    record = json.loads(marker.read_text())
    geometry = marker.with_suffix(".extxyz")
    if corruption == "missing":
        geometry.unlink()
    elif corruption == "truncated":
        geometry.write_text("truncated extxyz\n")
    elif corruption == "legacy":
        del record["geometry_artifact"]
    elif corruption == "wrong_count":
        record["geometry_artifact"]["frames"] += 1
    else:
        record["candidate_id"] = B
    marker.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="geometry|identity"):
        _refine(selection_dir, "--only", A)
    # A damaged cohort must not enter aggregates during an unrelated refinement.
    with pytest.raises(ValueError, match="geometry|identity"):
        _refine(selection_dir, "--only", B, "--redo")
    assert not _record_path(selection_dir, B).exists()
    assert _refine(selection_dir, "--only", A, "--redo") == 0


def test_interrupted_geometry_write_preserves_committed_result(
    selection_dir, factory_on_path, monkeypatch
):
    import pars.cli as cli

    assert _refine(selection_dir, "--only", A) == 0
    marker = _record_path(selection_dir, A)
    geometry = marker.with_suffix(".extxyz")
    before = (marker.read_bytes(), geometry.read_bytes())

    def interrupted(path, frames, **kwargs):
        Path(path).write_text("partial geometry\n")
        raise OSError("injected write interruption")

    with monkeypatch.context() as patch:
        patch.setattr(cli, "ase_write", interrupted)
        with pytest.raises(OSError, match="injected"):
            _refine(selection_dir, "--only", A, "--redo")
    assert (marker.read_bytes(), geometry.read_bytes()) == before
    assert not list(geometry.parent.glob("tmp*.extxyz"))
    with monkeypatch.context() as patch:
        patch.setattr(cli.AseNebBackend, "refine", lambda *a: pytest.fail("resume recomputed"))
        assert _refine(selection_dir, "--only", A) == 0


def test_interruption_between_geometry_and_json_is_detected(
    selection_dir, factory_on_path, monkeypatch
):
    import pars.cli as cli

    assert _refine(selection_dir, "--only", A) == 0
    marker = _record_path(selection_dir, A)
    old_marker = marker.read_bytes()
    real_ase_write = cli.ase_write
    real_json_write = cli._write_json

    def new_geometry(path, frames, **kwargs):
        real_ase_write(path, frames, **kwargs)
        with Path(path).open("a") as stream:
            stream.write("\n")

    def interrupted_commit(path, payload):
        if Path(path) == marker:
            raise OSError("injected commit interruption")
        real_json_write(path, payload)

    with monkeypatch.context() as patch:
        patch.setattr(cli, "ase_write", new_geometry)
        patch.setattr(cli, "_write_json", interrupted_commit)
        with pytest.raises(OSError, match="injected"):
            _refine(selection_dir, "--only", A, "--redo")
    assert marker.read_bytes() == old_marker
    with pytest.raises(ValueError, match="changed geometry"):
        _refine(selection_dir, "--only", A)
    assert _refine(selection_dir, "--only", A, "--redo") == 0


def test_calculator_failure_without_a_band_is_persisted_and_resumable(
    selection_dir, factory_on_path, monkeypatch
):
    import pars_test_factory

    def broken_calculator(*args, **kwargs):
        raise RuntimeError("injected calculator failure")

    monkeypatch.setattr(pars_test_factory.PerAtomDoubleWell, "calculate", broken_calculator)
    assert _refine(selection_dir, "--only", A) == 0
    marker = _record_path(selection_dir, A)
    record = json.loads(marker.read_text())
    assert record["status"] == "failed"
    assert record["barrier_eV"] is None
    assert "injected calculator failure" in record["message"]
    assert record["geometry_artifact"]["frames"] == 0
    assert marker.with_suffix(".extxyz").read_bytes() == b""
    assert _refine(selection_dir, "--only", A) == 0


def test_refined_record_points_back_to_the_run(run_dir, selection_dir, factory_on_path):
    assert _refine(selection_dir, "--only", A) == 0
    record = json.loads((selection_dir / "refined" / "refined_paths.json").read_text())[0]
    assert record["provenance"]["source"]["run_dir"] == str(run_dir)
    assert record["provenance"]["source"]["attempt_id"] == "t000001"
    assert record["provenance"]["source"]["extra"]["seed_path"].endswith(
        "rep-pc-rc0000-t000001.extxyz"
    )


def test_cli_select_on_real_run_smoke():
    """Read-only smoke against a real run, skipped when runs/ is absent.

    exercises the adapter against genuine PRRS artifacts rather than the
    synthetic contract."""
    repo_root = Path(__file__).resolve().parents[2]
    candidates = sorted((repo_root / "runs").glob("connect_g5_explore_detour/*/*"))
    if not candidates:
        pytest.skip("no connect runs available")
    run = candidates[-1]
    out = run.parent / "pars_smoke_should_not_exist"
    try:
        exit_code = main(["select", "--run-dir", str(run), "--out", str(out)])
        assert exit_code == 0
        selection = json.loads((out / "pars_selection.json").read_text())
        assert selection["counts"]["channels"] >= 1
        assert selection["units"]["energy"] == "eV"
    finally:
        import shutil

        shutil.rmtree(out, ignore_errors=True)
