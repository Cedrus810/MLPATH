"""Synthetic audit counterexamples; assertions describe defects, not acceptance.

Writes only under the supplied fresh evidence directory. Never changes source
or production runs. Run using openmm_dev, without a neural model or GPU.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / p) for p in ("src", "pars/src", "mature/src")]

import numpy as np  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402
from ase.io import write  # noqa: E402
from mature.io import OutputStore  # noqa: E402
from mature.schema import ShotRecord, sha256_object  # noqa: E402
from pars.adapters.prrs_run import PrrsRun  # noqa: E402
from pars.cli import _candidate_from_seed  # noqa: E402
from pars.config import ParsConfig  # noqa: E402
from pars.paths.ingest import ingest_run, skipped_attempts  # noqa: E402
from pars.provenance import implementation_sha256  # noqa: E402
from pars.refinement.ase_neb import AseNebBackend  # noqa: E402
from pars.schema import Provenance  # noqa: E402
from prrs.provenance import calculator_identity  # noqa: E402


def fixture_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Quartic(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        self.results = {
            "energy": float((x * x - 1) ** 2 + 0.5 * (y * y + z * z)),
            "forces": np.array([[-4 * x * (x * x - 1), -y, -z]]),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    pars_fixtures = fixture_module("audit_pars_fixtures", "pars/tests/conftest.py")
    mature_fixtures = fixture_module("audit_mature_fixtures", "mature/tests/conftest.py")
    report = {
        "workspace": str(ROOT),
        "pars_implementation_sha256": implementation_sha256(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "evidence_scope": "synthetic CPU counterexamples; no historical contamination claim",
    }

    # A well-formed extxyz chain changes species in one interior image. The
    # reader silently rebuilds it as H on every image and certifies that PES.
    provenance = Provenance(
        run_dir="synthetic",
        attempt_id="t000001",
        channel_id="rc0000",
        event_key=None,
        reactant_state_id=None,
        product_state_id=None,
        frame_start=0,
        frame_stop=7,
        electronic_state={"charge": 0, "spin": 1},
    )
    chain = []
    for index, x in enumerate(np.linspace(-1, 1, 7)):
        chain.append(
            Atoms(
                "He" if index == 3 else "H",
                positions=[[x, 0, 0]],
                info={
                    "candidate_id": "pc-species",
                    "channel_id": "rc0000",
                    "image": index,
                    "charge": 0,
                    "spin": 1,
                    "pars_provenance": json.dumps(provenance.to_dict()),
                },
            )
        )
    seed = out / "rep-pc-species.extxyz"
    write(seed, chain)
    candidate = _candidate_from_seed(seed)
    refined = AseNebBackend(Quartic()).refine(
        candidate, ParsConfig(refinement_remove_rigid_motion=False, refinement_fmax_neb=0.001)
    )
    assert candidate.numbers.tolist() == [1]
    assert refined.status.value == "ci_neb_converged", refined.message
    assert refined.provenance["saddle_check"]["connectivity"] == "matched"
    report["seed_species_mismatch"] = {
        "input_species": [a.get_chemical_symbols() for a in chain],
        "adopted_species": candidate.numbers.tolist(),
        "status": refined.status.value,
        "barrier_eV": refined.barrier_eV,
        "connectivity": refined.provenance["saddle_check"]["connectivity"],
    }

    # An earlier contradictory duplicate is replaced by the final row in the
    # adapter's ID map, then checked_attempts validates the final row twice.
    run_path = pars_fixtures.build_run(out / "duplicate-ledger")
    ledger_path = run_path / "attempts.jsonl"
    rows = [json.loads(line) for line in ledger_path.read_text().splitlines()]
    original = next(row for row in rows if row["attempt_id"] == "t000001")
    contradictory = {**original, "status": "failed", "classification": "elastic"}
    rows.insert(1, contradictory)
    ledger_path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    run = PrrsRun(run_path)
    checked = run.checked_attempts()
    traces = ingest_run(run, ParsConfig())
    assert [row["classification"] for row in checked if row["attempt_id"] == "t000001"] == [
        "reactive",
        "reactive",
    ]
    assert [t.trace_id for t in traces].count("t000001") == 2
    report["duplicate_attempt_id"] = {
        "ledger_classifications": [
            r["classification"] for r in rows if r["attempt_id"] == "t000001"
        ],
        "checked_classifications": [
            r["classification"] for r in checked if r["attempt_id"] == "t000001"
        ],
        "ingested_trace_ids": [t.trace_id for t in traces],
        "skipped": skipped_attempts(run, ParsConfig()),
    }

    def batch(name):
        source = mature_fixtures.make_source_run(out / name / "source")
        protocol = mature_fixtures.make_protocol(source["root"], source["files"])
        folder = out / name / "batch"
        implementation = mature_fixtures.fake_sha("actual-implementation")
        store = OutputStore.initialize(folder, protocol, implementation)
        proposal = mature_fixtures.make_proposal(folder, "p1")
        store.add_proposal(proposal)
        tasks = store.build_task_table(protocol, [proposal])
        record = ShotRecord(
            proposal_id="p1",
            shot_id="s0000",
            batch_role="explore",
            seed=tasks[0]["seed"],
            protocol_sha256=protocol.sha256,
            implementation_sha256=implementation,
            outcome="A",
            entry_side="A",
            entry_step=200,
            trial_dir="shots/p1/s0000",
        )
        return store, protocol, tasks, record

    # A different protocol, implementation, and seed still complete a task
    # allocated by the current batch. Statistics contain no cohort check.
    store, protocol, tasks, record = batch("foreign-record")
    wrong = dataclasses.replace(
        record,
        protocol_sha256=mature_fixtures.fake_sha("other-protocol"),
        implementation_sha256=mature_fixtures.fake_sha("other-implementation"),
        seed=999,
    )
    store.commit(wrong)
    estimate = store.estimates()["p1"]
    assert "p1#s0000" in store.done_task_ids()
    assert estimate.shots == 1
    report["foreign_shot_conditions"] = {
        "allocated_seed": tasks[0]["seed"],
        "accepted_seed": wrong.seed,
        "protocol_matches": wrong.protocol_sha256 == protocol.sha256,
        "implementation_matches": wrong.implementation_sha256 == record.implementation_sha256,
        "counted_shots": estimate.shots,
        "pending": len(store.pending_tasks()),
    }

    # A marker's payload changes A to B without changing its stored SHA-256.
    # Resume, estimates, and idempotent commit all accept the damaged marker.
    store, protocol, tasks, record = batch("damaged-marker")
    store.commit(record)
    marker_path = store.output_dir / "commits/p1/s0000.json"
    marker = json.loads(marker_path.read_text())
    marker["record"].update(outcome="B", entry_side="B")
    marker_path.write_text(json.dumps(marker))
    assert marker["record_sha256"] != sha256_object(marker["record"])
    resumed = OutputStore.resume(store.output_dir, protocol, record.implementation_sha256)
    resumed.commit(record)
    estimate = resumed.estimates()["p1"]
    assert estimate.p_B_given_AB == 1.0
    report["unchecked_commit_hash"] = {
        "hash_matches": marker["record_sha256"] == sha256_object(marker["record"]),
        "result_json_outcome": json.loads(
            (store.output_dir / "shots/p1/s0000/result.json").read_text()
        )["outcome"],
        "marker_outcome": marker["record"]["outcome"],
        "p_B_given_AB": estimate.p_B_given_AB,
        "resume_accepted": True,
        "recommit_accepted": True,
    }

    # Reconfirm the documented, deliberately deferred bound-method gap.
    class Factory:
        def __init__(self, scale):
            self.scale = scale

        def build(self):
            return self.scale

    one, two = Factory(1), Factory(2)
    id_one, id_two = calculator_identity(one.build, {}), calculator_identity(two.build, {})
    assert id_one == id_two and id_one["verifiable"]
    report["known_deferred_bound_method_identity"] = {
        "same_identity": id_one == id_two,
        "verifiable": id_one["verifiable"],
        "factory_outputs": [one.build(), two.build()],
    }
    (out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
