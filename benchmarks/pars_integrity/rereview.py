"""Reproduce remaining PARS electronic-state, Hessian and ledger audit gaps.

Uses analytic external fields and synthetic PRRS files in a fresh temporary
directory. Does not alter production runs or implementation sources. These
assertions pin observed defects, not the desired repaired behavior.
"""

# ruff: noqa: E402 -- standalone package and test fixture bootstrap

import argparse
import subprocess
import json
import sys
import tempfile
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "pars/src"), str(ROOT / "src"), str(ROOT / "pars/tests")]
from conftest import build_run, make_candidate
from pars.cli import main
from pars.adapters.prrs_run import PrrsRun
from pars.config import ParsConfig
from pars.paths.ingest import ingest_run
from pars.refinement.ase_neb import AseNebBackend

out = Path(tempfile.mkdtemp(prefix="mlpath-rereview-"))
args = argparse.ArgumentParser(description=__doc__)
args.add_argument("--output", type=Path, required=True)
args = args.parse_args()
from pars.provenance import implementation_sha256, numerical_versions

report = {
    "workspace": str(ROOT),
    "evidence_directory": str(out),
    "inspected_commit": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip(),
    "pars_implementation_sha256": implementation_sha256(),
    "numerical_versions": numerical_versions(),
}

# Full charged source -> representative -> backend -> committed CLI result.
run = build_run(out / "charge-source")
network = json.loads((run / "network.json").read_text())
network["reaction_channels"] = network["reaction_channels"][:1]
network["reactions"] = network["reactions"][:1]
(run / "network.json").write_text(json.dumps(network))
ledger = [json.loads(line) for line in (run / "attempts.jsonl").read_text().splitlines()]
ledger = [row for row in ledger if row["attempt_id"] == "t000001"]
(run / "attempts.jsonl").write_text("".join(json.dumps(row) + "\n" for row in ledger))
manifest = json.loads((run / "config.json").read_text())
manifest["config"] = {"charge_sensitive": True, "total_charge": -1, "multiplicity": 3}
(run / "config.json").write_text(json.dumps(manifest))
trial = run / "trials/t000001"
phases = ["source", "perturbed", "free_start", "free", "free", "free", "quench", "quenched"]
xs = [0, 0.1, 0.2, 0.7, 1.2, 1.8, 2, 2]
frames, rows = [], []
for i, (phase, x) in enumerate(zip(phases, xs)):
    atoms = Atoms("H", positions=[[x, 0, 0]], info={"phase": phase, "charge": -1, "spin": 3})
    energy = 0.5 * x * x * (x - 2) ** 2
    forces = np.array([[-2 * x * (x - 2) * (x - 1), 0, 0]])
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    frames.append(atoms)
    rows.append({"frame": i, "phase": phase, "potential_eV": energy, "edges": []})
write(str(trial / "response.traj"), frames)
write(str(trial / "endpoint.extxyz"), frames[-1])
(trial / "observations.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
factory_source = """import numpy as np
from ase.calculators.calculator import Calculator, all_changes
SEEN = []
class ChargedWell(Calculator):
    implemented_properties = ["energy", "forces"]
    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        q = self.atoms.info.get("charge", 0)
        SEEN.append({"charge": self.atoms.info.get("charge"), "spin": self.atoms.info.get("spin")})
        scale = 1 + abs(q)
        x, y, z = self.atoms.positions[0]
        self.results = {"energy": scale*(0.25*x*x*(x-2)**2 + 0.5*(y*y+z*z)), "forces": scale*np.array([[-x*(x-2)*(x-1), -y, -z]])}
def make():
    return ChargedWell()
"""
(out / "charged_review_factory.py").write_text(factory_source)
sys.path.insert(0, str(out))
selection = out / "charged-selection"
config_path = out / "pars-config.json"
config_path.write_text(
    json.dumps({"refinement_remove_rigid_motion": False, "refinement_fmax_neb": 0.001})
)
main(["select", "--run-dir", str(run), "--out", str(selection), "--config", str(config_path)])
main(
    [
        "refine",
        "--selection",
        str(selection),
        "--calculator-factory",
        "charged_review_factory:make",
    ]
)
import charged_review_factory as charged

record = json.loads((selection / "refined/refined-pc-rc0000-t000001.json").read_text())
seed = read(str(selection / "representatives/rep-pc-rc0000-t000001.extxyz"))
report["charge_spin_loss"] = {
    "source_charge": frames[0].info["charge"],
    "source_spin": frames[0].info["spin"],
    "seed_charge": seed.info.get("charge"),
    "seed_spin": seed.info.get("spin"),
    "backend_observed_states": sorted(
        {(v["charge"], v["spin"]) for v in charged.SEEN}, key=str
    ),
    "status": record["status"],
    "reported_barrier_eV": record["barrier_eV"],
    "correct_charged_barrier_eV": 0.5,
    "connectivity_status": record["provenance"]["connectivity_check"]["status"],
}
assert record["status"] == "ci_neb_converged"
assert abs(record["barrier_eV"] - 0.25) < 1e-5
assert all(v == {"charge": None, "spin": None} for v in charged.SEEN)


# Uncertainty is invalid only at transverse force-difference Hessian samples.
class HessianUncertainWell(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, bad):
        super().__init__()
        self.bad = bad
        self.bad_samples = 0

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        transverse = abs(y) > 1e-8 or abs(z) > 1e-8
        self.bad_samples += int(transverse)
        self.results = {
            "energy": 0.25 * x * x * (x - 2) ** 2 + 0.5 * (y * y + z * z),
            "forces": np.array([[-x * (x - 2) * (x - 1), -y, -z]]),
            "force_uncertainty_eV_A": self.bad if transverse else 0.005,
        }


positions = np.zeros((9, 1, 3))
positions[:, 0, 0] = np.linspace(0, 2, 9)
candidate = make_candidate(
    "uncertain-hessian", "channel", np.zeros(9), [frozenset()] * 9, positions=positions
)
config = ParsConfig(
    refinement_remove_rigid_motion=False,
    refinement_fmax_neb=0.001,
    refinement_uncertainty_refusal=0.01,
)
report["hessian_uncertainty_bypass"] = []
for value in [float("nan"), 1000.0]:
    calc = HessianUncertainWell(value)
    result = AseNebBackend(calc).refine(candidate, config)
    report["hessian_uncertainty_bypass"].append(
        {
            "bad_value": str(value),
            "bad_hessian_samples": calc.bad_samples,
            "status": result.status.value,
            "barrier_eV": result.barrier_eV,
            "recorded_uncertainty_max_eV_A": result.uncertainty_max_eV_A,
            "coverage": result.provenance.get("uncertainty_coverage"),
            "connectivity": result.provenance.get("connectivity_check", {}).get("status"),
        }
    )
    assert calc.bad_samples > 0 and result.status.value == "ci_neb_converged"

# Conflict with result.json bypasses validation when the ledger excludes the trial.
run2 = build_run(out / "ledger-source")
path = run2 / "attempts.jsonl"
rows = [json.loads(line) for line in path.read_text().splitlines()]
for row in rows:
    if row["attempt_id"] == "t000001":
        row["classification"] = "elastic"
path.write_text("".join(json.dumps(row) + "\n" for row in rows))
traces = ingest_run(run2, ParsConfig())
report["skipped_authoritative_conflict"] = {
    "ledger_classification": "elastic",
    "result_classification": "reactive",
    "ingest_succeeded": True,
    "ingested_ids": [t.trace_id for t in traces],
}
assert "t000001" not in report["skipped_authoritative_conflict"]["ingested_ids"]
try:
    PrrsRun(run2).load_trace("t000001")
except ValueError as exc:
    report["skipped_authoritative_conflict"]["explicit_load_error"] = str(exc)
else:
    raise AssertionError("expected direct load rejection")
(out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
print(json.dumps(report, indent=2, allow_nan=False))
