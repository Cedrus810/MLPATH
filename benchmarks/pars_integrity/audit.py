"""Reproduce endpoint identity, refinement durability and import alignment gaps.

Analytic potentials and temporary synthetic PRRS files only. No production run
is modified. The import/durability cases reuse the repository's test fixture;
their recorded refinement is a deliberate index-three refusal, not CI success.
"""

from __future__ import annotations

# ruff: noqa: E402 -- bootstrap the independently packaged module and fixture paths

import argparse
import contextlib
import io
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pars/src"))
sys.path.insert(0, str(ROOT / "pars/tests"))

import pars.cli as cli
from conftest import build_run
from test_cli import FACTORY_MODULE, _select
from pars.adapters.prrs_run import PrrsRun
from pars.config import ParsConfig
from pars.provenance import implementation_sha256, numerical_versions
from pars.refinement.ase_neb import AseNebBackend
from pars.schema import PathCandidate, Provenance


def endpoint_identity():
    n = 500
    fixed = np.zeros((n, 3))
    fixed[:, 0] = 3 * np.arange(n)

    class ThreeWells(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            d = self.atoms.positions - fixed
            x = d[0, 0]
            forces = -d.copy()
            forces[0, 0] = -2 * x * (x * x - 1) * (3 * x * x - 1)
            energy = x * x * (x * x - 1) ** 2
            energy += 0.5 * (np.sum(d[1:] ** 2) + np.sum(d[0, 1:] ** 2))
            self.results = {"energy": float(energy), "forces": forces}

    xs = np.linspace(-1, 1, 9)
    positions = np.repeat(fixed[None], len(xs), axis=0)
    positions[:, 0, 0] = xs
    candidate = PathCandidate(
        candidate_id="spectator-dilution",
        channel_id="analytic",
        trace_id="analytic",
        numbers=np.ones(n, dtype=int),
        positions=positions,
        energies=xs**2 * (xs**2 - 1) ** 2,
        uncertainty=None,
        signatures=None,
        provenance=Provenance(
            run_dir="analytic",
            attempt_id="analytic",
            channel_id="analytic",
            event_key=None,
            reactant_state_id="minus-one",
            product_state_id="plus-one",
            frame_start=0,
            frame_stop=len(xs),
        ),
    )
    config = ParsConfig(refinement_remove_rigid_motion=False, refinement_fmax_neb=0.001)
    result = AseNebBackend(ThreeWells()).refine(candidate, config)
    connection = result.provenance.get("connectivity_check", {})

    def branches(report):
        return [
            {
                "sign": b["sign"],
                "active_x": b["positions_A"][0][0],
                "matched_endpoint": b.get("matched_endpoint"),
                "endpoint_rmsd_A": b.get("endpoint_rmsd_A"),
                "endpoint_max_displacement_A": b.get("endpoint_max_displacement_A"),
                "active_atom_error_to_assigned_endpoint_A": abs(
                    b["positions_A"][0][0]
                    - (1 if b.get("matched_endpoint") == "product" else -1)
                )
                if b.get("matched_endpoint")
                else None,
            }
            for b in report.get("branches", [])
        ]

    fine = connection.get("resolution_check", {})
    return {
        "atoms": n,
        "potential": "x^2*(x^2-1)^2 + harmonic spectators/transverse DOFs",
        "true_minima_x": [-1, 0, 1],
        "outer_endpoint_rmsd_A": 2 / np.sqrt(n),
        "config": config.to_dict(),
        "status": result.status.value,
        "message": result.message,
        "barrier_eV": result.barrier_eV,
        "connectivity_status": connection.get("status"),
        "saddle_negative_modes": result.provenance.get("saddle_check", {}).get(
            "negative_modes"
        ),
        "branches": branches(connection),
        "resolution_status": fine.get("status"),
        "fine_branches": branches(fine.get("fine", {})),
    }


def file_integrity():
    def resume(args):
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                return {"exit_code": cli.main(args), "refused": False}
        except ValueError as error:
            return {"exit_code": None, "refused": True, "reason": str(error).split(": ", 1)[-1]}

    with tempfile.TemporaryDirectory(prefix="pars-integrity-audit-") as temp:
        root = Path(temp)
        run = build_run(root)
        out = root / "selection"
        (root / "audit_factory.py").write_text(FACTORY_MODULE)
        sys.path.insert(0, temp)
        args = [
            "refine",
            "--selection",
            str(out),
            "--calculator-factory",
            "audit_factory:make_calculator",
            "--only",
            "pc-rc0000-t000001",
        ]
        with contextlib.redirect_stdout(io.StringIO()):
            assert _select(run, out) == 0
            assert cli.main(args) == 0
        geometry = out / "refined/refined-pc-rc0000-t000001.extxyz"
        marker = geometry.with_suffix(".json")
        old_marker = marker.read_bytes()
        initial_status = json.loads(old_marker)["status"]
        geometry.unlink()
        missing = {**resume(args), "geometry_exists": geometry.exists()}
        # Restore the committed pair before independently testing redo interruption.
        with contextlib.redirect_stdout(io.StringIO()):
            assert cli.main([*args, "--redo"]) == 0
        before_geometry = geometry.read_bytes()
        real_write = cli.ase_write

        def interrupted_write(path, frames, **kwargs):
            Path(path).write_text("truncated extxyz\n")
            raise RuntimeError("simulated interruption after opening geometry")

        cli.ase_write = interrupted_write
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                cli.main([*args, "--redo"])
        except RuntimeError:
            pass
        finally:
            cli.ase_write = real_write
        crash = {
            **resume(args),
            "old_marker_preserved": marker.read_bytes() == old_marker,
            "old_geometry_preserved": geometry.read_bytes() == before_geometry,
        }
        observations = run / "trials/t000001/observations.jsonl"
        rows = [json.loads(line) for line in observations.read_text().splitlines()]
        rows.reverse()
        observations.write_text("".join(json.dumps(row) + "\n" for row in rows))
        reordered = {
            "observation_frame_indices": [r["frame"] for r in rows],
            "expected_first_phase": "source",
            "expected_first_energy_eV": -1.0,
        }
        try:
            trace = PrrsRun(run).load_trace("t000001")
            reordered.update(
                accepted=True,
                first_phase=trace.phases[0],
                first_energy_eV=float(trace.energies[0]),
            )
        except ValueError as error:
            reordered.update(accepted=False, reason=str(error))
        return {
            "initial_record_status": initial_status,
            "missing_geometry_resume": missing,
            "interrupted_redo_resume": crash,
            "reordered_observations": reordered,
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "inspected_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "pars_implementation_sha256": implementation_sha256(),
        "numerical_versions": numerical_versions(),
        "endpoint_identity": endpoint_identity(),
        "file_integrity": file_integrity(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
