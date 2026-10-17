"""Smoke: one constrained-relaxation scan on b04 (Cl- + CH3Cl, POLAR-1, charge -1).

Measures per-point forwards and wall time, the energy noise of a relaxed point, and
whether the probe reads a ridge. Not an acceptance run; the known TS is only compared
after the fact.
"""

import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, counted_factory  # noqa: E402
from explore.probe import Intent, ScanSettings, constrained_scan, read_scan, release  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.io import snapshot  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402
from prrs.state import encode  # noqa: E402

NOISE_METHOD = "perturb 0.02 A then pull q back to mid.q_target; energy of second point"
NOISE_DEFINITION = (
    "range (ptp) of converged repeat energies at mid point; used directly as read_scan noise_eV"
)


def _point_record(point):
    return {
        k: v
        for k, v in asdict(point).items()
        if k not in ("positions", "forces_pes")
    }


def _write_scan(out, record):
    (out / "scan.json").write_text(json.dumps(record, indent=2, default=str))


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
    atoms = _stamped_with_charge(read(ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"), config)
    graph = encode(atoms, config.bond_scale)
    bonded = {i for e in graph.edges for i in e}
    nucleophile = next(i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Cl" and i not in bonded)
    carbon = next(i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "C")
    ledger = Ledger(out / "ledger.jsonl", arm="smoke", budget=20000)
    factory = counted_factory(polar1, ledger)
    q0 = float(np.linalg.norm(atoms.positions[nucleophile] - atoms.positions[carbon]))
    steps = int(round((q0 - 1.9) / 0.1))
    settings = ScanSettings(
        deltas_A=tuple(round(-0.1 * k, 10) for k in range(steps + 1)),
        fmax_eV_A=0.01,
        max_steps_per_point=300,
        max_rise_eV=3.0,
    )
    intent = Intent(f"cmp_{nucleophile}_{carbon}", ((nucleophile, carbon),), (1.0,))
    result = constrained_scan(atoms, intent, settings, factory, config, ledger, task="scan")

    # Written before anything else, so a later crash keeps the scan.
    record = {
        "status": result.status,
        "q0": result.q0,
        "points": [_point_record(p) for p in result.points],
        "triggers": None,
        "release": None,
        "noise_status": None,
        "noise_eV_used": None,
        "noise_definition": NOISE_DEFINITION,
        "spent": ledger.spent,
        "error": None,
    }
    if not result.points:
        record["skipped"] = "empty scan: no points; noise, read_scan and release not run"
        ledger.close()
        record["spent"] = ledger.spent
        _write_scan(out, record)
        return
    _write_scan(out, record)

    triggers = None
    released = None
    error = None
    try:
        # noise: perturb the mid-scan point by N(0, 0.02 A), pull q back to mid.q_target,
        # relax, and take the energy of the second point (index 0 is the perturbed start
        # itself). 3 repeats; only converged repeats enter the range.
        mid = result.points[len(result.points) // 2]
        rng = np.random.default_rng(0)
        runs = []
        for r in range(3):
            start = atoms.copy()
            start.positions = mid.positions + rng.normal(scale=0.02, size=mid.positions.shape)
            q_perturbed = intent.value(start.positions)
            rep = constrained_scan(
                start,
                intent,
                ScanSettings(
                    (0.0, mid.q_target - q_perturbed),
                    settings.fmax_eV_A,
                    settings.max_steps_per_point,
                    settings.max_rise_eV,
                ),
                factory,
                config,
                ledger,
                task=f"noise{r}",
            )
            run = {
                "repeat": r,
                "status": rep.status,
                "q_perturbed": q_perturbed,
                "q_target": None,
                "q_actual": None,
                "converged": None,
                "steps": None,
                "energy_eV": None,
                "in_range": False,
            }
            if len(rep.points) >= 2:
                p = rep.points[1]
                run.update(
                    q_target=p.q_target,
                    q_actual=p.q_actual,
                    converged=p.converged,
                    steps=p.steps,
                    energy_eV=p.energy_eV,
                    in_range=rep.status == "completed" and bool(p.converged),
                )
            runs.append(run)
        in_range = [run["energy_eV"] for run in runs if run["in_range"]]
        noise_eV_used = float(np.ptp(in_range)) if len(in_range) >= 2 else None
        record["noise_eV_used"] = noise_eV_used
        record["noise_status"] = "measured" if noise_eV_used is not None else "noise_unmeasured"
        noise = {
            "method": NOISE_METHOD,
            "mid_delta_A": mid.delta_A,
            "reference_eV": mid.energy_eV,
            "runs": runs,
            "noise_status": record["noise_status"],
            "noise_eV_used": noise_eV_used,
            "noise_definition": NOISE_DEFINITION,
        }
        (out / "noise.json").write_text(json.dumps(noise, indent=2, default=str))

        if noise_eV_used is not None:
            triggers = read_scan(result.points, noise_eV=noise_eV_used)
        ridge = next((t for t in triggers or [] if t["kind"] == "ridge_suspected"), None)
        if ridge:
            source = snapshot(atoms, GuardedCalculator(factory(), config))
            rel = release(
                result.points[ridge["release"]].positions,
                atoms,
                source,
                factory,
                config,
                ledger,
                task="release",
                seed=0,
            )
            released = {
                k: rel[k] for k in ("status", "label", "fingerprint", "seed", "constraints_removed")
            }
    except BaseException as exc:
        error = repr(exc)
        raise
    finally:
        ledger.close()
        record.update(triggers=triggers, release=released, spent=ledger.spent, error=error)
        _write_scan(out, record)


if __name__ == "__main__":
    main(sys.argv[1])
