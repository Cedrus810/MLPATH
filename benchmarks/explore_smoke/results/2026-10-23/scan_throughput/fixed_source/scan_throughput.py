"""One resident POLAR-1: cold/warm curvature and serial/batched scan timings.

Only --dtype selects precision; no environment changes. Startup is written before
science imports, and each completed arm is saved. All scans and failures remain in
the report; one batch of m structures pays m ledger units, not one.
"""

import argparse
import cProfile
import faulthandler
import hashlib
import importlib.metadata
import io
import json
import os
import pstats
import signal
import socket
import sys
import time
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore/src"), str(ROOT)]


def plain(value):
    import numpy as np

    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def differences(left, right):
    import numpy as np

    if left["result"] is None or right["result"] is None:
        return {"unavailable": "at least one controller failed"}
    a, b = left["result"], right["result"]
    pairs = list(zip(a["points"], b["points"]))
    return {
        "statuses": [a["status"], b["status"]],
        "n_points": [len(a["points"]), len(b["points"])],
        "max_abs_dE_eV": max(
            (abs(p["energy_eV"] - q["energy_eV"]) for p, q in pairs), default=None
        ),
        "max_atom_displacement_A": max(
            (
                float(np.linalg.norm(p["positions"] - q["positions"], axis=1).max())
                for p, q in pairs
            ),
            default=None,
        ),
        "max_abs_reaction_force_diff": max(
            (abs(p["reaction_force"] - q["reaction_force"]) for p, q in pairs), default=None
        ),
        "edge_mismatch_points": [
            i for i, (p, q) in enumerate(pairs) if p["edges"] != q["edges"]
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--dtype", choices=("float32", "float64"), default="float32")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    record = {
        "host": socket.gethostname(),
        "dtype": args.dtype,
        "phase": "import_runtime",
        "arms": {},
        "complete": False,
        "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS")},
    }

    def save():
        # Before importing NumPy this record only contains ordinary JSON types.
        value = plain(record) if "numpy" in sys.modules else record
        temp = args.out / "measurement.json.part"
        temp.write_text(json.dumps(value, indent=2) + "\n")
        temp.replace(args.out / "measurement.json")

    def interrupted(signum, frame):
        raise TimeoutError(f"measurement interrupted by signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    faulthandler.enable()
    faulthandler.dump_traceback_later(60, repeat=True)
    save()
    print("phase=import_runtime", flush=True)
    started = time.perf_counter()
    try:
        import numpy as np
        import torch
        from ase.io import read

        from benchmarks.b04_sn2_chloride.calculator import MODELS, polar1
        from explore.ledger import Ledger
        from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many
        from explore.probe import Intent, ScanSettings
        from prrs.config import SearchConfig
        from prrs.search import _stamped_with_charge

        record["import_s"] = time.perf_counter() - started
        record["phase"] = "load_model"
        paths = [
            Path(__file__),
            ROOT / "explore/src/explore/probe.py",
            ROOT / "explore/src/explore/parallel.py",
            ROOT / "explore/src/explore/ledger.py",
            ROOT / "src/prrs/reliability.py",
        ]
        record["source_sha256"] = {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
        }
        record["versions"] = {
            k: importlib.metadata.version(k) for k in ("torch", "mace-torch", "ase", "numpy")
        }
        record["model"] = MODELS["polar1"]
        save()
        print("phase=load_model", flush=True)
        loaded = time.perf_counter()
        calc = polar1("cuda", args.dtype)
        record["load_s"] = time.perf_counter() - loaded
        record["gpu"] = torch.cuda.get_device_name()
        backend = MACEBatchEvaluator(calc)
        config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
        base = _stamped_with_charge(
            read(ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"), config
        )
        cation = base.copy()
        cation.info["charge"] = 1  # explicit diagnostic electronic-state variant
        cation = _stamped_with_charge(
            cation, SearchConfig(charge_sensitive=True, total_charge=1, multiplicity=1)
        )
        large = _stamped_with_charge(
            read(ROOT / "benchmarks/carbene_ambimodal/probe_b09/input.extxyz"),
            SearchConfig(charge_sensitive=True, total_charge=0, multiplicity=1),
        )
        # Same geometry, different electronic state, plus a different atom count.
        frames = [base, cation, large]
        record["parity_inputs"] = [
            {"n_atoms": len(a), "charge": a.info["charge"], "spin": a.info["spin"]}
            for a in frames
        ]
        reference = []
        for a in frames:
            calc.reset()
            a.calc = calc
            reference.append((a.get_potential_energy(), a.get_forces(apply_constraint=False)))
        record["parity"] = []
        for order in ((0,), (0, 1, 2), (2, 0, 1)):
            values = backend([frames[i] for i in order])
            record["parity"].append(
                {
                    "order": order,
                    "max_abs_dE_eV": [
                        abs(e - reference[i][0]) for i, (e, f) in zip(order, values)
                    ],
                    "max_abs_dF_eV_A": [
                        np.max(abs(f - reference[i][1])) for i, (e, f) in zip(order, values)
                    ],
                }
            )
        record["phase"] = "scan_arms"
        save()
        print("phase=scan_arms", flush=True)

        def arm(name, reuse, count, active):
            jobs = []
            for i in range(count):
                atoms = base.copy()
                atoms.rotate(17 * i, "z", center="COP")
                settings = ScanSettings(
                    tuple(round(-0.1 * k, 10) for k in range(18)),
                    0.01,
                    300,
                    3.0,
                    reuse_curvature=reuse,
                )
                jobs.append(
                    ScanJob(
                        atoms,
                        Intent(f"scan{i}", ((0, 5),), (1.0,)),
                        settings,
                        config,
                        Ledger(args.out / f"{name}_{i}.jsonl", name, budget=500),
                        str(i),
                        reserve=500,
                    )
                )
            first_call = len(backend.calls)
            profile = cProfile.Profile()
            profile.enable()
            try:
                outcomes, metrics = scan_many(jobs, backend, max_active=active)
            finally:
                profile.disable()
                buf = io.StringIO()
                pstats.Stats(profile, stream=buf).sort_stats("cumulative").print_stats(30)
                (args.out / f"{name}_controller_profile.txt").write_text(buf.getvalue())
                for job in jobs:
                    job.ledger.close()
            record["arms"][name] = {
                "reuse_curvature": reuse,
                "active": active,
                "metrics": metrics,
                "backend_calls": backend.calls[first_call:],
                "jobs": [
                    {
                        "spent": job.ledger.spent,
                        "error": outcome.error,
                        "result": asdict(outcome.result)
                        if outcome.result is not None
                        else None,
                    }
                    for job, outcome in zip(jobs, outcomes)
                ],
            }
            save()
            print(
                json.dumps(
                    {
                        "arm": name,
                        "wall_s": metrics["wall_s"],
                        "spent": [j.ledger.spent for j in jobs],
                    }
                ),
                flush=True,
            )

        arm("cold", False, 1, 1)
        arm("serial_warm", True, 4, 1)
        arm("batch_warm", True, 4, 4)
        arms = record["arms"]
        record["cold_vs_warm"] = differences(
            arms["cold"]["jobs"][0], arms["serial_warm"]["jobs"][0]
        )
        record["serial_vs_batch"] = [
            differences(a, b)
            for a, b in zip(arms["serial_warm"]["jobs"], arms["batch_warm"]["jobs"])
        ]
        record["batch_speedup"] = (
            arms["serial_warm"]["metrics"]["wall_s"] / arms["batch_warm"]["metrics"]["wall_s"]
        )
        record["complete"] = True
        record["phase"] = "complete"
    except BaseException as exc:
        record["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        record["total_wall_s"] = time.perf_counter() - started
        save()
        faulthandler.cancel_dump_traceback_later()


if __name__ == "__main__":
    main()
