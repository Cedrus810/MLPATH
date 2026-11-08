"""Bounded POLAR-1 scan performance probe; no production calculator changes.

Compare the installed calculator with a benchmark-only forward wrapper disabling
stress. One dtype per process. Fixed geometries are evaluated in alternating arm
order; cache is cleared on every evaluation. Optional b04 scans use the exact
production scan, settings and ledger, and report every point without filtering.
"""

import argparse
import contextlib
import faulthandler
import hashlib
import importlib.metadata
import json
import os
import platform
import socket
import statistics
import sys
import time
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore/src"), str(ROOT)]

# Make cold import/model-load stalls visible before the first result is written.
faulthandler.enable()
faulthandler.dump_traceback_later(45, repeat=True)
print("phase=import_science_runtime", flush=True)

import numpy as np  # noqa: E402
import torch  # noqa: E402
from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import MODELS, polar1  # noqa: E402
from explore.ledger import Ledger, reused_counted_factory  # noqa: E402
from explore.probe import Intent, ScanSettings, constrained_scan  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402


@contextlib.contextmanager
def stress_mode(calc, enabled):
    """Only change the forward keyword; preserve ASE batching/results/units."""
    originals = [model.forward for model in calc.models]
    try:
        if not enabled:
            for model, original in zip(calc.models, originals):

                def forward(*args, _original=original, **kwargs):
                    kwargs["compute_stress"] = False
                    return _original(*args, **kwargs)

                model.forward = forward
        calc.reset()
        yield
    finally:
        for model, original in zip(calc.models, originals):
            model.forward = original
        calc.reset()


def evaluate(calc, atoms):
    calc.reset()
    atoms.calc = calc
    torch.cuda.synchronize()
    start = time.perf_counter()
    energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces()
    torch.cuda.synchronize()
    return time.perf_counter() - start, energy, forces.copy()


def plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def write(out, record):
    (out / "measurement.json").write_text(json.dumps(plain(record), indent=2) + "\n")


def fixed_geometries(calc, atoms, repeats, out, record, name):
    samples = {"stress_on": [], "stress_off": []}
    differences = []
    # Warm both paths before timing; same seeded perturbations for both arms.
    for enabled in (True, False):
        with stress_mode(calc, enabled):
            evaluate(calc, atoms.copy())
    rng = np.random.default_rng(0)
    for k in range(repeats):
        frame = atoms.copy()
        frame.positions += rng.normal(scale=0.01, size=frame.positions.shape)
        values = {}
        for enabled in (True, False) if k % 2 == 0 else (False, True):
            label = "stress_on" if enabled else "stress_off"
            with stress_mode(calc, enabled):
                wall, energy, forces = evaluate(calc, frame.copy())
            samples[label].append(wall)
            values[label] = (energy, forces)
        on_e, on_f = values["stress_on"]
        off_e, off_f = values["stress_off"]
        differences.append(
            {"abs_dE_eV": abs(off_e - on_e), "max_abs_dF_eV_A": np.max(abs(off_f - on_f))}
        )
    medians = {k: statistics.median(v) for k, v in samples.items()}
    record["fixed_geometries"][name] = {
        "n_atoms": len(atoms),
        "charge": atoms.info["charge"],
        "wall_samples_s": samples,
        "median_s": medians,
        "speedup": medians["stress_on"] / medians["stress_off"],
        "differences": differences,
    }
    write(out, record)


def scan(calc, atoms, enabled, repeat, out):
    label = f"stress_{'on' if enabled else 'off'}_{repeat}"
    ledger = Ledger(out / f"{label}.jsonl", arm=label, budget=400)
    factory = reused_counted_factory(lambda: calc, ledger)
    settings = ScanSettings(
        tuple(round(-0.1 * k, 10) for k in range(18)),
        fmax_eV_A=0.01,
        max_steps_per_point=300,
        max_rise_eV=3.0,
    )
    try:
        with stress_mode(calc, enabled):
            torch.cuda.synchronize()
            start = time.perf_counter()
            result = constrained_scan(
                atoms,
                Intent("b04_C0_Cl5", ((0, 5),), (1.0,)),
                settings,
                factory,
                SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1),
                ledger,
                task=label,
                reserve=400,
            )
            torch.cuda.synchronize()
            wall = time.perf_counter() - start
    finally:
        ledger.close()
        # The same physical calculator is re-instrumented only after removing
        # this benchmark's wrappers. Restore methods saved by the caller.
    return {
        "arm": label,
        "wall_s_including_ledger": wall,
        "spent": ledger.spent,
        "status": result.status,
        "settings": asdict(settings),
        "points": [asdict(p) for p in result.points],
        "rejected_positions": result.rejected_positions,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dtype", choices=("float32", "float64"), required=True)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--scans", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    paths = [
        Path(__file__),
        ROOT / "explore/src/explore/probe.py",
        ROOT / "explore/src/explore/ledger.py",
        ROOT / "benchmarks/b04_sn2_chloride/calculator.py",
    ]
    record = {
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "dtype": args.dtype,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS")},
        "versions": {
            k: importlib.metadata.version(k) for k in ("torch", "mace-torch", "ase", "numpy")
        },
        "model": MODELS["polar1"],
        "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "fixed_geometries": {},
        "scans": [],
        "complete": False,
    }
    write(args.out, record)
    try:
        start = time.perf_counter()
        print("phase=load_calculator", flush=True)
        calc = polar1("cuda", args.dtype)
        record["load_s"] = time.perf_counter() - start
        print(f"phase=calculator_ready load_s={record['load_s']:.3f}", flush=True)
        record["gpu"] = torch.cuda.get_device_name()
        source = Path(sys.modules[type(calc).__module__].__file__)
        record["calculator_source"] = str(source)
        record["calculator_source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        systems = {
            "b04": ("benchmarks/b04_sn2_chloride/input.extxyz", -1),
            "b09": ("benchmarks/carbene_ambimodal/probe_b09/input.extxyz", 0),
        }
        atoms_by_system = {}
        for name, (path, charge) in systems.items():
            config = SearchConfig(charge_sensitive=True, total_charge=charge, multiplicity=1)
            atoms = _stamped_with_charge(read(ROOT / path), config)
            atoms_by_system[name] = atoms
            fixed_geometries(calc, atoms, args.repeats, args.out, record, name)
        if args.scans:
            # ABBA balances first/last position. Every started scan is retained.
            for k, enabled in enumerate((True, False, False, True)):
                original = {
                    key: getattr(calc, key)
                    for key in ("calculate", "get_property", "get_hessian")
                    if hasattr(calc, key)
                }
                try:
                    record["scans"].append(
                        scan(calc, atoms_by_system["b04"], enabled, k, args.out)
                    )
                    write(args.out, record)
                finally:
                    for key, method in original.items():
                        setattr(calc, key, method)
                    if hasattr(calc, "_explore_ledger"):
                        del calc._explore_ledger
        record["complete"] = True
    except BaseException as exc:
        record["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        write(args.out, record)
        faulthandler.cancel_dump_traceback_later()
    print(json.dumps({"complete": True, "output": str(args.out)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
