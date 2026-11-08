"""b04 float32 scan noise and float32/float64 coexistence (PLAN_EXPLORE_M5 Task 7).

  python fp32_noise.py --out DIR                 # (a) scan, (b) noise, (c) coexistence
  python fp32_noise.py --out DIR --mode fp64-only

Smoke on one GPU. Scan and noise forwards go through ledger.jsonl; the coexistence
evaluations (3 geometries x 4 readings) are plain calls, listed in coexist.npz. Not a
gate: numbers for the user to freeze noise_eV (M5-D5) and to judge one-process mixed
precision."""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many  # noqa: E402
from explore.probe import Intent, ScanSettings  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

INPUT = ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"
TS = ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"
INTENT = Intent("C0-Cl5", ((0, 5),), (1.0,))
NOISE_Q = (3.2, 2.7, 2.4, 2.3)
REPEATS = 5
SIGMA_A = 0.02
SEED = 20261023
FMAX, STEPS, RISE = 0.01, 300, 3.0


def config():
    return SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)


def source():
    atoms = _stamped_with_charge(read(INPUT), config())
    atoms.info["spin"] = 1
    return atoms


def settings(deltas):
    return ScanSettings(tuple(deltas), FMAX, STEPS, RISE, reuse_curvature=True)


def evaluate(calc, geometries):
    out = []
    for atoms in geometries:
        a = atoms.copy()
        a.calc = calc
        out.append((float(a.get_potential_energy()), a.get_forces().copy()))
    return out


def geometries(peak_positions=None):
    src = source()
    ts = read(TS)
    ts.info.update(src.info)
    geos = [src, ts]
    if peak_positions is not None:
        peak = src.copy()
        peak.positions = peak_positions
        geos.append(peak)
    return geos


def run_main(out):
    rec = {"host": socket.gethostname(), "seed": SEED, "sigma_A": SIGMA_A}
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    rec["load32_s"] = time.perf_counter() - t0
    batch = MACEBatchEvaluator(calc32)
    cfg = config()
    src = source()
    q0 = INTENT.value(src.positions)
    deltas = [0.0] + [round(-0.1 * k, 10) for k in range(1, int(round((q0 - 1.9) / 0.1)) + 1)]
    ledger_path = out / "ledger.jsonl"
    scan_job = ScanJob(
        src,
        INTENT,
        settings(deltas),
        cfg,
        Ledger(ledger_path, "fp32-scan", replay=False),
        "scan",
    )
    t0 = time.perf_counter()
    (scan,), timing = scan_many([scan_job], batch, 1)
    rec["scan_wall_s"] = time.perf_counter() - t0
    if scan.error is not None:
        rec["scan_error"] = scan.error
        (out / "noise.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
        return
    points = scan.result.points
    rec["scan"] = {
        "status": scan.result.status,
        "forwards": scan_job.ledger.spent,
        "q": [p.q_actual for p in points],
        "energy_eV": [p.energy_eV for p in points],
        "converged": [p.converged for p in points],
    }
    rng = np.random.default_rng(SEED)
    jobs, keys = [], []
    for q in NOISE_Q:
        base = min(points, key=lambda p: abs(p.q_target - q))
        for r in range(REPEATS):
            a = src.copy()
            a.positions = base.positions + rng.normal(0.0, SIGMA_A, base.positions.shape)
            back = settings((0.0, base.q_target - INTENT.value(a.positions)))
            jobs.append(
                ScanJob(
                    a,
                    INTENT,
                    back,
                    cfg,
                    Ledger(ledger_path, "fp32-noise", replay=False),
                    f"noise|{q}|{r}",
                )
            )
            keys.append((q, base))
    t0 = time.perf_counter()
    outcomes, _ = scan_many(jobs, batch, 4)
    rec["noise_wall_s"] = time.perf_counter() - t0
    noise = {}
    for (q, base), job, o in zip(keys, jobs, outcomes):
        row = noise.setdefault(
            str(q),
            {
                "q_target": base.q_target,
                "base_eV": base.energy_eV,
                "energies_eV": [],
                "failed": [],
                "forwards": 0,
            },
        )
        row["forwards"] += job.ledger.spent
        last = None if o.result is None or not o.result.points else o.result.points[-1]
        if last is not None and last.converged:
            row["energies_eV"].append(last.energy_eV)
        else:
            row["failed"].append(o.error or (o.result.status if o.result else None))
    for row in noise.values():
        e = np.asarray(row["energies_eV"])
        row["ptp_eV"] = float(np.ptp(e)) if len(e) >= 2 else None
        row["std_eV"] = float(e.std(ddof=1)) if len(e) >= 2 else None
        row["max_dev_from_base_eV"] = (
            float(np.max(np.abs(e - row["base_eV"]))) if len(e) else None
        )
    rec["noise"] = noise
    peak = max(points, key=lambda p: p.energy_eV)
    geos = geometries(peak.positions)
    before = evaluate(calc32, geos)
    t0 = time.perf_counter()
    calc64 = polar1("cuda", "float64")
    rec["load64_s"] = time.perf_counter() - t0
    after32 = evaluate(calc32, geos)
    same64 = evaluate(calc64, geos)
    np.savez(
        out / "coexist.npz",
        peak=peak.positions,
        e32_before=[v[0] for v in before],
        e32_after=[v[0] for v in after32],
        e64=[v[0] for v in same64],
        f32_before=np.stack([v[1] for v in before]),
        f32_after=np.stack([v[1] for v in after32]),
        f64=np.stack([v[1] for v in same64]),
    )
    (out / "noise.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")


def run_fp64_only(out):
    peak = np.load(out / "coexist.npz")["peak"]
    calc64 = polar1("cuda", "float64")
    values = evaluate(calc64, geometries(peak))
    np.savez(
        out / "fp64_only.npz", e64=[v[0] for v in values], f64=np.stack([v[1] for v in values])
    )
    a, b = np.load(out / "coexist.npz"), np.load(out / "fp64_only.npz")
    compare = {
        "fp64_same_process_minus_alone_eV": (a["e64"] - b["e64"]).tolist(),
        "fp64_force_maxdiff_eV_A": float(np.abs(a["f64"] - b["f64"]).max()),
        "fp32_after_minus_before_eV": (a["e32_after"] - a["e32_before"]).tolist(),
        "fp32_force_maxdiff_after_before_eV_A": float(
            np.abs(a["f32_after"] - a["f32_before"]).max()
        ),
        "fp32_minus_fp64_eV": (a["e32_after"] - a["e64"]).tolist(),
    }
    (out / "coexist.json").write_text(json.dumps(compare, indent=2), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=("main", "fp64-only"), default="main")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (run_main if args.mode == "main" else run_fp64_only)(out)


if __name__ == "__main__":
    main()
