"""Throughput of batched AFIR paths vs wave size (explore.afir_batch.run_paths), real model.
Measures how many concurrent paths one GPU takes before the time per step stops being flat
(user 2026-11-14: "我们可以大量并发 因为 分子量太小了"). Only the path stage, at most
--steps steps per path, the same random start structures for every wave size.

    python throughput.py OUT_JSON --system NAME --shell R --seed S [--waves 8 32 128 256]
        [--steps 20]

Per wave size: wall per step, forwards per second, and the batch evaluator's split into
graph preparation (CPU), model call (GPU) and return; the rest is the per-row Python update
(AFIR bias, bond graph, FIRE). Forwards go to a scratch ledger next to OUT_JSON.
"""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_mc_afir import prepare  # noqa: E402  (sets sys.path)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir import WeightedDistance  # noqa: E402
from explore.afir_batch import run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import sample_poses  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--system", required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--waves", type=int, nargs="+", default=[8, 32, 128, 256])
    ap.add_argument("--steps", type=int, default=20)
    args = ap.parse_args()
    out = Path(args.out)
    work = out.parent / (out.stem + "_work")
    work.mkdir(parents=True, exist_ok=True)
    system = SYSTEMS[args.system]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    frags, confs, _ = prepare(work, system, config, atoms, calc)
    n = max(args.waves)
    poses = [p for p in sample_poses(atoms, frags, confs[0], confs[1], [args.shell], 2 * n,
                                     args.seed) if p.status == "placed"][:n]  # fmt: skip
    settings = PathSettings(noise_eV=0.004, max_steps=args.steps)
    batch = MACEBatchEvaluator(calc)
    # warm-up (CUDA kernels, allocator)
    run_paths(atoms, [poses[0].positions], [WeightedDistance.build(atoms.numbers,
              frags[0].atoms, frags[1].atoms)], batch, config, settings,
              Ledger(work / "warm.ledger.jsonl", arm="M", replay=False), "warm")  # fmt: skip
    rows = []
    for w in args.waves:
        starts = [p.positions for p in poses[:w]]
        coords = [WeightedDistance.build(atoms.numbers, frags[0].atoms, frags[1].atoms)
                  for _ in starts]  # fmt: skip
        ledger = Ledger(work / f"w{w}.ledger.jsonl", arm="M", replay=False)
        batch.calls.clear()
        t0 = time.perf_counter()
        res = run_paths(atoms, starts, coords, batch, config, settings, ledger, f"w{w}")
        wall = time.perf_counter() - t0
        calls = list(batch.calls)
        steps = len(calls)
        row = {
            "wave": w,
            "paths": len(starts),
            "forwards": ledger.spent,
            "steps": steps,
            "wall_s": wall,
            "wall_per_step_s": wall / steps,
            "forwards_per_s": ledger.spent / wall,
            "prepare_s": sum(c["prepare_s"] for c in calls),
            "model_s": sum(c["model_s"] for c in calls),
            "return_s": sum(c["return_s"] for c in calls),
            "statuses": {s: sum(r.status == s for r in res) for s in {r.status for r in res}},
        }
        row["python_s"] = wall - row["prepare_s"] - row["model_s"] - row["return_s"]
        rows.append(row)
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in row.items()}), flush=True)  # fmt: skip
    out.write_text(json.dumps({"system": args.system, "atoms": len(atoms), "rows": rows},
                              indent=1))  # fmt: skip


if __name__ == "__main__":
    main()
