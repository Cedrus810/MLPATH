"""Serial against parallel HS on frozen brackets, with GPU utilisation
(docs/PLAN_HF_D2_V14_2026-11-07.md §3; user 2026-11-07 「可以」 to the plan).

Items: the 22 failed + 15 control intents of the 11-04 acceptance (items.json, same source
runs). Every arm runs each item's whole HS (`continue_branches_ef`) with the D2-HF-v1.3 HS
settings (`_ef_settings(hf, HSBracket(2e-3, 0.025))`, i.e. HF-EF-bracket-v7), each on its
own hs_cap ledger:
  serial  one item after another through one counted view of the resident calculator, as
          D2 runs HS up to v1.3;
  pK      all items through `run_many` with max_active K and the resident
          MACEBatchEvaluator, as D2-HF-v1.4 runs a wave.
Recorded per arm: wall time; per item status, forwards, candidate energy and positions;
`run_many` timing (batch sizes, evaluator wall, outside-evaluator wall); `nvidia-smi`
utilisation and memory sampled every second while the arm runs.

One JSON per arm; an existing file is skipped, so rerunning resumes.

    python hs_parallel_bench.py OUT_DIR ITEMS_JSON RESULTS_DIR LIMIT ARMS

LIMIT > 0 keeps the first LIMIT items per group (smoke). ARMS: comma list, e.g. serial,p4,p8.
"""

import json
import subprocess
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_cold_replay import context, load  # noqa: E402  (sets sys.path)

import explore.ts_branch as ts_branch  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import HSBracket, _classify, _ef_settings, _hf_bracket  # noqa: E402
from explore.d2 import _intent as intent_of  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.parallel import MACEBatchEvaluator, TaskJob, run_many  # noqa: E402

BRACKET = HSBracket(2e-3, 0.025)


class GPUSampler:
    def __init__(self):
        self.samples = []
        self.proc = subprocess.Popen(
            [
                "nvidia-smi",
                "--query-gpu=utilization.gpu,memory.used",
                "--format=csv,noheader,nounits",
                "-l",
                "1",
            ],
            stdout=subprocess.PIPE,
            text=True,
        )
        self.thread = threading.Thread(target=self._read, daemon=True)
        self.thread.start()

    def _read(self):
        for line in self.proc.stdout:
            try:
                util, mem = (float(x) for x in line.split(","))
            except ValueError:
                continue
            self.samples.append((time.time(), util, mem))

    def stop(self):
        self.proc.terminate()
        self.proc.wait()
        u = np.array([s[1] for s in self.samples]) if self.samples else np.array([])
        return {
            "n": int(u.size),
            "util_mean": float(u.mean()) if u.size else None,
            "util_p50": float(np.median(u)) if u.size else None,
            "util_p90": float(np.percentile(u, 90)) if u.size else None,
            "mem_max_MiB": max((s[2] for s in self.samples), default=None),
        }


def prepare(run, task):
    row, points, ceiling = context(run, task)
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    bracket = _hf_bracket(_classify(view, row), points)
    return row, bracket, ceiling


def hs(run, task, row, bracket, ceiling, factory, ledger):
    result = ts_branch.continue_branches_ef(
        run.poses[row["pose_id"]].copy(),
        intent_of(row),
        bracket,
        factory,
        run.config,
        ledger,
        task,
        branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
        settings=_ef_settings(run.hf, BRACKET),
        ceiling_eV=ceiling,
    )
    refine = result.refine or {}
    return {
        "status": result.status,
        "algorithm": (refine.get("provenance") or {}).get("algorithm"),
        "candidate_energy_eV": refine.get("energy_eV"),
        "candidate_positions": None
        if refine.get("positions") is None
        else np.asarray(refine["positions"]).tolist(),
    }


def main():
    out_dir, picks = Path(sys.argv[1]), json.loads(Path(sys.argv[2]).read_text())
    results, limit, arms = Path(sys.argv[3]), int(sys.argv[4]), sys.argv[5].split(",")
    if limit > 0:
        picks = [p for p in picks if p[0] == "failed"][:limit] + [
            p for p in picks if p[0] == "control"
        ][:limit]
    names = sorted({name for _, _, name in picks})
    runs = {n: load(next(results.glob(f"*/{n}"))) for n in names}
    items = [(g, t, runs[n]) + prepare(runs[n], t) for g, t, n in picks]
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    batch = MACEBatchEvaluator(calc)
    for arm in arms:
        out = out_dir / f"arm_{arm}.json"
        if out.exists():
            continue
        ledgers = [
            Ledger(out_dir / f"{arm}_{i:02d}.ledger.jsonl", arm="HF", budget=r.hf.hs_cap)
            for i, (_, _, r, *_rest) in enumerate(items)
        ]
        sampler, t0, timing = GPUSampler(), time.monotonic(), None
        records = []
        if arm == "serial":
            for (g, t, r, row, br, ceil), led in zip(items, ledgers):
                s0 = time.monotonic()
                try:
                    rec = hs(r, t, row, br, ceil, view_factory(calc, led), led)
                except Exception as exc:  # keep the arm going; record the failure
                    rec = {"status": f"error:{type(exc).__name__}", "error": str(exc)}
                records.append(rec | {"wall_s": time.monotonic() - s0})
        else:
            k = int(arm[1:])
            jobs = [
                TaskJob(
                    lambda factory, r=r, t=t, row=row, br=br, ceil=ceil, led=led: hs(
                        r, t, row, br, ceil, factory, led
                    ),
                    led,
                )
                for (g, t, r, row, br, ceil), led in zip(items, ledgers)
            ]
            outcomes, timing = run_many(jobs, batch, k)
            for o in outcomes:
                records.append(
                    o.result if o.error is None else {"status": "error", "error": o.error}
                )
        wall = time.monotonic() - t0
        gpu = sampler.stop()
        for led in ledgers:
            led.close()
        sizes = [rd["size"] for rd in (timing or {}).get("batch_rounds", [])]
        summary = {
            "arm": arm,
            "items": len(items),
            "wall_s": wall,
            "forwards": sum(led.spent for led in ledgers),
            "s_per_forward": wall / max(1, sum(led.spent for led in ledgers)),
            "gpu": gpu,
            "batch_size_mean": float(np.mean(sizes)) if sizes else None,
            "batch_rounds": len(sizes),
            "evaluator_wall_s": None
            if timing is None
            else sum(rd["wall_s"] for rd in timing["batch_rounds"]),
            "outside_evaluator_wall_s": None
            if timing is None
            else timing["outside_evaluator_wall_s"],
            "statuses": dict(
                zip(*np.unique([rec["status"] for rec in records], return_counts=True))
            ),
        }
        summary["statuses"] = {k: int(v) for k, v in summary["statuses"].items()}
        rows = [
            {"group": g, "task": t, "run": r.name, "spent": led.spent} | rec
            for (g, t, r, *_rest), led, rec in zip(items, ledgers, records)
        ]
        out.write_text(json.dumps({"summary": summary, "rows": rows}) + "\n", encoding="utf-8")
        print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
