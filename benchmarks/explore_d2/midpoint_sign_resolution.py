"""Sign-resolution bound for two-stage midpoints (docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md
§2; user 2026-11-02: "0.003 + 3σ 先作为待验证公式 ... 还要比较粗、精松弛的反力偏差").

Per sampled midpoint of finished D2-HF runs (read-only), REPS repeats: the recorded first
frame of the midpoint relaxation, perturbed by N(0, SIGMA_A) per coordinate (seeded), is
relaxed at the same q with `ts_ef._relax` to the fine threshold (run manifest, 1e-3 on
b05), cold BFGS, the run's float32 POLAR-1. One relaxation gives both stages: the reaction
force at its first frame with transverse fmax <= COARSE (what the coarse stage returns) and
at its end when it meets the fine threshold.
  sigma_c  spread of the coarse reaction force over the repeats of one point;
  bias     coarse minus fine reaction force within a repeat (where fine converged).

Points: converged midpoints in three strata of recorded |RF| (< 0.01, 0.01-0.05, > 0.05
eV/A), N_PER_STRATUM each (SEED), and one failed midpoint per intent (last run block).

One JSON per (point, repeat); an existing file is skipped, so rerunning resumes. GROUP
selects S (converged strata) or F (failed); LIMIT > 0 keeps the first LIMIT points (smoke).

    python midpoint_sign_resolution.py OUT_DIR GROUP REPS N_PER_STRATUM SEED LIMIT RUN_DIR ...
"""

import json
import random
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import explore.ts_ef as ts_ef  # noqa: E402
from midpoint_cold_replay import COLD_BUDGET, context, load  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _ef_settings, _intent  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402
from prrs.reliability import GateRejected  # noqa: E402

COARSE = 2e-3
SIGMA_A = 0.01
STRATA = ((0.0, 0.01), (0.01, 0.05), (0.05, np.inf))


def points(runs, group, n_per, seed):
    if group == "F":
        latest = {}
        for run in runs:
            for f in run.failed:
                latest[f["task"]] = f
        return [latest[k] for k in sorted(latest)]
    rng = random.Random(seed)
    pool = [f for r in runs for f in r.converged]
    picked = []
    for lo, hi in STRATA:
        stratum = [f for f in pool if lo <= abs(f["recorded_rf"]) < hi]
        picked += rng.sample(stratum, min(n_per, len(stratum)))
    return picked


def repeat(run, point, k, calc, out):
    row, _, ceiling = context(run, point["task"])
    rng = np.random.default_rng([int(abs(point["q"]) * 1e9) % 2**32, point["midpoint"], k])
    start = np.asarray(point["positions"], dtype=float)
    start = start + rng.normal(0.0, SIGMA_A, start.shape)
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=COLD_BUDGET)
    history, stop, t0 = [], None, time.monotonic()
    try:
        with ledger.segment(task=point["task"], stage="midpoint_sign_resolution") as seg:
            try:
                ts_ef._relax(
                    run.poses[row["pose_id"]].copy(),
                    _intent(row),
                    start,
                    point["q"],
                    point["origin"],
                    view_factory(calc, ledger),
                    run.config,
                    _ef_settings(run.hf),
                    history,
                    ceiling,
                    run.hf.midpoint_transverse_fmax_eV_A,
                    state=None,
                )
            except ts_ef._Stop as exc:
                stop = str(exc)
            except BudgetExhausted:
                stop = "budget_exhausted"
            except GateRejected as exc:
                stop = f"gate_rejected:{exc.code}"
            seg.outcome = stop or "converged"
    finally:
        ledger.close()
    first = next((f for f in history if f["fmax_projected"] <= COARSE), None)
    record = {
        k2: point[k2]
        for k2 in ("run", "task", "midpoint", "origin", "q", "recorded_rf", "recorded_steps")
    } | {
        "repeat": k,
        "fine_converged": stop is None,
        "stop": stop,
        "steps": history[-1]["relax_step"] if history else None,
        "rf_coarse": None if first is None else first["reaction_force"],
        "coarse_step": None if first is None else first["relax_step"],
        "fproj_coarse": None if first is None else first["fmax_projected"],
        "rf_fine": history[-1]["reaction_force"] if stop is None else None,
        "rf_trace": [f["reaction_force"] for f in history],
        "fproj_trace": [f["fmax_projected"] for f in history],
        "spent": ledger.spent,
        "wall_s": time.monotonic() - t0,
    }
    out.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def main():
    out_dir, group = Path(sys.argv[1]), sys.argv[2]
    reps, n_per, seed, limit = map(int, sys.argv[3:7])
    runs = [load(Path(p)) for p in sys.argv[7:]]
    by_name = {r.name: r for r in runs}
    picks = points(runs, group, n_per, seed)
    picks = picks[:limit] if limit > 0 else picks
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{group}_points.json").write_text(
        json.dumps(
            [{k: p[k] for k in ("run", "task", "midpoint", "q", "recorded_rf")} for p in picks],
            indent=1,
        )
        + "\n"
    )
    calc = polar1("cuda", "float32")
    for i, p in enumerate(picks):
        for k in range(reps):
            out = out_dir / f"{group}_{i:02d}_r{k}.json"
            if out.exists():
                continue
            rec = repeat(by_name[p["run"]], p, k, calc, out)
            print(
                json.dumps(
                    {"file": out.name}
                    | {
                        x: rec[x]
                        for x in ("rf_coarse", "rf_fine", "coarse_step", "steps", "spent")
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
