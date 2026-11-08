"""Acceptance of two-stage midpoints on frozen brackets, old and new on one machine
(docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md §4; user 2026-11-02: "验收用冻结括区的新旧同机对照
... 跑完整个后续二分，统计救回数量、后续卡点、候选一致性和计算成本").

Items (read-only from finished D2-HF runs, last run block per intent):
  failed   every intent whose recorded HS stopped at a midpoint (transverse_not_converged);
  control  a SEED sample of N_CONTROL intents whose recorded HS gave converged_candidate.
Each item's whole HS (branch walk and every later bisection step) is rerun as `_hs_task`
runs it, once per arm, with its own hs_cap ledger:
  v3  the run's EF settings (single-stage midpoints);
  v4  the same plus midpoint_coarse_transverse_fmax_eV_A = COARSE and
      midpoint_sign_bound_eV_A = BOUND.
Recorded per run: status, cost, candidate (positions, energy, metrics) and every midpoint
(stage, q, reaction forces, transverse force, steps, coarse decision).

One JSON per (item, arm); an existing file is skipped, so rerunning resumes. LIMIT > 0
keeps the first LIMIT items per group (smoke).

    python midpoint_coarse_validate.py OUT_DIR BOUND N_CONTROL SEED LIMIT RUN_DIR [RUN_DIR ...]
"""

import json
import random
import sys
import time
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import explore.ts_branch as ts_branch  # noqa: E402
from hs_transverse_diag import segments  # noqa: E402
from midpoint_cold_replay import context, load  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _classify, _ef_settings, _hf_bracket  # noqa: E402
from explore.d2 import _intent as intent_of  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402

COARSE = 2e-3
MIDPOINTS = ("from_left", "from_right", "from_left_coarse", "from_right_coarse")


def items(runs, n_control, seed):
    failed, control = {}, {}
    for run in runs:
        for f in run.failed:
            failed[f["task"]] = run.name
        for task, ev in run.settled.items():
            if task.startswith("HS|") and ev["payload"]["status"] == "converged_candidate":
                control[task] = run.name
    pool = sorted(t for t in control if t not in failed)
    picked = random.Random(seed).sample(pool, min(n_control, len(pool)))
    return [("failed", t, failed[t]) for t in sorted(failed)] + [
        ("control", t, control[t]) for t in sorted(picked)
    ]


def midpoints(history):
    out = []
    for s in segments(history):
        if s[0]["origin"] not in MIDPOINTS:
            continue
        last = s[-1]
        out.append(
            {
                "origin": s[0]["origin"],
                "q": s[0]["q_target"],
                "steps": last["relax_step"],
                "reaction_force": last["reaction_force"],
                "fmax_projected": last["fmax_projected"],
                "fmax": last["fmax_eV_A"],
                "coarse_sign": last.get("coarse_sign"),
                "coarse_candidate": bool(last.get("coarse_candidate")),
                "optimizer": last.get("optimizer"),
            }
        )
    return out


def settings_for(run, arm, bound):
    """v3: the run's settings; v4: two-stage midpoints; v5: v4 with a FIRE fine stage;
    v6: v5 with the coarse-stage candidate (docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md);
    v7: v6 with the endpoint candidate (docs/PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md)."""
    settings = _ef_settings(run.hf)
    if arm in ("v4", "v5", "v6", "v7"):
        settings = replace(
            settings,
            midpoint_coarse_transverse_fmax_eV_A=COARSE,
            midpoint_sign_bound_eV_A=bound,
            midpoint_fine_optimizer="bfgs" if arm == "v4" else "fire",
            midpoint_coarse_candidate=arm in ("v6", "v7"),
            endpoint_candidate=arm == "v7",
        )
    return settings


def search(run, task, arm, bound, calc, out):
    row, points, ceiling = context(run, task)
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    classified = _classify(view, row)
    settings = settings_for(run, arm, bound)
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=run.hf.hs_cap)
    t0 = time.monotonic()
    try:
        result = ts_branch.continue_branches_ef(
            run.poses[row["pose_id"]].copy(),
            intent_of(row),
            _hf_bracket(classified, points),
            view_factory(calc, ledger),
            run.config,
            ledger,
            task,
            branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
            settings=settings,
            ceiling_eV=ceiling,
        )
    finally:
        ledger.close()
    refine = result.refine or {}
    record = {
        "run": run.name,
        "task": task,
        "arm": arm,
        "bound": bound if arm in ("v4", "v5", "v6", "v7") else None,
        "status": result.status,
        "branch_status": result.branch_status,
        "spent": result.spent,
        "wall_s": time.monotonic() - t0,
        "recorded_status": run.settled[task]["payload"]["status"],
        "refine_spent": refine.get("spent"),
        "candidate": None
        if refine.get("positions") is None
        else {
            "positions": [list(map(float, r)) for r in refine["positions"]],
            "energy_eV": refine.get("energy_eV"),
            "metrics": {
                k: v
                for k, v in (refine.get("metrics") or {}).items()
                if isinstance(v, (int, float, str))
            },
        },
        "algorithm": (refine.get("provenance") or {}).get("algorithm"),
        "midpoints": midpoints(refine.get("history") or []),
    }
    out.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def main():
    out_dir = Path(sys.argv[1])
    bound = float(sys.argv[2])
    n_control, seed, limit = map(int, sys.argv[3:6])
    runs = [load(Path(p)) for p in sys.argv[6:]]
    by_name = {r.name: r for r in runs}
    picks = items(runs, n_control, seed)
    if limit > 0:
        picks = [p for p in picks if p[0] == "failed"][:limit] + [
            p for p in picks if p[0] == "control"
        ][:limit]
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "items.json").write_text(json.dumps(picks, indent=1) + "\n")
    calc = polar1("cuda", "float32")
    for i, (group, task, name) in enumerate(picks):
        for arm in ("v3", "v4"):
            out = out_dir / f"{group}_{i:02d}_{arm}.json"
            if out.exists():
                continue
            rec = search(by_name[name], task, arm, bound, calc, out)
            print(
                json.dumps(
                    {"file": out.name, "status": rec["status"], "spent": rec["spent"]}
                    | {"mids": [(m["origin"][-6:], m["coarse_sign"]) for m in rec["midpoints"]]}
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
