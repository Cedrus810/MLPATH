"""Frozen-input, same-machine comparison of two branch-walk policies for a non-converging
trial relaxation (docs/PLAN_HF_WALK_RETRY_2026-10-28.md; user 2026-10-29):

  v3  halve the step and retry, end the side only below min_step_q (legacy_walk_v3.py)
  v4  end the side at the first failure (explore.ts_branch, current)

Inputs are frozen: the HS brackets of a finished D2-HF run, rebuilt with the shipped
`_classify` / `_hf_bracket` from its journal (read-only). Per run block, a seeded sample of
items whose recorded HS met at least one non-converging trial ("failed"), and a seeded
sample of items whose HS did not ("control": with no failed trial both policies run the
same code, so a difference there measures same-machine non-reproducibility). Every item is
searched with both policies on this machine, as `_hs_task` calls `continue_branches_ef`,
each with its own ledger capped at hs_cap. One JSON per (item, policy); an existing file
is skipped, so rerunning the same command resumes.

    python hs_policy_compare.py OUT_DIR N_FAILED N_CONTROL SEED RUN_DIR [RUN_DIR ...]
"""

import json
import random
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase import Atoms  # noqa: E402

import explore.ts_branch as ts_branch  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import (  # noqa: E402
    D2Settings,
    HFSettings,
    _classify,
    _ef_settings,
    _hf_bracket,
    _hs_payload,
    _intent,
    _points,
    _task,
)
from explore.ledger import Ledger, view_factory  # noqa: E402
from legacy_walk_v3 import _walk as walk_v3  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

WALKS = {"v3": walk_v3, "v4": ts_branch._walk}


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def load(run_dir):
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    raw = dict(m["settings"])
    raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    settled, failed, clean = {}, [], []
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") != "settle" or ev["task"][:3] not in ("L1|", "L2|", "HS|"):
            continue
        if ev["task"].startswith("HS|"):
            sides = ev["payload"].get("sides") or []
            hit = any(sd.get("unconverged_trials") for sd in sides)
            (failed if hit else clean).append(ev["task"])
            ev["payload"] = {"status": ev["status"], "spent": ev["spent"]}
        settled[ev["task"]] = ev
    rows = {r["intent_id"]: r for r in _lines(run_dir / "population.jsonl")}
    return SimpleNamespace(
        settings=D2Settings(**raw),
        hf=HFSettings(**m["hf"]),
        config=SearchConfig.from_dict(m["config"]),
        poses=poses,
        rows=rows,
        settled=settled,
        failed=sorted(failed),
        clean=sorted(clean),
    )


def search(run, task, policy, calc, out):
    row = run.rows[task[3:]]
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    classified = _classify(view, row)
    points = _points(run.settled[_task("L1", row)]["payload"])
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=run.hf.hs_cap)
    ts_branch._walk = WALKS[policy]
    t0 = time.monotonic()
    try:
        result = ts_branch.continue_branches_ef(
            run.poses[row["pose_id"]].copy(),
            _intent(row),
            _hf_bracket(classified, points),
            view_factory(calc, ledger),
            run.config,
            ledger,
            task,
            branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
            settings=_ef_settings(run.hf),
            ceiling_eV=points[0].energy_eV + run.settings.certify_max_rise_eV,
        )
    finally:
        ts_branch._walk = WALKS["v4"]
        ledger.close()
    payload = _hs_payload(result)
    record = {
        "task": task,
        "policy": policy,
        "class": classified["class"],
        "status": result.status,
        "branch_status": result.branch_status,
        "spent": result.spent,
        "wall_s": time.monotonic() - t0,
        "sides": [
            {
                "status": sd["status"],
                "termination": sd.get("termination"),
                "points": len(sd["points"]),
                "unconverged_trials": len(sd.get("unconverged_trials", [])),
            }
            for sd in result.sides
        ],
        "recorded": run.settled[task]["payload"],
        "candidate": None if result.refine is None else result.refine.get("positions"),
        "payload": payload,
    }
    out.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def main():
    out_dir = Path(sys.argv[1])
    n_failed, n_control, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    for run_dir in map(Path, sys.argv[5:]):
        run = load(run_dir)
        rng = random.Random(f"{seed}:{run_dir.name}")
        picks = [("failed", t) for t in rng.sample(run.failed, min(n_failed, len(run.failed)))]
        picks += [("control", t) for t in rng.sample(run.clean, min(n_control, len(run.clean)))]
        (out_dir / f"{run_dir.name}_picks.json").write_text(json.dumps(picks, indent=1) + "\n")
        for i, (group, task) in enumerate(picks):
            for policy in ("v3", "v4"):
                out = out_dir / f"{run_dir.name}_{i:02d}_{group}_{policy}.json"
                if out.exists():
                    continue
                rec = search(run, task, policy, calc, out)
                print(
                    json.dumps(
                        {k: rec[k] for k in ("task", "policy", "status", "spent", "wall_s")}
                        | {"group": group}
                    ),
                    flush=True,
                )


if __name__ == "__main__":
    main()
