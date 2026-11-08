"""Smoke of the walk retry (docs/PLAN_HF_WALK_RETRY_2026-10-28.md) on the real model:
rerun HS, as `_hs_task` calls it, for the first N items of a finished D2-HF run whose HS
ended `transverse_not_converged` during the branch walk (hs_transverse_diag.json: stage
"branch", walked > 0), and set the new outcome next to the recorded one. Read-only on
the run directory; one ledger per item, capped at hs_cap.

    python hs_walk_retry_smoke.py OUT_DIR RUN_DIR DIAG_JSON N
"""

import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402

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
from explore.ts_branch import BranchSettings, continue_branches_ef  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def main():
    out_dir, run_dir, diag, n = (
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        sys.argv[3],
        int(sys.argv[4]),
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    rows_diag = next(
        r for r in json.loads(Path(diag).read_text()) if Path(r["run"]).name == run_dir.name
    )
    picked = [x for x in rows_diag["failed"] if x["stage"] == "branch" and x["walked"] > 0][:n]
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    raw = dict(m["settings"])
    raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
    s, hf = D2Settings(**raw), HFSettings(**m["hf"])
    config = SearchConfig.from_dict(m["config"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    rows = {r["intent_id"]: r for r in _lines(run_dir / "population.jsonl")}
    settled = {}
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") == "settle" and ev["task"][:3] in ("L1|", "L2|", "HS|"):
            if ev["task"].startswith("HS|"):
                ev["payload"] = {}
            settled[ev["task"]] = ev
    run = SimpleNamespace(settled=settled, settings=s)
    calc = polar1("cuda", "float32")
    summary = []
    for i, old in enumerate(picked):
        row = rows[old["task"][3:]]
        classified = _classify(run, row)
        points = _points(settled[_task("L1", row)]["payload"])
        ledger = Ledger(out_dir / f"item{i}_ledger.jsonl", arm="HF", budget=hf.hs_cap)
        t0 = time.monotonic()
        result = continue_branches_ef(
            poses[row["pose_id"]].copy(),
            _intent(row),
            _hf_bracket(classified, points),
            view_factory(calc, ledger),
            config,
            ledger,
            old["task"],
            branch_settings=BranchSettings(geometry_noise_A=hf.geometry_noise_A),
            settings=_ef_settings(hf),
            ceiling_eV=points[0].energy_eV + s.certify_max_rise_eV,
        )
        wall = time.monotonic() - t0
        ledger.close()
        (out_dir / f"item{i}.json").write_text(
            json.dumps(_hs_payload(result)) + "\n", encoding="utf-8"
        )
        summary.append(
            {
                "task": old["task"],
                "class": classified["class"],
                "old": {"status": settled[old["task"]]["status"], "spent": old["spent"]},
                "new": {
                    "status": result.status,
                    "branch_status": result.branch_status,
                    "spent": result.spent,
                    "wall_s": wall,
                    "sides": [
                        {
                            "status": sd["status"],
                            "termination": sd.get("termination"),
                            "points": len(sd["points"]),
                            "unconverged_trials": [
                                t["step_q"] for t in sd.get("unconverged_trials", [])
                            ],
                        }
                        for sd in result.sides
                    ],
                },
            }
        )
        print(json.dumps(summary[-1]), flush=True)
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
