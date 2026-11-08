"""Smoke of the 2026-10-28 HF bracket fix on the real model: HS (continue_branches_ef,
as `_hs_task` calls it) on the forward jump interval of branch_bracketed items that
failed with `bracket must have positive width` in a finished D2-HF run. Read-only on the
run directory; each item has its own ledger capped at hs_cap.

    python hf_bracket_smoke.py OUT_DIR RUN_DIR N
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
    out_dir, run_dir, n = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
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
    settled = {}
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") == "settle" and ev["task"][:3] in ("L1|", "L2|", "HS|"):
            if ev["task"].startswith("HS|"):
                ev["payload"] = {"error": ev["payload"].get("error")}
            settled[ev["task"]] = ev
    run = SimpleNamespace(settled=settled, settings=s)
    failed = [
        row
        for row in _lines(run_dir / "population.jsonl")
        if settled.get(_task("HS", row), {}).get("status") == "error:ValueError"
    ][:n]
    calc = polar1("cuda", "float32")
    summary = []
    for i, row in enumerate(failed):
        classified = _classify(run, row)
        points = _points(settled[_task("L1", row)]["payload"])
        bracket = _hf_bracket(classified, points)
        ledger = Ledger(out_dir / f"item{i}_ledger.jsonl", arm="HF", budget=hf.hs_cap)
        t0 = time.monotonic()
        result = continue_branches_ef(
            poses[row["pose_id"]].copy(),
            _intent(row),
            bracket,
            view_factory(calc, ledger),
            config,
            ledger,
            f"HS|{row['intent_id']}",
            branch_settings=BranchSettings(geometry_noise_A=hf.geometry_noise_A),
            settings=_ef_settings(hf),
            ceiling_eV=points[0].energy_eV + s.certify_max_rise_eV,
        )
        wall = time.monotonic() - t0
        ledger.close()
        payload = _hs_payload(result)
        (out_dir / f"item{i}.json").write_text(json.dumps(payload) + "\n", encoding="utf-8")
        summary.append(
            {
                "intent": row["intent_id"],
                "class": classified["class"],
                "bracket_q": [float(p.q_target) for p in bracket],
                "status": result.status,
                "forwards": result.spent,
                "wall_s": wall,
            }
        )
        print(json.dumps(summary[-1]), flush=True)
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
