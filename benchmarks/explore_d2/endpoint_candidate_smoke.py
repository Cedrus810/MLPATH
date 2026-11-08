"""Targeted smoke of the endpoint candidate (HF-EF-bracket-v7;
docs/PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md §4; user 2026-11-07 「可以」 to R1).

For every HS of the given D2-HF run that stopped `unresolved:no_resolved_maximum_bracket`,
the walk bracket (last two points of the side that bracketed a maximum) is refined cold with
`refine_ts_ef` under v6 and v7 settings (midpoint_coarse_validate.settings_for), on its own
hs_cap ledger. Recorded: status, cost, the re-relaxed ends, the candidate (origin, energy,
positions). Output names follow candidate_identity.py (`failed_NN_vK.json`).

One JSON per (item, arm); an existing file is skipped, so rerunning resumes.

    python endpoint_candidate_smoke.py OUT_DIR BOUND RUN_DIR
"""

import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_coarse_validate import settings_for  # noqa: E402  (sets sys.path)
from midpoint_cold_replay import context, load  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _intent as intent_of  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.ts_ef import refine_ts_ef  # noqa: E402
from hs_transverse_diag import segments  # noqa: E402

STATUS = "unresolved:no_resolved_maximum_bracket"
ARMS = ("v6", "v7")


def brackets(run_dir):
    out = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if STATUS not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle" or ev["status"] != STATUS:
                continue
            (side,) = [s for s in ev["payload"]["sides"] if s["status"] == "maximum_bracketed"]
            pts = side["points"][-2:]
            out.append(
                {
                    "task": ev["task"],
                    "pair": [(p["q_actual"], np.asarray(p["positions"], float)) for p in pts],
                }
            )
    return out


def main():
    out_dir, bound, run_dir = Path(sys.argv[1]), float(sys.argv[2]), Path(sys.argv[3])
    items = brackets(run_dir)
    run = load(run_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    for i, item in enumerate(items):
        row, _, ceiling = context(run, item["task"])
        for arm in ARMS:
            out = out_dir / f"failed_{i:02d}_{arm}.json"
            if out.exists():
                continue
            ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=run.hf.hs_cap)
            t0 = time.monotonic()
            try:
                res = refine_ts_ef(
                    run.poses[row["pose_id"]].copy(),
                    intent_of(row),
                    [SimpleNamespace(q_target=q, positions=x) for q, x in item["pair"]],
                    view_factory(calc, ledger),
                    run.config,
                    ledger,
                    item["task"],
                    settings=settings_for(run, arm, bound),
                    ceiling_eV=ceiling,
                )
            finally:
                ledger.close()
            ends = [s[-1] for s in segments(res.history) if s[0]["origin"] in ("left", "right")]
            rec = {
                "run": run.name,
                "task": item["task"],
                "arm": arm,
                "status": res.status,
                "spent": res.spent,
                "wall_s": time.monotonic() - t0,
                "algorithm": res.provenance.get("algorithm"),
                "ends": [
                    {
                        k: f[k]
                        for k in ("origin", "reaction_force", "fmax_projected", "fmax_eV_A")
                    }
                    | {"energy_eV": f["energy_eV"], "relax_step": f["relax_step"]}
                    for f in ends
                ],
                "candidate": None
                if res.positions is None
                else {
                    "positions": res.positions.tolist(),
                    "energy_eV": res.energy_eV,
                    "metrics": {
                        k: v
                        for k, v in res.metrics.items()
                        if isinstance(v, (int, float, str, bool))
                    },
                },
            }
            out.write_text(json.dumps(rec) + "\n", encoding="utf-8")
            print(
                json.dumps(
                    {"file": out.name, "status": rec["status"], "spent": rec["spent"]}
                    | {"origin": (rec["candidate"] or {}).get("metrics", {}).get("origin")}
                    | {"ends_rf": [round(e["reaction_force"], 4) for e in rec["ends"]]}
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
