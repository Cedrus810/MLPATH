"""Independent property check of connected HF saddles (PLAN_HF4 §5: offline, the run's
journal and report are not touched). For each HC `connected` task of a finished D2-HF run
whose intent ends with one of the given suffixes, the saddle geometry is read from the
journal and checked with the same `check` used for the earlier b04 / b05 candidates
(hf_ef_bounded.py): full Hessian order via `confirm_minimum`, then, if index one, descent
on both sides along the unstable mode and the event of the two minima. Model: the run's
single float32 POLAR-1. One ledger per saddle, capped at CHECK_CAP.

    python hf_saddle_validate.py OUT_DIR RUN_DIR SUFFIX [SUFFIX ...]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_ef_bounded import CHECK_CAP, check  # noqa: E402  (also sets sys.path)

from ase import Atoms  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def main():
    out_dir, run_dir, suffixes = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    out_dir.mkdir(parents=True, exist_ok=True)
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    config = SearchConfig.from_dict(m["config"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    pose_of = {r["intent_id"]: r["pose_id"] for r in _lines(run_dir / "population.jsonl")}
    found = []
    for ev in _lines(run_dir / "journal.jsonl"):
        if (
            ev.get("event") == "settle"
            and ev["task"].startswith("HC|")
            and ev["status"] == "connected"
            and any(ev["task"].endswith(s) for s in suffixes)
        ):
            found.append(ev)
    calc = polar1("cuda", "float32")
    summary = []
    for i, ev in enumerate(found):
        name = f"saddle{i}"
        out = out_dir / f"{name}.json"
        if out.exists():
            summary.append(json.loads(out.read_text()))
            continue
        p = ev["payload"]
        cand = poses[pose_of[ev["task"][3:]]].copy()
        cand.set_constraint()
        cand.positions = np.asarray(p["saddle_positions"], float)
        ledger = Ledger(out_dir / f"{name}_ledger.jsonl", arm="HF", budget=CHECK_CAP)
        factory = view_factory(calc, ledger)
        cand.calc = factory()
        t0 = time.monotonic()
        checks = check(cand, factory, config, ledger, name)
        ledger.close()
        record = {
            "task": ev["task"],
            "saddle_energy_eV_run": p["saddle_energy_eV"],
            "run_event": {
                "label": p.get("label"),
                "broken": (p.get("fingerprint") or {}).get("broken"),
                "formed": (p.get("fingerprint") or {}).get("formed"),
            },
            "checks": checks,
            "forwards": ledger.spent,
            "wall_s": time.monotonic() - t0,
        }
        out.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
        summary.append(record)
        prop = checks["property"] or {}
        conn = checks["connection"] or {}
        print(
            json.dumps(
                {
                    "task": ev["task"][-30:],
                    "saddle_order": prop.get("saddle_order"),
                    "reason": prop.get("reason"),
                    "connection": conn.get("status"),
                    "label": conn.get("label"),
                    "forwards": ledger.spent,
                }
            ),
            flush=True,
        )
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
