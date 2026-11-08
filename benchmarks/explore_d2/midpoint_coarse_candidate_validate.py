"""Acceptance of the coarse-stage candidate on frozen brackets: v5 against v6 (v5 + coarse
candidate), paired on one machine (docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md §4;
user 2026-11-06 「可以」 to P1).

Same items and recheck rule as midpoint_fire_validate.py (11-04): the 22 failed + 15 control
intents of midpoint_coarse_validate.items (same SEED / N_CONTROL); each item's whole HS runs
once per arm; items whose arms differ in candidate presence or beyond the run's duplicate
criteria run REPEATS more times per arm (files `_rK`).

One JSON per (item, arm, repeat); an existing file is skipped, so rerunning resumes.
LIMIT > 0 keeps the first LIMIT items per group (smoke).

    python midpoint_coarse_candidate_validate.py OUT_DIR BOUND N_CONTROL SEED LIMIT REPEATS RUN_DIR ...
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_coarse_validate import items, search  # noqa: E402  (sets sys.path)
from midpoint_cold_replay import load  # noqa: E402
from midpoint_fire_validate import changed, log  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402

ARMS = ("v5", "v6")


def main():
    out_dir = Path(sys.argv[1])
    bound = float(sys.argv[2])
    n_control, seed, limit, repeats = map(int, sys.argv[3:7])
    runs = [load(Path(p)) for p in sys.argv[7:]]
    by_name = {r.name: r for r in runs}
    picks = items(runs, n_control, seed)
    if limit > 0:
        picks = [p for p in picks if p[0] == "failed"][:limit] + [
            p for p in picks if p[0] == "control"
        ][:limit]
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "items.json").write_text(json.dumps(picks, indent=1) + "\n")
    calc = polar1("cuda", "float32")

    def one(i, group, task, name, arm, k):
        out = out_dir / f"{group}_{i:02d}_{arm}{'' if k == 0 else f'_r{k}'}.json"
        if out.exists():
            return json.loads(out.read_text(encoding="utf-8"))
        rec = search(by_name[name], task, arm, bound, calc, out)
        log(rec, out)
        return rec

    first = {}
    for i, (group, task, name) in enumerate(picks):
        first[i] = {arm: one(i, group, task, name, arm, 0) for arm in ARMS}
    flagged = [
        i
        for i, (group, task, name) in enumerate(picks)
        if changed(first[i]["v5"], first[i]["v6"], by_name[name].config)
    ]
    (out_dir / "recheck.json").write_text(json.dumps(flagged) + "\n")
    print(json.dumps({"recheck": flagged}), flush=True)
    for i in flagged:
        group, task, name = picks[i]
        for k in range(1, repeats + 1):
            for arm in ARMS:
                one(i, group, task, name, arm, k)


if __name__ == "__main__":
    main()
