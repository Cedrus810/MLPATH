"""Read-only summary of midpoint_coarse_candidate_validate.py output (v5 vs v6): first-run
status transitions, per-item outcomes over all repeats, candidate agreement where both arms
found one (input-order Kabsch RMSD and symmetric RMSD over the pose-graph automorphisms, as
candidate_identity.py), how each v6 candidate was reached (coarse or fine stage), fine-stage
forwards after a coarse pair that already held a candidate, where each remaining failure
stopped, and cost / wall time per group and arm (first runs). Run on the host that wrote the
files; names are de-duplicated (183 NFS listing repeats entries).

    python midpoint_coarse_candidate_validate_summary.py OUT_JSON RESULT_DIR RUN_DIR
"""

import collections
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402

from explore.ts_ef import _rmsd  # noqa: E402
from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

NAME = re.compile(r"^(failed|control)_(\d+)_(v[56])(?:_r(\d))?\.json$")
FINE = ("from_left", "from_right")


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def wasted_fine(rec):
    """Fine-stage forwards after a coarse pair that held a candidate and passed the gate."""
    mids, total, i = rec["midpoints"], 0, 0
    while i < len(mids) - 1:
        a, b = mids[i], mids[i + 1]
        if a["origin"] == "from_left_coarse" and b["origin"] == "from_right_coarse":
            fine = [m for m in mids[i + 2 : i + 4] if m["origin"] in FINE]
            if not a["coarse_sign"] and not a.get("coarse_candidate") and fine:
                if min(a["fmax"], b["fmax"]) <= 0.01:
                    total += sum(m["steps"] + 1 for m in fine)
            i += 2 + len(fine)
        else:
            i += 1
    return total


def main():
    out, d, run_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    c = SearchConfig.from_dict(manifest["config"])
    pose = json.loads((run_dir / "poses.jsonl").read_text().splitlines()[0])
    template = Atoms(numbers=pose["numbers"], positions=pose["positions"])
    group, _ = automorphisms(template, c.bond_scale, c.active_atoms, c.automorphism_limit)

    def srmsd(x, y):
        a, b = template.copy(), template.copy()
        a.positions, b.positions = np.asarray(x), np.asarray(y)
        return float(symmetric_rmsd(a, b, group, c.active_atoms))

    items = json.loads((d / "items.json").read_text())
    runs = collections.defaultdict(list)
    for name in sorted({Path(f).name for f in glob.glob(str(d / "*.json"))}):
        m = NAME.match(name)
        if m:
            runs[(int(m[2]), m[3])].append((int(m[4] or 0), json.loads((d / name).read_text())))
    for v in runs.values():
        v.sort(key=lambda x: x[0])
    rows = []
    transitions, cost, reached = (
        collections.Counter(),
        collections.Counter(),
        collections.Counter(),
    )
    for i, (grp, task, run) in enumerate(items):
        a = [r for _, r in runs[(i, "v5")]]
        b = [r for _, r in runs[(i, "v6")]]
        if not a or not b:
            continue
        transitions[f"{grp}: {a[0]['status']} -> {b[0]['status']}"] += 1
        for arm, r in (("v5", a[0]), ("v6", b[0])):
            cost[f"{grp}_{arm}_forward"] += r["spent"]
            cost[f"{grp}_{arm}_wall_s"] += r["wall_s"]
        cost[f"{grp}_v5_fine_after_coarse_candidate"] += wasted_fine(a[0])
        for r in b:
            if r["candidate"]:
                origin = r["candidate"]["metrics"].get("origin", "")
                reached[f"{grp}: {'coarse' if origin.endswith('_coarse') else 'fine'}"] += 1
        agree = None
        if a[0]["candidate"] and b[0]["candidate"]:
            pa, pb = a[0]["candidate"]["positions"], b[0]["candidate"]["positions"]
            agree = {
                "dE_eV": b[0]["candidate"]["energy_eV"] - a[0]["candidate"]["energy_eV"],
                "rmsd_A": _rmsd(np.array(pa), np.array(pb)),
                "symmetric_rmsd_A": srmsd(pa, pb),
            }
        stop = None
        if b[0]["candidate"] is None and b[0]["midpoints"]:
            last = b[0]["midpoints"][-1]
            stop = {k: last[k] for k in ("origin", "optimizer", "q", "reaction_force", "steps")}
            stop["fmax_projected"] = last["fmax_projected"]
        rows.append(
            {
                "item": i,
                "group": grp,
                "intent": short(task),
                "run": run,
                "v5_status": [r["status"] for r in a],
                "v6_status": [r["status"] for r in b],
                "v6_candidate_origin": [
                    (r["candidate"] or {}).get("metrics", {}).get("origin") for r in b
                ],
                "first_run_agreement": agree,
                "v6_first_run_stop": stop,
            }
        )
    summary = {
        "runs": sum(len(v) for v in runs.values()),
        "items": dict(collections.Counter(g for g, _, _ in items)),
        "recheck": json.loads((d / "recheck.json").read_text())
        if (d / "recheck.json").exists()
        else None,
        "first_run_transitions": dict(transitions),
        "first_run_cost": dict(cost),
        "v6_candidates_reached_at": dict(reached),
        "rows": rows,
    }
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
