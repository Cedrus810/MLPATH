"""Which saddles are the rescued midpoint candidates? (TODO item 9b; user 2026-11-04:
"已有两个救回候选也应先比较几何身份，再做 HC——多个意图可能找到同一个 TS，新增候选和新增事件要分别计数")

Read-only, no model. References: the saddle of every HC `connected` event in the given
D2-HF runs (journal). Candidates: every HS candidate in the given acceptance result
directories (midpoint_coarse_validate / midpoint_fire_validate JSONs). Identity uses the
run's own duplicate rule for HF events: symmetric RMSD over the automorphisms of the pose
graph (as hf_event_symmetry.py, |Aut| 1296 on b05) within basin_rmsd_A, and energy within
basin_energy_eV. References are first grouped into classes with that rule (greedy, in
journal order); candidates are matched to the nearest class, and candidates matching no
class are grouped among themselves with the same rule.

    python candidate_identity.py OUT_JSON --runs RUN_DIR ... --candidates RESULT_DIR ...
"""

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from ase import Atoms  # noqa: E402

from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

CANDIDATE_FILE = re.compile(r"^(failed|control)_(\d+)_(v[3-7])(?:_r(\d))?\.json$")


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def references(run_dirs):
    out = []
    for run_dir in map(Path, run_dirs):
        with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
            for line in fh:
                if '"HC|' not in line or '"settle"' not in line:
                    continue
                ev = json.loads(line)
                if ev.get("event") != "settle" or ev["status"] != "connected":
                    continue
                p = ev["payload"]
                out.append(
                    {
                        "source": f"{run_dir.parent.name}/{run_dir.name}",
                        "intent": short(ev["task"]),
                        "energy_eV": p["saddle_energy_eV"],
                        "broken": p["fingerprint"]["broken"],
                        "formed": p["fingerprint"]["formed"],
                        "positions": np.asarray(p["saddle_positions"], float),
                    }
                )
    return out


def candidates(result_dirs):
    out = []
    for d in map(Path, result_dirs):
        for name in sorted({Path(f).name for f in glob.glob(str(d / "*.json"))}):
            m = CANDIDATE_FILE.match(name)
            if not m:
                continue
            r = json.loads((d / name).read_text(encoding="utf-8"))
            if not r.get("candidate"):
                continue
            out.append(
                {
                    "source": f"{d.parent.name}/{d.name}/{name}",
                    "group": m[1],
                    "arm": m[3],
                    "intent": short(r["task"]),
                    "energy_eV": r["candidate"]["energy_eV"],
                    "positions": np.asarray(r["candidate"]["positions"], float),
                }
            )
    return out


def cluster(items, same):
    reps = []
    for x in items:
        for k, r in enumerate(reps):
            if same(x, r):
                x["cluster"] = k
                break
        else:
            x["cluster"] = len(reps)
            reps.append(x)
    return reps


def main():
    args = sys.argv[1:]
    out = Path(args[0])
    runs = args[args.index("--runs") + 1 : args.index("--candidates")]
    result_dirs = args[args.index("--candidates") + 1 :]
    manifest = json.loads((Path(runs[0]) / "manifest.json").read_text(encoding="utf-8"))
    c = SearchConfig.from_dict(manifest["config"])
    pose = json.loads((Path(runs[0]) / "poses.jsonl").read_text().splitlines()[0])
    template = Atoms(numbers=pose["numbers"], positions=pose["positions"])
    group, _ = automorphisms(template, c.bond_scale, c.active_atoms, c.automorphism_limit)

    def rmsd(x, y):
        a, b = template.copy(), template.copy()
        a.positions, b.positions = x["positions"], y["positions"]
        return symmetric_rmsd(a, b, group, c.active_atoms)

    def same(x, y):
        return (
            abs(x["energy_eV"] - y["energy_eV"]) <= c.basin_energy_eV
            and rmsd(x, y) <= c.basin_rmsd_A
        )

    refs = references(runs)
    classes = cluster(refs, same)
    cands = candidates(result_dirs)
    for x in cands:
        d = [rmsd(x, r) for r in classes]
        k = int(np.argmin(d))
        x["nearest_class"] = k
        x["nearest_rmsd_A"] = float(d[k])
        x["nearest_dE_eV"] = float(x["energy_eV"] - classes[k]["energy_eV"])
        x["matches_class"] = bool(
            d[k] <= c.basin_rmsd_A and abs(x["nearest_dE_eV"]) <= c.basin_energy_eV
        )
    unmatched = [x for x in cands if not x["matches_class"]]
    new_groups = cluster(unmatched, same)

    def strip(x):
        return {k: v for k, v in x.items() if k != "positions"}

    summary = {
        "automorphisms": len(group),
        "basin_rmsd_A": c.basin_rmsd_A,
        "basin_energy_eV": c.basin_energy_eV,
        "reference_events": len(refs),
        "reference_classes": [
            strip(r)
            | {
                "members": sum(1 for x in refs if x["cluster"] == k),
                "sources": sorted({x["source"] for x in refs if x["cluster"] == k}),
            }
            for k, r in enumerate(classes)
        ],
        "class_to_class_rmsd_A": [[round(rmsd(a, b), 4) for b in classes] for a in classes],
        "candidates": len(cands),
        "candidates_by_class": {
            str(k): sorted(
                {x["intent"] for x in cands if x["matches_class"] and x["nearest_class"] == k}
            )
            for k in range(len(classes))
        },
        "unmatched_groups": [
            strip(r) | {"members": [x["source"] for x in unmatched if x["cluster"] == k]}
            for k, r in enumerate(new_groups)
        ],
        "rows": [strip(x) for x in cands],
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
