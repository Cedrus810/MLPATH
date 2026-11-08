"""Replay the duplicate rule, as fixed on 2026-10-28 (saddle + endpoint symmetries,
`explore.certify.pair_symmetries`), over the certified events of finished D2 runs.

In journal order, every event the run kept as distinct (HC `connected` / L4 `c2`) is
re-judged against the earlier events kept by the replay. HF runs call the fixed
`explore.d2._hf_duplicate` itself; Hessian runs apply certify step 5 as written (same
chemical pair, |dE| <= basin_energy_eV, symmetric RMSD <= basin_rmsd_A under
`pair_symmetries(saddle, endpoints)`). Reads journals only; no model.

    python replay_duplicates.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from ase import Atoms  # noqa: E402
from explore.certify import pair_symmetries  # noqa: E402
from explore.d2 import _hf_duplicate  # noqa: E402
from prrs.chemistry import symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402


def events(run_dir, hf):
    head, kept = ("HC|", "connected") if hf else ("L4|", "c2")
    out = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if f'"{head}' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle" and ev["task"].startswith(head) and ev["status"] == kept:
                out.append(ev)
    return out


def replay(run_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    c = SearchConfig.from_dict(manifest["config"])
    poses = {}
    for line in (run_dir / "poses.jsonl").read_text(encoding="utf-8").splitlines():
        r = json.loads(line)
        poses[r["pose_id"]] = Atoms(numbers=r["numbers"], positions=r["positions"])
    pose_of = {}
    for line in (run_dir / "population.jsonl").read_text(encoding="utf-8").splitlines():
        r = json.loads(line)
        pose_of[r["intent_id"]] = r["pose_id"]
    hf = manifest["rule_version"].startswith("D2-HF")
    kept, verdicts = [], []
    for ev in events(run_dir, hf):
        p = ev["payload"]
        template = poses[pose_of[ev["task"][3:]]]
        pair = tuple(p["chemical_pair"])
        if hf:
            ends = [sd["descent"]["positions"] for sd in p["sides"]]
            run = SimpleNamespace(config=c, settled={k["task"]: k for k in kept})
            dup = _hf_duplicate(run, template, p["saddle_positions"], p["saddle_energy_eV"], pair, ends)
        else:
            saddle = template.copy()
            saddle.positions = np.asarray(p["saddle_positions"], float)
            minima = []
            for side in p["sides"]:
                end = template.copy()
                end.positions = np.asarray(side["positions"], float)
                minima.append(end)
            group, dup = None, None
            for k in kept:
                q = k["payload"]
                if tuple(q["chemical_pair"]) != pair:
                    continue
                if abs(q["saddle_energy_eV"] - p["saddle_energy_eV"]) > c.basin_energy_eV:
                    continue
                if group is None:
                    group, _ = pair_symmetries(saddle, minima, c)
                other = saddle.copy()
                other.positions = np.asarray(q["saddle_positions"], float)
                if symmetric_rmsd(saddle, other, group, c.active_atoms) <= c.basin_rmsd_A:
                    dup = k["task"]
                    break
        verdicts.append({"task": ev["task"], "energy_eV": p["saddle_energy_eV"], "duplicate_of": dup})
        if dup is None:
            kept.append(ev)
    return {
        "run": str(run_dir),
        "rule_version": manifest["rule_version"],
        "events_kept_by_run": len(verdicts),
        "events_after_replay": len(kept),
        "verdicts": verdicts,
    }


def main():
    out_path, run_dirs = Path(sys.argv[1]), sys.argv[2:]
    results = [replay(r) for r in run_dirs]
    out_path.write_text(json.dumps(results, indent=1), encoding="utf-8")
    for r in results:
        print(r["run"], r["rule_version"], "distinct events:", r["events_kept_by_run"], "->", r["events_after_replay"])


if __name__ == "__main__":
    main()
