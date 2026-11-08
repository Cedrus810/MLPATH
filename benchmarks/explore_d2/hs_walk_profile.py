"""Branch-walk profiles of settled HS tasks (TODO items 7 / 9c'' R3; user 2026-11-07 「可」 to
an options document for hysteresis and walk_relax_unconverged). Read-only, from finished
D2-HF journals; no model.

Per HS: run, intent, HS class (payload `class`), status, forwards. Per walked side: origin,
status / termination, the reaction force at every accepted point, the forwards of every
accepted trial (relax_step + 1) in order, and, where the side ended in a jump or an
unconverged trial, the covalent bond change between the last accepted point and the
geometry after the jump / the failed trial's last frame (bonded when the distance is below
bond_scale times the sum of covalent radii, bond_scale from the run config).

    python hs_walk_profile.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase.data import covalent_radii

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def edges(numbers, positions, scale):
    x = np.asarray(positions, float)
    r = covalent_radii[np.asarray(numbers)]
    d = np.linalg.norm(x[:, None] - x[None], axis=2)
    i, j = np.where(np.triu(d < scale * (r[:, None] + r[None]), 1))
    return {(int(a), int(b)) for a, b in zip(i, j)}


def bond_change(numbers, before, after, scale):
    a, b = edges(numbers, before, scale), edges(numbers, after, scale)
    return {"broken": sorted(a - b), "formed": sorted(b - a)}


def scan(run_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    scale = manifest["config"].get("bond_scale", 1.2)
    numbers = json.loads((run_dir / "poses.jsonl").read_text().splitlines()[0])["numbers"]
    out = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HS|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle" or not ev["task"].startswith("HS|"):
                continue
            p = ev["payload"]
            sides = []
            for sd in p.get("sides") or []:
                pts = sd["points"]
                side = {
                    "origin": sd["origin"],
                    "status": sd["status"],
                    "termination": sd.get("termination"),
                    "q": [pt["q_actual"] for pt in pts],
                    "rf": [pt["reaction_force"] for pt in pts],
                    "energy": [pt["energy_eV"] for pt in pts],
                    "forwards": [pt["relax_step"] + 1 for pt in pts],
                }
                if sd.get("jump"):
                    j = sd["jump"]
                    side["jump"] = {
                        k: j[k] for k in ("rmsd_A", "dE_eV", "reaction_force_after", "step_q")
                    } | bond_change(numbers, pts[-1]["positions"], j["positions_after"], scale)
                if sd.get("relax_unconverged", {}).get("positions") is not None:
                    r = sd["relax_unconverged"]
                    side["relax_unconverged"] = {
                        k: r.get(k) for k in ("rmsd_A", "dE_eV", "relax_step")
                    } | bond_change(numbers, pts[-1]["positions"], r["positions"], scale)
                sides.append(side)
            out.append(
                {
                    "run": run_dir.name,
                    "intent": short(ev["task"]),
                    "class": p.get("class"),
                    "status": ev["status"],
                    "branch_status": p.get("branch_status"),
                    "spent": ev["spent"],
                    "refine_spent": (p.get("refine") or {}).get("spent"),
                    "sides": sides,
                }
            )
    return out


def main():
    out, runs = Path(sys.argv[1]), sys.argv[2:]
    rows = [r for run in runs for r in scan(run)]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows) + "\n", encoding="utf-8")
    print(len(rows), "HS rows")


if __name__ == "__main__":
    main()
