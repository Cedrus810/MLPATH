"""Are the distinct HF events of a D2-HF run symmetry-equivalent saddles?

Pairwise symmetric RMSD (graph automorphisms of the pose graph) between the saddles of all
HC `connected` events, next to their energies and bond changes, and the hf_duplicate rule's
own radius (basin_rmsd_A). Reads the journal only; no model.

    python hf_event_symmetry.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from ase import Atoms  # noqa: E402
from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402


def analyse(run_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    c = SearchConfig.from_dict(manifest["config"])
    pose = json.loads((run_dir / "poses.jsonl").read_text(encoding="utf-8").splitlines()[0])
    template = Atoms(numbers=pose["numbers"], positions=pose["positions"])
    group, _ = automorphisms(template, c.bond_scale, c.active_atoms, c.automorphism_limit)
    events = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HC|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle" and ev["status"] == "connected":
                p = ev["payload"]
                events.append(
                    {
                        "task": ev["task"],
                        "energy_eV": p["saddle_energy_eV"],
                        "broken": p["fingerprint"]["broken"],
                        "formed": p["fingerprint"]["formed"],
                        "positions": np.asarray(p["saddle_positions"], float),
                    }
                )
    n = len(events)
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = template.copy(), template.copy()
            a.positions, b.positions = events[i]["positions"], events[j]["positions"]
            d[i, j] = d[j, i] = symmetric_rmsd(a, b, group, c.active_atoms)
    return {
        "run": str(run_dir),
        "basin_rmsd_A": c.basin_rmsd_A,
        "basin_energy_eV": c.basin_energy_eV,
        "events": [{k: v for k, v in e.items() if k != "positions"} for e in events],
        "pairwise_symmetric_rmsd_A": d.round(4).tolist(),
    }


def main():
    out_path, run_dirs = Path(sys.argv[1]), sys.argv[2:]
    results = [analyse(r) for r in run_dirs]
    out_path.write_text(json.dumps(results, indent=1), encoding="utf-8")
    for r in results:
        print("==", r["run"], "basin_rmsd_A", r["basin_rmsd_A"])
        for e in r["events"]:
            print("  ", e["energy_eV"], "broken", e["broken"], "formed", e["formed"])
        print(np.array(r["pairwise_symmetric_rmsd_A"]))


if __name__ == "__main__":
    main()
