"""Summary of hf6_dimer_saddle_identity.py output (read-only, no model calls).

    python hf6_dimer_saddle_identity_summary.py DIR
"""

import collections
import json
import sys
from pathlib import Path

from ase.data import chemical_symbols

ROOT = Path(__file__).resolve().parents[2]
POSES = {
    "b04": "benchmarks/explore_d2/results/2026-11-09/b04_hf5_v14_sA/poses.jsonl",
    "b05": "benchmarks/explore_d2/results/2026-11-08/b05_hf_v14_sA/poses.jsonl",
    "n1": "benchmarks/explore_d2/results/2026-11-09/n1_hf5_v12_sA_n0020/poses.jsonl",
}


def main():
    d = Path(sys.argv[1])
    out = {}
    for name, poses in POSES.items():
        nums = json.loads((ROOT / poses).read_text().splitlines()[0])["numbers"]

        def lab(e):
            return f"{chemical_symbols[nums[e[0]]]}{e[0]}-{chemical_symbols[nums[e[1]]]}{e[1]}"

        rows = []
        for f in sorted(d.glob(f"{name}_*.json")):
            r = json.loads(f.read_text())
            row = {
                "group": r["group"],
                "members": r["members"],
                "energy_eV": r["energy_eV"],
                "order": r["order"],
                "order_reason": r["order_reason"],
                "same_as_d2_event": (r.get("nearest_d2_event") or {}).get("same"),
                "d2_rmsd_A": (r.get("nearest_d2_event") or {}).get("rmsd_A"),
            }
            sides = r.get("sides") or []
            row["sides"] = [
                {
                    "converged": s["converged"],
                    "dE_eV": s["energy_eV"] - r["energy_eV"],
                    "vs_pose": [f"-{lab(e)}" for e in s["vs_pose_broken"]]
                    + [f"+{lab(e)}" for e in s["vs_pose_formed"]],
                    "components": s["components"],
                }
                for s in sides
            ]
            if r.get("ends_label"):
                fp = r["ends_fingerprint"]
                row["ends"] = r["ends_label"]
                row["ends_delta"] = [f"-{lab(e)}" for e in fp["broken"]] + [
                    f"+{lab(e)}" for e in fp["formed"]
                ]
                # one end in the pose graph = elementary step out of the reactant basin
                row["one_end_reactant"] = any(
                    s["converged"] and not s["vs_pose"] for s in row["sides"]
                )
                kind = (
                    "no_bond_change"
                    if not row["ends_delta"]
                    else "from_reactant"
                    if row["one_end_reactant"]
                    else "between_other_basins"
                )
            elif r["order"] != 1:
                kind = "order_not_1_on_recheck"
            else:
                kind = "descent_incomplete"
            row["kind"] = kind
            rows.append(row)
        out[name] = {
            "groups": len(rows),
            "kinds": dict(collections.Counter(r["kind"] for r in rows)),
            "order_reasons": dict(
                collections.Counter(str(r["order_reason"]) for r in rows if r["order"] != 1)
            ),
            "same_as_d2_event": sum(bool(r["same_as_d2_event"]) for r in rows),
            "rows": rows,
        }
    (d / "summary.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    for name, s in out.items():
        print(name, s["groups"], s["kinds"], s["order_reasons"], "same_as_d2", s["same_as_d2_event"])
        for r in s["rows"]:
            if r["kind"] in ("from_reactant", "between_other_basins"):
                print(
                    "  ",
                    r["group"],
                    len(r["members"]),
                    r["kind"],
                    r["ends_delta"],
                    [sd["vs_pose"] for sd in r["sides"]],
                    "D2" if r["same_as_d2_event"] else "",
                )


if __name__ == "__main__":
    main()
