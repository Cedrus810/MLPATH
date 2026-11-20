"""Post-search checks of AFIR runs and comparison with D2-HF-v1.5 / DB at the same budget
(same checks and references as benchmarks/explore_d2/v15_validate.py): saddle order by
`confirm_minimum` (validation config; an optional, separately accounted check -- design
§0: the MLP Hessian is not reliable, its verdict does not override the search) and
reference match (symmetric RMSD <= basin_rmsd_A, |dE| <= basin_energy_eV under
`match_group`, mirror image also tried). Coverage R(b) counts order-1 reference matches,
on the grid of B* = the D2-HF-v1.5 ledger of the same system and seed. The cost of a TS is
the cumulative spend at the end of the wave that produced it (an upper bound: the runner
records spend per wave). Also: distinct TS candidates (any order), distinct connected
reaction events.

    python afir_validate.py OUT_JSON RUN_DIR:SYSTEM:SEEDTAG [...]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks" / "explore_d2"))

from hf_event_validate import match_group, validation_config  # noqa: E402  (sets sys.path)

from ase import Atoms  # noqa: E402

import hf6_e2e_summary as e2e  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

V15 = ROOT / "benchmarks/explore_d2/results/2026-11-13/v15_validate.json"
SEED = 20261114


def main():
    out = Path(sys.argv[1])
    specs = [a.split(":") for a in sys.argv[2:]]
    v15 = json.loads(V15.read_text())
    calc = polar1("cuda", "float32")
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="V", replay=False)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    report = {}
    for run_s, name, tag in specs:
        run = Path(run_s)
        sysspec = e2e.systems(tag)[name]
        d2run = sysspec["d2"]
        config = validation_config(
            SearchConfig.from_dict(json.loads((d2run / "manifest.json").read_text())["config"])
        )
        numbers = json.loads((d2run / "poses.jsonl").read_text().splitlines()[0])["numbers"]

        def atoms_at(pos):
            a = _stamped_with_charge(Atoms(numbers=numbers, positions=np.asarray(pos, float)),
                                     config)  # fmt: skip
            a.info["spin"] = config.multiplicity
            a.calc = factory()
            return a

        refs = []
        if "refs" in sysspec:
            from ase.io import read

            for p in sysspec["refs"]:
                nm = p.parent.name if p.name == "geometry.extxyz" else p.stem
                refs.append({"name": nm, "atoms": atoms_at(read(p).positions)})
        else:
            for k, ev in enumerate(e2e.hc_connected(sysspec["ref_events"])):
                refs.append({"name": f"b05_event{k}", "atoms": atoms_at(ev["positions"])})
        for r in refs:
            r["energy"] = float(r["atoms"].get_potential_energy())
            r["group"] = match_group(r["atoms"], d2run, config)
        ts_rows, events = [], {}
        layer_events = {"single": set(), "double": set()}
        # single-pair waves, then the two-pair layer (dwave_*, 2026-11-15); spent_total is
        # cumulative over both, so costs stay comparable
        waves = [("single", f) for f in sorted(set(run.glob("wave_*.json")))] + [
            ("double", f) for f in sorted(set(run.glob("dwave_*.json")))
        ]
        for layer, wf in waves:
            wave = json.loads(wf.read_text())
            cost = wave["spent_total"]
            for path in wave["paths"]:
                for k in path.get("new_events", []):
                    events.setdefault(k, cost)
                    layer_events[layer].add(k)
                for t in path.get("ts", []):
                    if t.get("status") != "converged_candidate":
                        continue
                    a = atoms_at(t["positions"])
                    e = float(a.get_potential_energy())
                    with ledger.segment(
                        task=f"{name}|{tag}|{path['k']}", stage="hessian"
                    ) as seg:
                        _, chk = confirm_minimum(a, factory, config, SEED)
                        seg.outcome = str(chk.get("saddle_order", chk.get("reason")))
                    plain = mirror = None
                    near = []
                    for r in refs:
                        for flip in (False, True):
                            b = Atoms(numbers=numbers,
                                      positions=a.positions * (np.array([1, 1, -1]) if flip else 1))  # fmt: skip
                            d = float(
                                symmetric_rmsd(b, r["atoms"], r["group"], config.active_atoms)
                            )
                            if not flip:
                                near.append((d, e - r["energy"], r["name"]))
                            if (
                                d <= config.basin_rmsd_A
                                and abs(e - r["energy"]) <= config.basin_energy_eV
                            ):
                                if flip:
                                    mirror = mirror or r["name"]
                                else:
                                    plain = plain or r["name"]
                    best = min(near)
                    ts_rows.append({
                        "path": path["k"], "label": path.get("label"), "cost": cost,
                        "layer": layer,
                        "energy_eV": e, "order": chk.get("saddle_order"),
                        "order_reason": chk.get("reason"), "match": plain,
                        "mirror_match": None if plain else mirror,
                        "nearest_ref": {"name": best[2], "rmsd_A": best[0], "dE_eV": best[1]},
                        "link": (t.get("link") or {}).get("status"),
                        "fingerprint": {k: ((t.get("link") or {}).get("fingerprint") or {}).get(k)
                                        for k in ("broken", "formed")},
                        "positions": a.positions.tolist(),
                    })  # fmt: skip
        # distinct TS candidates (any order)
        groups = []
        for row in ts_rows:
            a = Atoms(numbers=numbers, positions=np.asarray(row["positions"]))
            g = match_group(a, d2run, config)
            for grp in groups:
                b = Atoms(numbers=numbers, positions=np.asarray(grp[0]["positions"]))
                if abs(row["energy_eV"] - grp[0]["energy_eV"]) <= config.basin_energy_eV and \
                        symmetric_rmsd(a, b, g, config.active_atoms) <= config.basin_rmsd_A:  # fmt: skip
                    grp.append(row)
                    break
            else:
                groups.append([row])
        b_star = v15[name][tag]["budget"]

        def coverage(key):
            first = {}
            for row in sorted(ts_rows, key=lambda x: x["cost"]):
                m = key(row)
                if row["order"] == 1 and m and m not in first:
                    first[m] = row["cost"]
            grid = {str(g): sum(v <= g * b_star for v in first.values()) / len(refs)
                    for g in e2e.GRID + (2.0, 5.0)}  # fmt: skip
            return first, grid

        first, R = coverage(lambda x: x["match"])
        first_m, R_m = coverage(lambda x: x["match"] or x["mirror_match"])
        spent = max((json.loads(w.read_text())["spent_total"] for _, w in waves),
                    default=0)  # fmt: skip
        report[f"{name}_{tag}"] = {
            "run": str(run), "B_star_d2_v15": b_star, "spent": spent,
            "events_single_layer": len(layer_events["single"]),
            "events_double_layer": len(layer_events["double"]),
            "events_only_double_layer": len(layer_events["double"] - layer_events["single"]),
            "references": [r["name"] for r in refs],
            "first_hit": first, "R": R, "first_hit_mirror_allowed": first_m, "R_mirror_allowed": R_m,
            "d2_v15_first_hit": v15[name][tag]["first_hit"], "d2_v15_R": v15[name][tag]["R"],
            "ts_candidates": len(ts_rows), "distinct_ts": len(groups),
            "distinct_ts_order1": sum(1 for g in groups if g[0]["order"] == 1),
            "events": len(events), "events_within_B_star": sum(c <= b_star for c in events.values()),
            "ts_rows": [{k: v for k, v in r.items() if k != "positions"} for r in ts_rows],
            "distinct_ts_energies": [g[0]["energy_eV"] for g in groups],
        }  # fmt: skip
        brief = {k: v for k, v in report[f"{name}_{tag}"].items() if k not in ("ts_rows",)}
        print(json.dumps(brief, default=str), flush=True)
    ledger.close()
    report["validation_forwards"] = ledger.spent
    report["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(report, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
