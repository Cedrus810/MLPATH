"""D2-HF-v1.5 three-system rerun: post-search checks and comparison with the previous D2 runs
(docs/PLAN_HF_D2_V15_2026-11-13.md §3). Same checks and references as hf6_e2e_summary.py:
saddle order by `confirm_minimum` (validation config, analytic cross-check 0.1) and reference
match (symmetric RMSD <= basin_rmsd_A, |dE| <= basin_energy_eV under `match_group`). Also a
mirror-allowed match (z -> -z; DONE_LOG 2026-11-12: the N1 references keep one enantiomer),
reported separately. Coverage R(b) on the grid of B* = this run's ledger total, counting only
order-1 reference matches (strict, as HF-6). Previous runs: the HF-6 summaries' D2 arm.

    python v15_validate.py OUT_JSON
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_event_validate import _jsonable, match_group, validation_config  # noqa: E402  (sets sys.path)

from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

import hf6_e2e_summary as e2e  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

RES = e2e.RES
NEW = RES / "2026-11-13"
PREV = {
    "sA": RES / "2026-11-11/hf6_e2e_summary_sA_v3.json",
    "sB": RES / "2026-11-11/hf6_e2e_summary_sB_v2.json",
}
SEED = 20261113


def main():
    out = Path(sys.argv[1])
    calc = polar1("cuda", "float32")
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="V", replay=False)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    report = {}
    for tag in ("sA", "sB"):
        prev = json.loads(PREV[tag].read_text())
        for name, spec in e2e.systems(tag).items():
            run = NEW / f"{name}_v15_{tag}"
            manifest = json.loads((run / "manifest.json").read_text())
            config = validation_config(SearchConfig.from_dict(manifest["config"]))
            numbers = json.loads((run / "poses.jsonl").read_text().splitlines()[0])["numbers"]
            budget = json.loads((run / "report.json").read_text())["ledger_total"]

            def atoms_at(positions):
                a = _stamped_with_charge(
                    Atoms(numbers=numbers, positions=np.asarray(positions, float)), config
                )
                a.info["spin"] = config.multiplicity
                a.calc = factory()
                return a

            refs = []
            if "refs" in spec:
                for p in spec["refs"]:
                    name_ = p.parent.name if p.name == "geometry.extxyz" else p.stem
                    refs.append({"name": name_, "atoms": atoms_at(read(p).positions)})
            else:
                for k, ev in enumerate(e2e.hc_connected(spec["ref_events"])):
                    refs.append({"name": f"b05_event{k}", "atoms": atoms_at(ev["positions"])})
            for r in refs:
                r["energy"] = float(r["atoms"].get_potential_energy())
                r["group"] = match_group(r["atoms"], run, config)

            rows = []
            for k, ev in enumerate(e2e.hc_connected(run)):
                a = atoms_at(ev["positions"])
                e = float(a.get_potential_energy())
                with ledger.segment(
                    task=f"{name}|{tag}|{k}", stage="hessian", attempt=0
                ) as seg:
                    _, chk = confirm_minimum(a, factory, config, SEED)
                    seg.outcome = str(chk.get("saddle_order", chk.get("reason")))
                plain = mirror = None
                near = []
                for r in refs:
                    for flip in (False, True):
                        q = a.positions * (np.array([1, 1, -1]) if flip else 1)
                        b = Atoms(numbers=numbers, positions=q)
                        d = float(
                            symmetric_rmsd(b, r["atoms"], r["group"], config.active_atoms)
                        )
                        de = e - r["energy"]
                        if not flip:
                            near.append((d, de, r["name"]))
                        if d <= config.basin_rmsd_A and abs(de) <= config.basin_energy_eV:
                            if flip:
                                mirror = mirror or r["name"]
                            else:
                                plain = plain or r["name"]
                best = min(near)
                rows.append(
                    {
                        "task": ev["task"],
                        "cost": ev["cost"],
                        "energy_eV": e,
                        "order": chk.get("saddle_order"),
                        "order_reason": chk.get("reason"),
                        "imaginary_icm": _jsonable(chk.get("imaginary_wavenumbers_icm")),
                        "match": plain,
                        "mirror_match": mirror if not plain else None,
                        "nearest_ref": {"name": best[2], "rmsd_A": best[0], "dE_eV": best[1]},
                    }
                )

            def coverage(key_fn):
                first = {}
                for row in sorted(rows, key=lambda r: r["cost"]):
                    m = key_fn(row)
                    if row["order"] == 1 and m and m not in first:
                        first[m] = row["cost"]
                return first, {
                    str(g): sum(v <= g * budget for v in first.values()) / len(refs)
                    for g in e2e.GRID
                }

            first, R = coverage(lambda r: r["match"])
            first_m, R_m = coverage(lambda r: r["match"] or r["mirror_match"])
            pd2 = prev[name]["arms"]["d2"]
            report.setdefault(name, {})[tag] = {
                "budget": budget,
                "references": [r["name"] for r in refs],
                "events": rows,
                "first_hit": {r["name"]: first.get(r["name"]) for r in refs},
                "R": R,
                "first_hit_mirror_allowed": {r["name"]: first_m.get(r["name"]) for r in refs},
                "R_mirror_allowed": R_m,
                "previous": {
                    "run": str(spec["d2"].relative_to(e2e.ROOT)),
                    "budget": prev[name]["budget"],
                    "first_hit": pd2["first_hit"],
                    "R": pd2["R"],
                    "order1_not_reference": pd2["order1_not_reference"],
                    "events": len(pd2["candidates_detail"]),
                },
            }
            print(
                name,
                tag,
                budget,
                json.dumps(
                    {
                        "R": R,
                        "first": first,
                        "Rm": R_m,
                        "events": [
                            (
                                r["order"],
                                r["match"],
                                r["mirror_match"],
                                round(r["nearest_ref"]["rmsd_A"], 3),
                            )
                            for r in rows
                        ],
                        "prev_R": pd2["R"],
                    }
                ),
                flush=True,
            )
    ledger.close()
    report["validation_forwards"] = ledger.spent
    report["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(report, indent=1, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
