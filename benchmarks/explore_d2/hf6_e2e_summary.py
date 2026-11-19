"""HF-6 phase 2 summary: D2 against DB on the same L0 and budget (docs/PLAN_HF6_DIMER_E2E_2026-11-11.md
§3). After the search, both arms get the same checks as hf6_local_summary.py (saddle order by
`confirm_minimum` on float32 POLAR-1; reference match by symmetric RMSD <= basin_rmsd_A and
|dE| <= basin_energy_eV). Only order-1 saddles matching a reference count towards coverage
R(b) = |found refs with hit cost <= b| / |S|; the hit cost of an item is the arm's cumulative
forward count when it was produced (D2: the HC settle's `cumulative_spent`; DB: the attempt's
`cumulative_forwards`). Reported on the budget grid 0, 1/8, 1/4, 1/2, 3/4, 1 of B*; first-hit
cost per reference (None = right-censored at B*).

    python hf6_e2e_summary.py OUT_JSON SEED_TAG          (SEED_TAG: sA or sB)
"""

import collections
import json
import sys
import time
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_event_validate import _jsonable, match_group, validation_config  # noqa: E402  (sets sys.path)

from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "benchmarks/explore_d2/results"
N1REF = ROOT / "benchmarks/n1_diels_alder/results/2026-11-09"
GRID = (0.0, 0.125, 0.25, 0.5, 0.75, 1.0)
SEED = 20261111


def systems(tag):
    db = RES / "2026-11-11" / ("hf6_e2e_db" if tag == "sA" else "hf6_e2e_db_sB")
    d2 = {
        "sA": {
            "b04": RES / "2026-11-09/b04_hf5_v14_sA",
            "b05": RES / "2026-11-08/b05_hf_v14_sA",
            "n1": RES / "2026-11-09/n1_hf5_v12_sA_n0020",
        },
        "sB": {
            "b04": RES / "2026-11-09/b04_hf5_v14_sB",
            "b05": RES / "2026-11-08/b05_hf_v14_sB",
            "n1": RES / "2026-11-11/n1_hf5_v12_sB_n0020",
        },
    }[tag]
    return {
        "b04": {
            "d2": d2["b04"],
            "db": db / "b04",
            "refs": [ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"],
        },
        "b05": {
            "d2": d2["b05"],
            "db": db / "b05",
            "ref_events": RES / "2026-11-08/b05_hf_v14_sA",
        },
        "n1": {
            "d2": d2["n1"],
            "db": db / "n1",
            "refs": [
                N1REF / f"ts_{k}/geometry.extxyz"
                for k in ("endo_s-cis", "endo_s-trans", "exo_s-cis", "exo_s-trans")
            ],
        },
    }


def hc_connected(run_dir):
    out = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HC|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle" and ev["status"] == "connected":
                p = ev["payload"]
                out.append(
                    {
                        "positions": p["saddle_positions"],
                        "cost": ev["cumulative_spent"],
                        "task": ev["task"],
                    }
                )
    return out


def main():
    out, tag = Path(sys.argv[1]), sys.argv[2]
    calc = polar1("cuda", "float32")
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="V", replay=False)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    report = {"seed_tag": tag}
    for name, spec in systems(tag).items():
        manifest = json.loads((spec["d2"] / "manifest.json").read_text())
        config = validation_config(SearchConfig.from_dict(manifest["config"]))
        numbers = json.loads((spec["d2"] / "poses.jsonl").read_text().splitlines()[0])[
            "numbers"
        ]
        budget = json.loads((spec["d2"] / "report.json").read_text())["ledger_total"]

        def atoms_at(positions):
            a = _stamped_with_charge(Atoms(numbers=numbers, positions=positions), config)
            a.info["spin"] = config.multiplicity
            a.calc = factory()
            return a

        refs = []
        if "refs" in spec:
            for p in spec["refs"]:
                refs.append(
                    {
                        "name": p.parent.name if p.name == "geometry.extxyz" else p.stem,
                        "atoms": atoms_at(read(p).positions),
                    }
                )
        else:
            for k, ev in enumerate(hc_connected(spec["ref_events"])):
                refs.append({"name": f"b05_event{k}", "atoms": atoms_at(ev["positions"])})
        for r in refs:
            r["energy"] = float(r["atoms"].get_potential_energy())
            r["group"] = match_group(r["atoms"], spec["d2"], config)

        def check(positions, label):
            a = atoms_at(positions)
            e = float(a.get_potential_energy())
            with ledger.segment(task=f"{name}|{label}", stage="hessian", attempt=0) as seg:
                _, chk = confirm_minimum(a, factory, config, SEED)
                seg.outcome = str(chk.get("saddle_order", chk.get("reason")))
            match, near = None, []
            for r in refs:
                d = float(symmetric_rmsd(a, r["atoms"], r["group"], config.active_atoms))
                de = e - r["energy"]
                near.append((d, de, r["name"]))
                if d <= config.basin_rmsd_A and abs(de) <= config.basin_energy_eV:
                    match = r["name"]
            best = min(near)
            return a, {
                "energy_eV": e,
                "order": chk.get("saddle_order"),
                "order_reason": chk.get("reason"),
                "imaginary_icm": _jsonable(chk.get("imaginary_wavenumbers_icm")),
                "match": match,
                "nearest_ref": {"name": best[2], "rmsd_A": best[0], "dE_eV": best[1]},
            }

        arms = {}
        # D2: HC-connected events
        d2c = []
        for k, ev in enumerate(hc_connected(spec["d2"])):
            a, cell = check(ev["positions"], f"d2:{k}")
            d2c.append((ev["cost"], a, cell))
        # DB: converged attempts
        dbc, status = [], collections.Counter()
        for f in sorted(spec["db"].glob("[0-9][0-9][0-9]_[0-9][0-9][0-9].json")):
            rec = json.loads(f.read_text())
            status[rec["status"]] += 1
            if rec["status"] == "converged":
                a, cell = check(rec["positions"], f"db:{f.stem}")
                dbc.append(
                    (
                        rec["cumulative_forwards"],
                        a,
                        cell | {"start_displacement_A": rec.get("start_displacement_A")},
                    )
                )
        db_done = json.loads((spec["db"] / "done.json").read_text())

        for arm, cands in (("d2", d2c), ("db", dbc)):
            first = {}
            for cost, a, cell in sorted(cands, key=lambda t: t[0]):
                if cell["order"] == 1 and cell["match"] and cell["match"] not in first:
                    first[cell["match"]] = cost
            groups = []
            for cost, a, cell in cands:
                if cell["order"] != 1:
                    continue
                g = match_group(a, spec["d2"], config)
                for grp in groups:
                    b, eb = grp[0]
                    if (
                        abs(cell["energy_eV"] - eb) <= config.basin_energy_eV
                        and symmetric_rmsd(a, b, g, config.active_atoms) <= config.basin_rmsd_A
                    ):
                        grp.append((a, cell["energy_eV"]))
                        break
                else:
                    groups.append([(a, cell["energy_eV"])])
            arms[arm] = {
                "candidates": len(cands),
                "order_1": sum(c["order"] == 1 for _, _, c in cands),
                "order_other": dict(
                    collections.Counter(
                        str(c["order"] if c["order"] is not None else c["order_reason"])
                        for _, _, c in cands
                        if c["order"] != 1
                    )
                ),
                "distinct_order1": len(groups),
                "order1_not_reference": sum(
                    c["order"] == 1 and not c["match"] for _, _, c in cands
                ),
                "first_hit": {r["name"]: first.get(r["name"]) for r in refs},
                "R": {
                    str(g): sum(1 for v in first.values() if v <= g * budget) / len(refs)
                    for g in GRID
                },
                "candidates_detail": [{"cost": cost} | cell for cost, _, cell in cands],
            }
        arms["db"]["attempt_status"] = dict(status)
        arms["db"]["done"] = db_done
        report[name] = {"budget": budget, "references": [r["name"] for r in refs], "arms": arms}
        brief = {
            arm: {k: v for k, v in s.items() if k not in ("candidates_detail",)}
            for arm, s in arms.items()
        }
        print(name, budget, json.dumps(brief, default=str), flush=True)
    ledger.close()
    report["validation_forwards"] = ledger.spent
    report["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(report, indent=1, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
