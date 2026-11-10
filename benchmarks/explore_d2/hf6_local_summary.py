"""HF-6 phase 1 summary (docs/PLAN_HF6_DIMER_2026-11-10.md §3.3). After the search: every
candidate of both arms (hf `converged_candidate`, dimer `converged`) gets the same checks as
hf_event_validate.py -- saddle order by `prrs.runner.confirm_minimum` on float32 POLAR-1, and
match to each reference saddle (symmetric RMSD under the reference's own graph automorphisms
<= basin_rmsd_A and |dE| <= basin_energy_eV, the runs' duplicate rule). Order-1 candidates of
one arm are grouped into distinct saddles with the same rule. Per item, the two arms' outcomes
are cross-tabulated.

    python hf6_local_summary.py OUT_JSON
"""

import collections
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_event_validate import _jsonable, events  # noqa: E402  (sets sys.path)

from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "benchmarks/explore_d2/results"
N1 = ROOT / "benchmarks/n1_diels_alder/results/2026-11-09"
SYSTEMS = {
    "b04": {
        "result": RES / "2026-11-10/hf6_local/b04",
        "run": RES / "2026-11-08/b04_hf5_v12_sA",
        "refs": [ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"],
    },
    "b05": {
        "result": RES / "2026-11-10/hf6_local/b05",
        "run": RES / "2026-11-08/b05_hf_v14_sA",
        "ref_events": RES / "2026-11-08/b05_hf_v14_sA",
    },
    "n1": {
        "result": RES / "2026-11-10/hf6_local/n1",
        "run": RES / "2026-11-09/n1_hf5_v12_sA_n0020",
        "refs": [
            N1 / f"ts_{k}/geometry.extxyz"
            for k in ("endo_s-cis", "endo_s-trans", "exo_s-cis", "exo_s-trans")
        ],
    },
}
SEED = 20261111


def main():
    out = Path(sys.argv[1])
    calc = polar1("cuda", "float32")
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="V", replay=False)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    report = {}
    for name, spec in SYSTEMS.items():
        manifest = json.loads((spec["run"] / "manifest.json").read_text())
        config = SearchConfig.from_dict(manifest["config"])
        numbers = json.loads((spec["run"] / "poses.jsonl").read_text().splitlines()[0])[
            "numbers"
        ]

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
            for k, ev in enumerate(events(spec["ref_events"])):
                refs.append(
                    {
                        "name": f"event{k}:{ev['broken']}->{ev['formed']}",
                        "atoms": atoms_at(ev["positions"]),
                    }
                )
        for r in refs:
            r["energy"] = float(r["atoms"].get_potential_energy())
            r["group"], _ = automorphisms(
                r["atoms"], config.bond_scale, config.active_atoms, config.automorphism_limit
            )

        def rmsd(a, b, group):
            return float(symmetric_rmsd(a, b, group, config.active_atoms))

        items = json.loads((spec["result"] / "items.json").read_text())
        rows, arms = [], {"hf": [], "dimer": []}
        for i in range(len(items)):
            row = {"item": i}
            for arm in ("hf", "dimer"):
                rec = json.loads((spec["result"] / f"{i:03d}_{arm}.json").read_text())
                ok = rec["status"] == ("converged_candidate" if arm == "hf" else "converged")
                cell = {
                    "status": rec["status"],
                    "forwards": rec["forwards"],
                    "candidate": ok and rec.get("positions") is not None,
                }
                if cell["candidate"]:
                    a = atoms_at(rec["positions"])
                    e = float(a.get_potential_energy())
                    with ledger.segment(
                        task=f"{name}|{i}|{arm}", stage="hessian", attempt=0
                    ) as seg:
                        _, chk = confirm_minimum(a, factory, config, SEED)
                        seg.outcome = str(chk.get("saddle_order", chk.get("reason")))
                    match = None
                    near = []
                    for r in refs:
                        d, de = rmsd(a, r["atoms"], r["group"]), e - r["energy"]
                        near.append((d, de, r["name"]))
                        if d <= config.basin_rmsd_A and abs(de) <= config.basin_energy_eV:
                            match = r["name"]
                    best = min(near)
                    cell |= {
                        "energy_eV": e,
                        "order": chk.get("saddle_order"),
                        "order_reason": chk.get("reason"),
                        "imaginary_icm": _jsonable(chk.get("imaginary_wavenumbers_icm")),
                        "match": match,
                        "nearest_ref": {"name": best[2], "rmsd_A": best[0], "dE_eV": best[1]},
                        "start_displacement_A": rec.get("start_displacement_A"),
                    }
                    arms[arm].append((i, a, e, cell))
                row[arm] = cell
            rows.append(row)

        def distinct(cands):
            groups = []
            for i, a, e, cell in cands:
                if cell["order"] != 1:
                    continue
                g, _ = automorphisms(
                    a, config.bond_scale, config.active_atoms, config.automorphism_limit
                )
                for grp in groups:
                    b, eb = grp[0]
                    if (
                        abs(e - eb) <= config.basin_energy_eV
                        and rmsd(a, b, g) <= config.basin_rmsd_A
                    ):
                        grp.append((a, e))
                        break
                else:
                    groups.append([(a, e)])
            return [{"energy_eV": grp[0][1], "members": len(grp)} for grp in groups]

        summ = {}
        for arm in ("hf", "dimer"):
            cells = [r[arm] for r in rows]
            cand = [c for c in cells if c["candidate"]]
            succ = [c for c in cand if c["order"] == 1]
            summ[arm] = {
                "items": len(cells),
                "status": dict(collections.Counter(c["status"] for c in cells)),
                "candidates": len(cand),
                "order_1": len(succ),
                "order_other": dict(
                    collections.Counter(
                        str(c["order"] if c["order"] is not None else c["order_reason"])
                        for c in cand
                        if c["order"] != 1
                    )
                ),
                "matched": dict(collections.Counter(c["match"] for c in cand if c["match"])),
                "refs_found": sorted({c["match"] for c in cand if c["match"]}),
                "distinct_order1_saddles": distinct(arms[arm]),
                "forwards_total": sum(c["forwards"] for c in cells),
                "forwards_per_order1_median": float(np.median([c["forwards"] for c in succ]))
                if succ
                else None,
            }
        cross = collections.Counter(
            (
                "hf:"
                + (
                    "ref"
                    if r["hf"].get("match")
                    else "saddle"
                    if r["hf"].get("order") == 1
                    else "no"
                ),
                "dimer:"
                + (
                    "ref"
                    if r["dimer"].get("match")
                    else "saddle"
                    if r["dimer"].get("order") == 1
                    else "no"
                ),
            )
            for r in rows
        )
        report[name] = {
            "references": [{"name": r["name"], "energy_eV": r["energy"]} for r in refs],
            "arms": summ,
            "cross": {f"{a} / {b}": n for (a, b), n in sorted(cross.items())},
            "rows": rows,
        }
        print(
            name,
            json.dumps(
                {
                    "arms": {
                        k: {kk: v for kk, v in s.items() if kk != "distinct_order1_saddles"}
                        | {"distinct": len(s["distinct_order1_saddles"])}
                        for k, s in summ.items()
                    },
                    "cross": report[name]["cross"],
                },
                default=str,
            ),
            flush=True,
        )
    ledger.close()
    report["validation_forwards"] = ledger.spent
    report["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(report, indent=1, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
