"""Saddles of a b09 AFIR run (any locator): identity against the four POLAR-1 references
and what each side connects to. Geometry only, no model call (2026-11-26).

For every ts record with status converged_candidate in RUN/{wave,dwave}_*.json:
  - symmetric heavy-atom RMSD to TS-2OMe / -cis / -C2 / -2 (automorphisms of both
    structures, hydrogen-only permutations collapsed, mirror included; = classify_b09);
  - carbene-C distances to Ca / Cb / aldehyde O;
  - the two link ends named by the carbene-C bonds (b09_saddle_check.name);
then distinct saddles: greedy clusters at symmetric RMSD <= 0.15 A (basin_rmsd_A).

  python b09_ts_summary.py RUN_DIR OUT_JSON
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [
    str(HERE),
    str(ROOT / "benchmarks/carbene_ambimodal"),
    str(ROOT / "src"),
    str(ROOT),
]
from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

import classify_b09 as C  # noqa: E402
from b09_saddle_check import name  # noqa: E402
from prrs.chemistry import automorphisms  # noqa: E402
from prrs.state import centred_kabsch_rmsd  # noqa: E402

OPT = ROOT / "benchmarks/carbene_ambimodal/results/2026-10-08/opt"
INPUT = ROOT / "benchmarks/carbene_ambimodal/results/2026-11-20/b09_afir_input/input.extxyz"


def main():
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    refs = {s: read(OPT / s / "geometry.extxyz") for s in C.DFT}
    nums = refs["TS-2OMe"].numbers
    heavy = np.array([i for i, z in enumerate(nums) if z > 1])
    idx = C.atoms_of_interest(refs["TS-2OMe"])
    inp = read(INPUT)
    from explore.afir_batch import edges_of

    ref_edges = edges_of(nums, inp.positions, 1.2)

    def perms(m):
        g, _ = automorphisms(m, 1.2)
        return {tuple(np.asarray(p)[heavy]) for p in g} | {tuple(heavy)}

    ref_p = {s: perms(r) for s, r in refs.items()}

    def rmsd(a, pa, b, pb):
        x = a.positions[heavy] - a.positions[heavy].mean(axis=0)
        best = np.inf
        for sign in (1.0, -1.0):
            for pm in pa | pb:
                y = sign * b.positions[list(pm)]
                best = min(best, centred_kabsch_rmsd(x, y - y.mean(axis=0)))
        return float(best)

    rows = []
    for fn in sorted({n for n in os.listdir(run) if n.endswith(".json") and "wave_" in n}):
        for p in json.loads((run / fn).read_text())["paths"]:
            for t in p.get("ts", []):
                if t.get("status") != "converged_candidate" or t.get("positions") is None:
                    continue
                a = Atoms(numbers=nums, positions=np.asarray(t["positions"]))
                pa = perms(a)
                link = t.get("link") or {}
                sides = {k: name(edges_of(nums, np.asarray(v["positions"]), 1.2), ref_edges)
                         for k, v in (link.get("ends") or {}).items()}  # fmt: skip
                rows.append({"wave": fn, "label": p["label"], "energy_eV": t.get("energy_eV"),
                             "link": link.get("status"), "sides": sides,
                             "rmsd_A": {s: rmsd(refs[s], ref_p[s], a, pa) for s in C.DFT},
                             "markers": C.markers(a, idx), "_a": a, "_p": pa})  # fmt: skip
    clusters = []
    for i, r in enumerate(rows):
        for c in clusters:
            if rmsd(rows[c[0]]["_a"], rows[c[0]]["_p"], r["_a"], r["_p"]) <= 0.15:
                c.append(i)
                break
        else:
            clusters.append([i])
    summary = []
    for c in clusters:
        r = rows[c[0]]
        summary.append({"n": len(c), "first": f"{r['wave']}:{r['label']}",
                        "energy_eV": r["energy_eV"], "link": r["link"], "sides": r["sides"],
                        "nearest_ref": min(r["rmsd_A"], key=r["rmsd_A"].get),
                        "rmsd_A": r["rmsd_A"], "markers": r["markers"]})  # fmt: skip
        print(json.dumps(summary[-1]), flush=True)
    for r in rows:
        r.pop("_a")
        r.pop("_p")
    out.write_text(json.dumps({"rows": rows, "clusters": summary}, indent=1) + "\n")


if __name__ == "__main__":
    main()
