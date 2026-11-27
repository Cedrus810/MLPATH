"""Does the String chain of a b09 path pass near TS-2OMe? (user 2026-11-24 "做诊断").
Geometry only, no model call. For each checkpoint chain (chains.npz of e2_string_climb.py)
and for the raw replayed AFIR path (frames/*.npz): per node the stored energy (when saved),
the carbene-C distances to Ca / Cb / aldehyde O, and the symmetric heavy-atom RMSD to the four
POLAR-1 saddles (classify_b09.heavy_rmsd, automorphisms + mirror).

  python b09_chain_diag.py E2_OUT_DIR LABEL OUT_JSON
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "benchmarks/carbene_ambimodal"), str(ROOT / "src"), str(ROOT)]
from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

import classify_b09 as C  # noqa: E402

OPT = ROOT / "benchmarks/carbene_ambimodal/results/2026-10-08/opt"


def main():
    d, label, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    refs = {s: read(OPT / s / "geometry.extxyz") for s in C.DFT}
    nums = refs["TS-2OMe"].numbers
    idx = C.atoms_of_interest(refs["TS-2OMe"])
    key = label.replace("|", "_").replace("+", "_")
    z = np.load(d / "chains.npz", allow_pickle=True)
    rows = {}

    from prrs.chemistry import automorphisms
    from prrs.state import centred_kabsch_rmsd

    heavy = np.array([i for i, zz in enumerate(nums) if zz > 1])

    def heavy_perms(atoms):
        """Distinct heavy-atom images of the automorphisms (hydrogen-only permutations
        collapse, which is what makes this fast); same minimum as classify_b09.heavy_rmsd."""
        group, _ = automorphisms(atoms, 1.2)
        return {tuple(np.asarray(g)[heavy]) for g in group} | {tuple(heavy)}

    ref_perms = {s: heavy_perms(refs[s]) for s in C.DFT}

    def rmsd(ref, perms, b):
        x = refs[ref].positions[heavy] - refs[ref].positions[heavy].mean(axis=0)
        best = np.inf
        for sign in (1.0, -1.0):
            for pm in perms:
                y = sign * b.positions[list(pm)]
                best = min(best, centred_kabsch_rmsd(x, y - y.mean(axis=0)))
        return float(best)

    def describe(pos, energies=None):
        res = []
        for i, x in enumerate(pos):
            a = Atoms(numbers=nums, positions=x)
            own = heavy_perms(a)
            r = {s: rmsd(s, ref_perms[s] | own, a) for s in C.DFT}
            res.append({"node": i, "E_eV": None if energies is None else float(energies[i]),
                        **{m: v for m, v in C.markers(a, idx).items()},
                        "rmsd": r})  # fmt: skip
        return res

    fr = np.load(d / "frames" / f"{key}.npz", allow_pickle=True)
    raw = fr["positions"]
    step = max(1, len(raw) // 40)
    sel = list(range(0, len(raw), step)) + [len(raw) - 1]
    rows["afir_raw"] = describe(raw[sel], fr["energies"][sel])
    for i, row in zip(sel, rows["afir_raw"]):
        row["frame"] = int(i)
    for name in sorted(z.files):
        if name.startswith("x_" + key + "_k"):
            k = name.rsplit("_k", 1)[1]
            e = z.get("e_" + key + "_k" + k)
            rows["k" + k] = describe(z[name], e)
    summary = {}
    for name, res in rows.items():
        best = min(res, key=lambda r: r["rmsd"]["TS-2OMe"])
        summary[name] = {
            "n": len(res),
            "closest_to_TS-2OMe": {
                "node": best.get("frame", best["node"]),
                "rmsd_A": best["rmsd"]["TS-2OMe"],
                "C:-Cb": best["C:-Cb"],
                "C:-Ca": best["C:-Ca"],
                "C:-O_ald": best["C:-O_ald"],
            },
        }
        print(name, summary[name], flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"summary": summary, "rows": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main()
