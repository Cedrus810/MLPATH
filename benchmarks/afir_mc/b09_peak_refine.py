"""Discriminating test for the b09 locator failure (DONE_LOG 2026-11-24; user "keyi"):
the String chains of the C:-Cb paths jump over the TS-2OMe region (C:-Cb 2.2 -> 1.6 A between
adjacent nodes). If the failure is only missing resolution there, refining the peak window
and climbing again should reach TS-2OMe.

For each (label, k) chain saved by e2_string_climb.py (chains.npz):
  1. window = nodes [top-2, top+2] around the highest interior node (clipped to the chain);
  2. resample the window to M nodes at equal arc length (ends = the window's two outer
     nodes, kept fixed; they are not minima -- a local window, as in adaptive NEB);
  3. short simplified String on the window (afir_string.string_climb, string_steps S), then
     the fixed-direction climb from its highest node (same climb, same 0.01 eV/A, ceiling =
     chain start energy + 3 eV);
  4. converged points: symmetric heavy-atom RMSD to the four POLAR-1 saddles and markers.
Same rule for every chain; no chemistry in the window choice. E/F only.

  python b09_peak_refine.py E2_OUT_DIR OUT_DIR [--labels L1,L2] [--ks 0,10,200] [--m 16] [--s 20]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from run_mc_afir_batch import CEILING_RISE_EV  # noqa: E402  (sets sys.path)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402

sys.path.insert(0, str(ROOT / "benchmarks/carbene_ambimodal"))
import classify_b09 as C  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.chemistry import automorphisms  # noqa: E402
from prrs.state import centred_kabsch_rmsd  # noqa: E402
from explore.afir_batch import edges_of  # noqa: E402
from explore.afir_string import resample, string_climb  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402

OPT = ROOT / "benchmarks/carbene_ambimodal/results/2026-10-08/opt"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("e2_out")
    ap.add_argument("out")
    ap.add_argument("--labels", default="n25|2-15,n28|2-15")
    ap.add_argument("--ks", default="0,10,200")
    ap.add_argument("--m", type=int, default=16)
    ap.add_argument("--s", type=int, default=20)
    # window rule: "top" = nodes [top-2, top+2] (2026-11-24 first test); "bond" = one window
    # [j-1, j+2] per interval j -> j+1 where the bond graph changes (2026-11-25, user
    # "可以 顺着往下做"): generic, the saddle of a bond-forming step sits where the graph jumps
    ap.add_argument("--window", choices=("top", "bond"), default="top")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    e2, out = Path(args.e2_out), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    system = SYSTEMS["b09"]
    config = make_config(system)
    atoms, _ = make_atoms(system, config)
    calc = polar1("cuda", "float32")
    batch = MACEBatchEvaluator(calc)
    refs = {s: read(OPT / s / "geometry.extxyz") for s in C.DFT}
    idx = C.atoms_of_interest(refs["TS-2OMe"])
    ref_e = {}
    for s_, r_ in refs.items():
        m = r_.copy()
        m.info.update({"charge": 0, "spin": 1})
        m.calc = calc
        ref_e[s_] = float(m.get_potential_energy())
    heavy = np.array([i for i, zz in enumerate(atoms.numbers) if zz > 1])

    def heavy_perms(m):
        group, _ = automorphisms(m, config.bond_scale)
        return {tuple(np.asarray(g)[heavy]) for g in group} | {tuple(heavy)}

    ref_perms = {s_: heavy_perms(r_) for s_, r_ in refs.items()}

    def rmsd(s_, b):
        """classify_b09.heavy_rmsd, with hydrogen-only permutations collapsed (fast)."""
        perms = ref_perms[s_] | heavy_perms(b)
        x = refs[s_].positions[heavy] - refs[s_].positions[heavy].mean(axis=0)
        best = np.inf
        for sign in (1.0, -1.0):
            for pm in perms:
                y = sign * b.positions[list(pm)]
                best = min(best, centred_kabsch_rmsd(x, y - y.mean(axis=0)))
        return float(best)

    z = np.load(e2 / "chains.npz", allow_pickle=True)
    rows = []
    for label in args.labels.split(","):
        key = label.replace("|", "_").replace("+", "_")
        for k in args.ks.split(","):
            x, e = z[f"x_{key}_k{k}"], z[f"e_{key}_k{k}"]
            top = int(np.argmax(e[1:-1])) + 1
            if args.window == "top":
                windows = [(max(top - 2, 0), min(top + 2, len(x) - 1), None)]
            else:
                g = [edges_of(atoms.numbers, xi, config.bond_scale) for xi in x]
                windows = [(max(j - 1, 0), min(j + 2, len(x) - 1), j)
                           for j in range(len(x) - 1) if g[j] != g[j + 1]]  # fmt: skip
            for lo, hi, jump in windows:
                nodes = resample(list(x[lo : hi + 1]), args.m)
                tag = f"{key}_k{k}_w{lo}-{hi}{args.tag}"
                led = Ledger(out / f"{tag}.ledger.jsonl", arm="E2-refine", replay=False)
                t0 = time.monotonic()
                res = string_climb(atoms, nodes, batch, led, f"refine|{label}|k{k}|{lo}-{hi}",
                                   string_steps=args.s,
                                   ceiling_eV=float(e[0]) + CEILING_RISE_EV)  # fmt: skip
                led.close()
                row = {"label": label, "k": int(k), "rule": args.window, "window": [lo, hi],
                       "jump": jump, "top": top, "string_status": res["string_status"],
                       "string_iterations": res["string_iterations"], "ci_index": res["ci_index"],
                       "ci_status": res["ci_status"], "ci_steps": res["ci_steps"],
                       "ci_rel_eV": res["ci_energy_eV"] - float(e[0]),
                       "ci_fmax_eV_A": res["ci_fmax_eV_A"], "forwards": res["forwards"],
                       "wall_s": time.monotonic() - t0}  # fmt: skip
                a = Atoms(numbers=atoms.numbers, positions=res["ci_positions"])
                row["rmsd_A"] = {s_: rmsd(s_, a) for s_ in C.DFT}
                row["dE_vs_ref_eV"] = {s_: res["ci_energy_eV"] - ref_e[s_] for s_ in C.DFT}
                row["markers"] = C.markers(a, idx)
                row["ci_positions"] = res["ci_positions"].tolist()
                rows.append(row)
                print(json.dumps({kk: v for kk, v in row.items() if kk != "ci_positions"}),
                      flush=True)  # fmt: skip
    (out / f"refine{args.tag}.json").write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
