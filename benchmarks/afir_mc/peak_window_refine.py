"""Bond-change peak windows on saved String chains, any system (2026-11-25; user "可以 顺着往下
做", generality test of the b09 result in b09_peak_refine.py --window bond).

For each (label, k) chain in CHAINS (chains.npz of e2_string_climb.py: x_<key>_k<k>,
e_<key>_k<k>): one window [j-1, j+2] per interval j -> j+1 where the bond graph changes;
resample to M nodes (window ends fixed), short simplified String (S iterations), fixed-
direction climb from its highest node (ceiling = chain start + 3 eV). Converged points are
written per k as OUT/k<k>/wave_0000.json in the e2_string_climb format, so afir_validate.py
identifies them (order, reference match) unchanged. E/F only in the search.

  python peak_window_refine.py CHAINS_NPZ SYSTEM OUT --labels L1,L2 [--ks 0,10] [--m 16] [--s 20]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_mc_afir_batch import CEILING_RISE_EV  # noqa: E402  (sets sys.path)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import edges_of  # noqa: E402
from explore.afir_string import resample, string_climb  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import _jsonable  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chains")
    ap.add_argument("system")
    ap.add_argument("out")
    ap.add_argument("--labels", required=True)
    ap.add_argument("--ks", default="0,10")
    ap.add_argument("--m", type=int, default=16)
    ap.add_argument("--s", type=int, default=20)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    system = SYSTEMS[args.system]
    config = make_config(system)
    atoms, _ = make_atoms(system, config)
    batch = MACEBatchEvaluator(polar1("cuda", "float32"))
    z = np.load(args.chains, allow_pickle=True)
    rows, waves, spent = [], {}, {}
    for label in args.labels.split(","):
        key = label.replace("|", "_").replace("+", "_")
        for k in args.ks.split(","):
            if f"x_{key}_k{k}" not in z.files:
                rows.append({"label": label, "k": int(k), "status": "no_chain"})
                continue
            x, e = z[f"x_{key}_k{k}"], z[f"e_{key}_k{k}"]
            g = [edges_of(atoms.numbers, xi, config.bond_scale) for xi in x]
            jumps = [j for j in range(len(x) - 1) if g[j] != g[j + 1]]
            for j in jumps:
                lo, hi = max(j - 1, 0), min(j + 2, len(x) - 1)
                nodes = resample(list(x[lo : hi + 1]), args.m)
                led = Ledger(out / f"{key}_k{k}_w{lo}-{hi}.ledger.jsonl", arm="E2-refine",
                             replay=False)  # fmt: skip
                t0 = time.monotonic()
                res = string_climb(atoms, nodes, batch, led, f"refine|{label}|k{k}|{lo}-{hi}",
                                   string_steps=args.s,
                                   ceiling_eV=float(e[0]) + CEILING_RISE_EV)  # fmt: skip
                led.close()
                spent[k] = spent.get(k, 0) + res["forwards"]
                row = {"label": label, "k": int(k), "window": [lo, hi], "jump": j,
                       "ci_status": res["ci_status"], "ci_steps": res["ci_steps"],
                       "ci_rel_eV": res["ci_energy_eV"] - float(e[0]),
                       "ci_fmax_eV_A": res["ci_fmax_eV_A"], "forwards": res["forwards"],
                       "wall_s": time.monotonic() - t0}  # fmt: skip
                rows.append(row)
                print(json.dumps(row), flush=True)
                if res["ci_status"] == "converged":
                    waves.setdefault(k, []).append({
                        "k": len(waves.get(k, [])), "label": f"{label}#w{lo}-{hi}",
                        "new_events": [],
                        "ts": [{"status": "converged_candidate",
                                "positions": res["ci_positions"],
                                "energy_eV": res["ci_energy_eV"], "link": {}}],
                    })  # fmt: skip
    for k in args.ks.split(","):
        d = out / f"k{k}"
        d.mkdir(parents=True, exist_ok=True)
        (d / "wave_0000.json").write_text(json.dumps(_jsonable(
            {"wave": 0, "paths": waves.get(k, []), "spent_total": spent.get(k, 0)}),
            default=str) + "\n")  # fmt: skip
    (out / "refine.json").write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
