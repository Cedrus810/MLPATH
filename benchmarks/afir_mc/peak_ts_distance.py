"""How far is a path's peak frame from the saddle its bracket converged to, and from the
other saddles of the run? (docs/PLAN_AFIR_PEAK_DEDUP_2026-11-20.md §3; user 2026-11-20
"keyi"). No model call: geometry only.

Saddles: HF-kernel candidates (`ts[].positions`, status converged_candidate) of PBS 29729
(results/2026-11-15/pbs_pair_double/<sys>_<seed>/{wave,dwave}_*.json). Saddles of one system
(both seeds, same atom order) are clustered by symmetric RMSD <= basin_rmsd_A (0.15 A, the
runner's MIN de-duplication tolerance), mirror images included.

Peak frames: the original runs did not save path frames; the replayed AFIR paths cached by
E2 (b04 / b05, results/2026-11-19/e2_pbs/<sys>_<seed>/frames/) and by the N1 String runs
(results/2026-11-19/n1_string_{smoke180,rest181}/<seed>/frames_*.npz) are used, frame index =
the original record's `top`. A replay is only used when its first frame matches the
original start energy within 0.005 eV (E2 field first_energy_diff_eV; b04 frames predate
that field and are marked unchecked).

Distances (symmetric RMSD over the union of the automorphism groups of both structures,
mirror included):
  own   = peak -> the saddle its own bracket converged to;
  other = peak -> nearest saddle of every OTHER cluster of the same system.

  python peak_ts_distance.py OUT_JSON
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]
from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402

RUN = HERE / "results/2026-11-15/pbs_pair_double"
E2 = HERE / "results/2026-11-19/e2_pbs"
N1 = [HERE / "results/2026-11-19/n1_string_smoke180", HERE / "results/2026-11-19/n1_string_rest181"]
INPUT = {
    "b04": ROOT / "benchmarks/b04_sn2_chloride/input.extxyz",
    "b05": ROOT / "benchmarks/b05_tertiary_halide/input.extxyz",
    "n1": ROOT / "benchmarks/n1_diels_alder/input.extxyz",
}
BASIN_RMSD_A = 0.15
BOND_SCALE = 1.2
FIRST_E_TOL = 0.005


def uniq(d, prefix):
    return sorted({n for n in os.listdir(d) if n.startswith(prefix) and n.endswith(".json")})


def saddles(sys_name):
    out = []
    for seed in ("sA", "sB"):
        d = RUN / f"{sys_name}_{seed}"
        for name in uniq(d, "wave_") + uniq(d, "dwave_"):
            for p in json.loads((d / name).read_text())["paths"]:
                for t in p.get("ts", []):
                    if t.get("status") == "converged_candidate" and t.get("positions") is not None:
                        out.append({"seed": seed, "label": p["label"], "top": t["top"],
                                    "link": (t.get("link") or {}).get("status"),
                                    "positions": np.asarray(t["positions"], float)})  # fmt: skip
    return out


def frames_for(sys_name, seed, label):
    """(positions[n_frames], first_energy_ok or None) of a cached replay, or None."""
    key = label.replace("|", "_").replace("+", "_")
    if sys_name == "n1":
        for base in N1:
            f = base / seed / f"frames_{key}.npz"
            if f.exists():
                return np.load(f, allow_pickle=True)["positions"], None
        return None
    f = E2 / f"{sys_name}_{seed}" / "frames" / f"{key}.npz"
    if not f.exists():
        return None
    ok = None
    e2 = json.loads((E2 / f"{sys_name}_{seed}" / "e2.json").read_text())
    for p in e2["paths"]:
        if p["label"] == label and p.get("first_energy_diff_eV") is not None:
            ok = abs(p["first_energy_diff_eV"]) <= FIRST_E_TOL
    return np.load(f, allow_pickle=True)["positions"], ok


def make_rmsd(numbers):
    cache = {}

    def group(pos):
        k = pos.tobytes()
        if k not in cache:
            g, _ = automorphisms(Atoms(numbers=numbers, positions=pos), BOND_SCALE)
            cache[k] = g
        return cache[k]

    def rmsd(p, q):
        a = Atoms(numbers=numbers, positions=p)
        perms = list(group(p)) + list(group(q))
        best = np.inf
        for sign in (1.0, -1.0):
            b = Atoms(numbers=numbers, positions=sign * q)
            best = min(best, symmetric_rmsd(a, b, perms))
        return float(best)

    return rmsd


def main():
    out_path = Path(sys.argv[1])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    report = {}
    for sys_name in ("b04", "b05", "n1"):
        numbers = read(INPUT[sys_name]).numbers
        rmsd = make_rmsd(numbers)
        ts = saddles(sys_name)
        clusters = []  # list of member indices
        for i, t in enumerate(ts):
            for c in clusters:
                if rmsd(t["positions"], ts[c[0]]["positions"]) <= BASIN_RMSD_A:
                    c.append(i)
                    break
            else:
                clusters.append([i])
        cid = {i: k for k, c in enumerate(clusters) for i in c}
        rows = []
        for i, t in enumerate(ts):
            fr = frames_for(sys_name, t["seed"], t["label"])
            if fr is None:
                continue
            pos, ok = fr
            if t["top"] >= len(pos):
                rows.append({"label": t["label"], "seed": t["seed"], "skip": "top_out_of_range"})
                continue
            peak = pos[t["top"]]
            own = rmsd(peak, t["positions"])
            other = sorted(
                (min(rmsd(peak, ts[j]["positions"]) for j in c), k)
                for k, c in enumerate(clusters)
                if k != cid[i]
            )
            rows.append({"seed": t["seed"], "label": t["label"], "top": t["top"],
                         "cluster": cid[i], "link": t["link"], "replay_start_ok": ok,
                         "own_A": own, "other_min_A": other[0][0] if other else None,
                         "other_cluster": other[0][1] if other else None})  # fmt: skip
        own = [r["own_A"] for r in rows if "own_A" in r]
        oth = [r["other_min_A"] for r in rows if r.get("other_min_A") is not None]
        report[sys_name] = {
            "n_saddles": len(ts),
            "n_clusters": len(clusters),
            "cluster_sizes": [len(c) for c in clusters],
            "n_with_frames": len(own),
            "own_max_A": max(own) if own else None,
            "own_median_A": float(np.median(own)) if own else None,
            "other_min_A": min(oth) if oth else None,
            "gap": (min(oth) - max(own)) if own and oth else None,
            "rows": rows,
        }
        print(sys_name, {k: v for k, v in report[sys_name].items() if k != "rows"}, flush=True)
    out_path.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
