"""Read-only summary of midpoint_sign_resolution.py output: per-point coarse reaction-force
jitter sigma_c, coarse-minus-fine bias, and the candidate bounds 0.003 + 3 sigma_c with
their sign checks. File names are de-duplicated: on 183 the NFS listing of this directory
repeated every entry (2026-11-02).

    python midpoint_sign_resolution_summary.py OUT_JSON RESULT_DIR
"""

import collections
import glob
import json
import sys
from pathlib import Path

import numpy as np

TOL = 0.003


def quant(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], dtype=float)
    if not v.size:
        return None
    return {f"p{p}": float(np.percentile(v, p)) for p in (0, 5, 50, 90, 95, 100)} | {
        "n": int(v.size),
        "mean": float(v.mean()),
        "sd": float(v.std()),
    }


def stratum(group, rf):
    if group == "F":
        return "F"
    a = abs(rf)
    return "S<0.01" if a < 0.01 else "S0.01-0.05" if a < 0.05 else "S>0.05"


def main():
    out, d = Path(sys.argv[1]), Path(sys.argv[2])
    files = sorted(set(glob.glob(str(d / "[SF]_[0-9][0-9]_r[0-9].json"))))
    recs = {Path(f).name: json.loads(Path(f).read_text()) for f in files}
    points = collections.defaultdict(list)
    for name, r in recs.items():
        points[name[:4]].append(r)
    rows = []
    for key, reps in sorted(points.items()):
        coarse = [r["rf_coarse"] for r in reps if r["rf_coarse"] is not None]
        paired = [r for r in reps if r["rf_coarse"] is not None and r["rf_fine"] is not None]
        rows.append(
            {
                "point": key,
                "task": reps[0]["task"],
                "q": reps[0]["q"],
                "recorded_rf": reps[0]["recorded_rf"],
                "stratum": stratum(key[0], reps[0]["recorded_rf"]),
                "repeats": len(reps),
                "no_coarse_frame": len(reps) - len(coarse),
                "coarse_mean": float(np.mean(coarse)) if coarse else None,
                "sigma_c": float(np.std(coarse, ddof=1)) if len(coarse) > 1 else None,
                "fine_converged": len(paired),
                "bias_mean": float(np.mean([r["rf_coarse"] - r["rf_fine"] for r in paired]))
                if paired
                else None,
                "spent": sum(r["spent"] for r in reps),
            }
        )
    sig = [r["sigma_c"] for r in rows if r["sigma_c"] is not None]
    summary = {
        "records": len(recs),
        "points": len(rows),
        "spent": sum(r["spent"] for r in recs.values()),
        "sigma_c_by_stratum": {
            s: quant([r["sigma_c"] for r in rows if r["stratum"] == s])
            for s in ("S<0.01", "S0.01-0.05", "S>0.05", "F")
        },
        "sigma_c_all": quant(sig),
        "sigma_c_pooled": float(np.sqrt(np.mean(np.square(sig)))),
        "bias_repeat_level": quant(
            [
                r["rf_coarse"] - r["rf_fine"]
                for r in recs.values()
                if r["rf_coarse"] is not None and r["rf_fine"] is not None
            ]
        ),
        "bounds": {},
        "rows": rows,
    }
    for name, s in (
        ("median", float(np.median(sig))),
        ("pooled", summary["sigma_c_pooled"]),
        ("p90", float(np.percentile(sig, 90))),
    ):
        bound = TOL + 3 * s
        beyond = [
            r
            for r in recs.values()
            if r["rf_coarse"] is not None
            and r["rf_fine"] is not None
            and abs(r["rf_coarse"]) > bound
        ]
        f_points = [
            k
            for k, reps in points.items()
            if k[0] == "F"
            and all(r["rf_coarse"] is not None and abs(r["rf_coarse"]) > bound for r in reps)
            and len({r["rf_coarse"] > 0 for r in reps}) == 1
        ]
        summary["bounds"][name] = {
            "sigma_c": s,
            "bound": bound,
            "fine_converged_beyond_bound": len(beyond),
            "beyond_bound_opposite_fine_sign": sum(
                r["rf_coarse"] * r["rf_fine"] < 0 for r in beyond
            ),
            "F_points_all_repeats_beyond_same_sign": len(f_points),
        }
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
