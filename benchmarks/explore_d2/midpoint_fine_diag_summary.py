"""Summary of midpoint_fine_diag.py results (TODO item 9c); read-only. Run on the host that
wrote the files (NFS attribute cache); file names are deduplicated.

    python midpoint_fine_diag_summary.py RESULT_DIR [OUT_JSON]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from explore.ts_ef import _rmsd  # noqa: E402

FMAX, TRANSVERSE = 0.01, 0.003  # candidate contract of the runs (EFSettings defaults)


def q(xs, p):
    return float(np.percentile(xs, p)) if len(xs) else None


def fire_stats(m):
    d = m["fire"] or []
    fr = m["frames"]
    fp = np.array([f["fmax_projected"] for f in fr])
    tail = fp[-100:]
    slope = (
        float(np.polyfit(np.arange(len(tail)), np.log(tail), 1)[0]) if len(tail) > 10 else None
    )
    return {
        "steps": len(fr) - 1,
        "converged": m["converged"],
        "fproj_min": float(fp.min()),
        "fproj_last50_median": float(np.median(fp[-50:])),
        "frames_le_1p5e-3": int((fp <= 1.5e-3).sum()),
        "log_fproj_slope_last100": slope,
        "rf_last50_range": float(np.ptp([f["reaction_force"] for f in fr[-50:]])),
        "v_normal_fraction_median": q([x["v_normal_fraction"] for x in d], 50),
        "v_normal_fraction_p90": q([x["v_normal_fraction"] for x in d], 90),
        "step_removed_fraction_median": q([x["step_removed_fraction"] for x in d], 50),
        "step_removed_fraction_p90": q([x["step_removed_fraction"] for x in d], 90),
        "step_normal_fraction_old_median": q([x["step_normal_fraction_old"] for x in d], 50),
        "clipped_fraction": float(np.mean([x["step_clipped"] for x in d])) if d else None,
        "resets": int(sum(x["reset"] for x in d)),
        "dt_last": d[-1]["dt"] if d else None,
        "step_actual_max_atom_median_A": q([x["step_actual_max_atom_A"] for x in d], 50),
    }


def coarse_pairs(mids):
    out = []
    for a, b in zip(mids, mids[1:]):
        if a["origin"] != "from_left_coarse" or b["origin"] != "from_right_coarse":
            continue
        if not (a["frames"] and b["frames"] and a["converged"] and b["converged"]):
            continue
        fa, fb = a["frames"][-1], b["frames"][-1]
        out.append(
            {
                "q": a["q"],
                "coarse_sign": a["coarse_sign"],
                "rf": [fa["reaction_force"], fb["reaction_force"]],
                "fmax": [fa["fmax_eV_A"], fb["fmax_eV_A"]],
                "candidate_in_hand": [
                    f["fmax_eV_A"] <= FMAX and f["fmax_projected"] <= TRANSVERSE
                    for f in (fa, fb)
                ],
                "pair_rmsd_A": _rmsd(np.array(fa["positions"]), np.array(fb["positions"])),
            }
        )
    return out


def noise_view(x):
    keep = (
        "recorded_fproj",
        "fproj_of_mean",
        "component_std_max",
        "component_std_median",
        "energy_spread_eV",
        "transverse_eigs_lowest",
        "mean_force_fraction_soft",
        "mean_force_soft_only_atom_max",
        "mean_force_rest_atom_max",
        "soft_position_error_A",
        "soft_noise_position_A",
        "noise_std_on_modes_lowest",
        "roughness_on_soft_std",
    )
    out = {k: x[k] for k in keep if k in x}
    for k in ("fproj_each", "fproj_noise_only", "rf_each", "roughness_fproj"):
        if k in x:
            out[k + "_median_max"] = [float(np.median(x[k])), float(np.max(x[k]))]
    if "rf_each" in x:
        out["rf_std"] = float(np.std(x["rf_each"]))
    return out


def main():
    src = Path(sys.argv[1])
    rows = []
    for name in sorted({p.name for p in src.glob("item_*.json")}):
        if name.endswith(".ledger.jsonl"):
            continue
        rec = json.loads((src / name).read_text(encoding="utf-8"))
        mids = rec["midpoints"]
        rows.append(
            {
                "file": name,
                "intent": rec["task"].split(":")[-6] + ":" + rec["task"].split("|")[-1][:14],
                "status": rec["status"],
                "spent": rec["spent"],
                "refine_spent": rec["refine_spent"],
                "midpoints": [
                    [
                        m["origin"],
                        round(m["q"], 4),
                        m["optimizer"],
                        m["converged"],
                        len(m["frames"]) - 1,
                    ]
                    + (
                        [
                            round(m["frames"][-1]["reaction_force"], 4),
                            round(m["frames"][-1]["fmax_projected"], 5),
                            round(m["frames"][-1]["fmax_eV_A"], 4),
                        ]
                        if m["frames"]
                        else []
                    )
                    for m in mids
                ],
                "coarse_pairs": coarse_pairs(mids),
                "fire": [
                    {"origin": m["origin"], "q": m["q"]} | fire_stats(m)
                    for m in mids
                    if m["optimizer"] == "fire" and m["frames"]
                ],
                "diagnosis": [
                    {k: d[k] for k in ("midpoint", "origin", "best_frame")}
                    | {"best": noise_view(d["best"]), "last": noise_view(d["last"])}
                    for d in rec["diagnosis"]
                ],
            }
        )
    text = json.dumps(rows, indent=1)
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
