"""Read-only summary of midpoint_optimizer_compare.py output: per group and arm the
convergence count, steps and cost; per near-root frame the outcome of every repeat; the
reaction force at convergence compared between repeats and between arms. File names are
de-duplicated (183 NFS listing repeated entries, 2026-11-02).

    python midpoint_optimizer_compare_summary.py OUT_JSON RESULT_DIR
"""

import collections
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

ARM_FILE = re.compile(r"^([NFS])_(\d+)(?:_m(\d+))?_(A\d)(?:_r(\d))?\.json$")


def quant(v):
    v = np.asarray(v, dtype=float)
    if not v.size:
        return None
    return {f"p{p}": float(np.percentile(v, p)) for p in (0, 50, 90, 100)} | {"n": int(v.size)}


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def main():
    out, d = Path(sys.argv[1]), Path(sys.argv[2])
    names = sorted({Path(f).name for f in glob.glob(str(d / "*.json"))})
    runs = collections.defaultdict(list)
    for name in names:
        m = ARM_FILE.match(name)
        if m:
            group, item, mid, arm, _ = m.groups()
            r = json.loads((d / name).read_text())
            runs[(group, int(item), int(mid or 0), arm)].append(r)
    hs = [json.loads((d / n).read_text()) for n in names if re.match(r"^N_\d+_hs\.json$", n)]
    arms = {}
    for group in "NFS":
        for arm in ("A0", "A1", "A2"):
            rs = [r for (g, _, _, a), v in runs.items() if g == group and a == arm for r in v]
            ok = [r for r in rs if r["converged"]]
            arms[f"{group}_{arm}"] = {
                "converged": len(ok),
                "runs": len(rs),
                "steps": quant([r["steps"] for r in ok]),
                "spent": sum(r["spent"] for r in rs),
            }
    frames, repeat_diff = [], []
    for (g, item, mid, arm), v in sorted(runs.items()):
        if g != "N":
            continue
        rf = [r["reaction_force"][-1] for r in v if r["converged"]]
        if len(rf) == 2:
            repeat_diff.append(abs(rf[0] - rf[1]))
        frames.append(
            {
                "item": item,
                "midpoint": mid,
                "arm": arm,
                "intent": short(v[0]["task"]),
                "outcomes": [
                    {
                        "converged": r["converged"],
                        "steps": r["steps"],
                        "rf": r["reaction_force"][-1],
                    }
                    for r in v
                ],
            }
        )
    arm_diff = {}
    for other in ("A1", "A2"):
        diffs = []
        for (g, item, mid, arm), v in runs.items():
            if g != "S" or arm != "A0" or not v[0]["converged"]:
                continue
            w = runs.get((g, item, mid, other))
            if w and w[0]["converged"]:
                diffs.append(abs(v[0]["reaction_force"][-1] - w[0]["reaction_force"][-1]))
        arm_diff[f"S_A0_vs_{other}"] = quant(diffs)
    summary = {
        "arm_runs": sum(len(v) for v in runs.values()),
        "N_hs_status": dict(collections.Counter(h["status"] for h in hs)),
        "arms": arms,
        "N_repeat_rf_diff": {"all_arms": quant(repeat_diff)},
        "converged_rf_diff_between_arms": arm_diff,
        "N_frames": frames,
    }
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "N_frames"}, indent=1))


if __name__ == "__main__":
    main()
