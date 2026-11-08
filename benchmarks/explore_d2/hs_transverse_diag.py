"""Why do HS searches stop with `unresolved:transverse_not_converged`? Offline, from the
HS histories in finished D2-HF journals (positions and fmax per relaxation step; forces
are not stored). No model.

For every HS that stopped this way, the failing relaxation is the last segment (one
`_relax` call: same origin and q, relax_step increasing) of the refine history when the
branch walk had handed a maximum bracket to `refine_ts_ef`, else of the branch history.
Each segment is judged at the threshold `_relax` used for it: bisection midpoints
(refine origins "from_left" / "from_right") at the run manifest's
hf.midpoint_transverse_fmax_eV_A (1e-3 on b05); the refine end re-relaxations ("left" /
"right") and every branch relaxation at EFSettings.transverse_fmax_eV_A (0.003).
Reported per failing segment: steps, first / min / last projected fmax, energy change
(also in fp32 energy quanta at this energy, np.spacing), Kabsch RMSD and unaligned max
atom move first -> last frame, how far the minimum fmax still fell in the last 50 steps,
and the last 50 steps themselves: projected / full fmax quantiles, frames within 1.5x and
2x the threshold, per-step max atom move, energy range in quanta, reaction-force range.
Relaxations whose last frame meets their own threshold count as converged and give the
steps they needed (cap 200), per kind. No classification threshold is applied; the raw
values are written for reading.

    python hs_transverse_diag.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from explore.certify import _kabsch_rmsd  # noqa: E402
from explore.ts_ef import EFSettings  # noqa: E402

CAP = 200
TAIL = 50
MIDPOINT_ORIGINS = ("from_left", "from_right")


def segments(history):
    out, prev = [], None
    for f in history:
        if (
            prev is None
            or f["relax_step"] <= prev["relax_step"]
            or f["q_target"] != prev["q_target"]
            or f["origin"] != prev["origin"]
        ):
            out.append([])
        out[-1].append(f)
        prev = f
    return out


def quantiles(v):
    return {k: float(np.quantile(v, q)) for k, q in (("p10", 0.1), ("p50", 0.5), ("p90", 0.9))}


def describe(seg, target):
    fp = np.array([f["fmax_projected"] for f in seg])
    ff = np.array([f["fmax_eV_A"] for f in seg])
    e = np.array([f["energy_eV"] for f in seg])
    rf = np.array([f["reaction_force"] for f in seg])
    pos = np.array([f["positions"] for f in seg])
    quantum = float(np.spacing(np.float32(abs(e[0]))))
    n = min(TAIL, len(fp))
    tail = slice(len(fp) - n, None)
    before_tail = fp[: len(fp) - n]
    moves = np.linalg.norm(np.diff(pos[tail], axis=0), axis=2).max(axis=1)
    return {
        "origin": seg[0]["origin"],
        "threshold": target,
        "steps": int(seg[-1]["relax_step"]),
        "fproj_first": float(fp[0]),
        "fproj_min": float(fp.min()),
        "fproj_last": float(fp[-1]),
        "min_over_threshold": float(fp.min() / target),
        "min_fell_in_last_50": None
        if before_tail.size == 0
        else float(before_tail.min() - fp[tail].min()),
        "dE_eV": float(e[-1] - e[0]),
        "dE_quanta": float((e[-1] - e[0]) / quantum),
        "rmsd_first_last_A": _kabsch_rmsd(pos[-1], pos[0]),
        "max_atom_move_A": float(np.linalg.norm(pos[-1] - pos[0], axis=1).max()),
        "tail": {
            "frames": n,
            "fproj": quantiles(fp[tail]),
            "fmax_full": quantiles(ff[tail]),
            "within_1.5x": int((fp[tail] <= 1.5 * target).sum()),
            "within_2x": int((fp[tail] <= 2 * target).sum()),
            "step_max_atom_move_A": quantiles(moves) if moves.size else None,
            "energy_range_quanta": float(np.ptp(e[tail]) / quantum),
            "reaction_force_range": float(np.ptp(rf[tail])),
        },
    }


def thresholds(run_dir):
    manifest = json.loads((Path(run_dir) / "manifest.json").read_text(encoding="utf-8"))
    base = EFSettings().transverse_fmax_eV_A
    mid = manifest["hf"]["midpoint_transverse_fmax_eV_A"]
    return base, base if mid is None else mid


def scan(run_dir):
    base, mid = thresholds(run_dir)
    failed = []
    steps = {"midpoint": [], "refine_end": [], "branch": []}
    attempts = dict.fromkeys(steps, 0)
    with open(Path(run_dir) / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HS|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle" or not ev["task"].startswith("HS|"):
                continue
            segs = segments(ev["payload"].get("history") or [])
            rsegs = segments((ev["payload"].get("refine") or {}).get("history") or [])
            tagged = [("branch", s, base) for s in segs] + [
                ("midpoint", s, mid)
                if s[0]["origin"] in MIDPOINT_ORIGINS
                else ("refine_end", s, base)
                for s in rsegs
            ]
            for kind, s, target in tagged:
                attempts[kind] += 1
                if s[-1]["fmax_projected"] <= target:
                    steps[kind].append(s[-1]["relax_step"])
            if ev["status"] == "unresolved:transverse_not_converged" and tagged:
                # Refine segments follow the branch ones: the last is the failing one.
                kind, s, target = tagged[-1]
                row = describe(s, target)
                row.update(
                    task=ev["task"],
                    spent=ev["spent"],
                    stage="refine" if rsegs else "branch",
                    kind=kind,
                    walked=len(segs) - 1,
                )
                failed.append(row)
    converged = {}
    for kind, v in steps.items():
        v = np.array(v)
        converged[kind] = {
            "threshold": mid if kind == "midpoint" else base,
            "attempted": attempts[kind],
            "converged": int(v.size),
            "steps_quantiles": None
            if v.size == 0
            else {
                q: float(np.quantile(v, x))
                for q, x in (("p50", 0.5), ("p90", 0.9), ("p99", 0.99), ("max", 1.0))
            },
            "steps_ge_150": int((v >= 150).sum()),
        }
    return {"run": str(run_dir), "failed": failed, "converged": converged}


def main():
    out, runs = Path(sys.argv[1]), sys.argv[2:]
    results = [scan(r) for r in runs]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=1), encoding="utf-8")
    for r in results:
        f = r["failed"]
        print("==", r["run"], "failed", len(f))
        for kind, c in r["converged"].items():
            print(
                f"   converged {kind} (<= {c['threshold']}): {c['converged']} / "
                f"{c['attempted']}, steps {c['steps_quantiles']}, >=150: {c['steps_ge_150']}"
            )
        for kind in ("midpoint", "refine_end", "branch"):
            g = [x for x in f if x["kind"] == kind]
            if not g:
                continue
            print(
                f"   failed {kind}: n={len(g)}, endpoint (walked 0): "
                f"{sum(x['walked'] == 0 for x in g)}"
            )
            cols = {
                "steps": [x["steps"] for x in g],
                "fproj_min": [x["fproj_min"] for x in g],
                "fproj_last": [x["fproj_last"] for x in g],
                "tail_fproj_p50": [x["tail"]["fproj"]["p50"] for x in g],
                "tail_within_2x": [x["tail"]["within_2x"] for x in g],
                "tail_step_p50_A": [x["tail"]["step_max_atom_move_A"]["p50"] for x in g],
                "dE_eV": [x["dE_eV"] for x in g],
                "rmsd_first_last_A": [x["rmsd_first_last_A"] for x in g],
            }
            for k, v in cols.items():
                v = np.array(v, dtype=float)
                print(
                    f"      {k}: min {v.min():.3g} p50 {np.quantile(v, 0.5):.3g} max {v.max():.3g}"
                )


if __name__ == "__main__":
    main()
