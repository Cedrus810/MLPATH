"""HS stops before the first bisection midpoint (TODO item 9c'' / item 10; user 2026-11-07
「可以」 to diagnosing them). Read-only, from finished D2-HF journals with the current branch
continuation (HF-EF-branches-v4); no model.

Per run: every settled HS by status, with forwards. Then, per stop kind:
  walk_relax_unconverged  the failing trial relaxation of the terminated side (history
      segment at that side's origin and failed q): steps, transverse force first / min /
      last, energy change total and in the first / last 50 steps, aligned RMSD and max atom
      move from the trial start, last-50-step per-step max atom move; plus the side's walk:
      points accepted, reaction force at the last accepted point, smallest |RF| on the
      side, and how the other side ended. Converged walk trials give the step distribution
      for comparison.
  hysteresis  per side: termination (reached_end / jump / relax_unconverged), jump size
      against its limit, energy change and reaction force across the jump; sign pattern of
      the reaction force along each side.
  no_resolved_maximum_bracket  walk bracket pair reaction forces against the refine end
      re-relaxations.
  transverse_not_converged before any walk point (end relaxations).
No classification threshold is applied; raw values are written.

    python hs_prewalk_census.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import collections
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from explore.certify import _kabsch_rmsd  # noqa: E402
from hs_transverse_diag import segments  # noqa: E402

TOL = 0.003  # EFSettings.reaction_force_tol and transverse_fmax_eV_A of these runs


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def q(v, p):
    return float(np.quantile(v, p)) if len(v) else None


def trial_stats(seg):
    fp = np.array([f["fmax_projected"] for f in seg])
    e = np.array([f["energy_eV"] for f in seg])
    pos = np.array([f["positions"] for f in seg])
    moves = np.linalg.norm(np.diff(pos[-50:], axis=0), axis=2).max(axis=1)
    return {
        "steps": int(seg[-1]["relax_step"]),
        "fproj_first": float(fp[0]),
        "fproj_min": float(fp.min()),
        "fproj_last": float(fp[-1]),
        "dE_eV": float(e[-1] - e[0]),
        "dE_first50_eV": float(e[min(50, len(e) - 1)] - e[0]),
        "dE_last50_eV": float(e[-1] - e[max(0, len(e) - 51)]),
        "rmsd_A": float(_kabsch_rmsd(pos[-1], pos[0])),
        "max_atom_move_A": float(np.linalg.norm(pos[-1] - pos[0], axis=1).max()),
        "last50_step_move_p50_A": q(moves, 0.5) if moves.size else None,
        "last50_fproj_p50": float(np.median(fp[-50:])),
        "rf_first": float(seg[0]["reaction_force"]),
        "rf_last": float(seg[-1]["reaction_force"]),
    }


def side_view(sd):
    rf = [p["reaction_force"] for p in sd["points"]]
    signs = "".join("+" if r > TOL else "-" if r < -TOL else "0" for r in rf)
    out = {
        "origin": sd["origin"],
        "status": sd["status"],
        "termination": sd.get("termination"),
        "points": len(sd["points"]),
        "rf_last_point": rf[-1],
        "min_abs_rf": float(min(abs(r) for r in rf)),
        "rf_signs": signs,
    }
    if sd.get("jump"):
        j = sd["jump"]
        out["jump"] = {
            k: j[k] for k in ("step_q", "rmsd_A", "limit_A", "dE_eV", "reaction_force_after")
        }
    if sd.get("relax_unconverged"):
        r = sd["relax_unconverged"]
        out["relax_unconverged"] = {k: v for k, v in r.items() if k != "positions"}
    return out


def scan(run_dir):
    status = collections.Counter()
    spent = collections.Counter()
    walk_steps = []
    rows = collections.defaultdict(list)
    with open(Path(run_dir) / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HS|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle" or not ev["task"].startswith("HS|"):
                continue
            p, st = ev["payload"], ev["status"]
            status[st] += 1
            spent[st] += ev["spent"]
            segs = segments(p.get("history") or [])
            for s in segs:
                if s[0]["origin"] in ("low", "high") and s[-1]["fmax_projected"] <= TOL:
                    walk_steps.append(s[-1]["relax_step"])
            base = {"task": short(ev["task"]), "spent": ev["spent"]}
            sides = [side_view(sd) for sd in p.get("sides") or []]
            if st == "unresolved:walk_relax_unconverged":
                for sd, view in zip(p["sides"], sides):
                    if sd.get("termination") != "relax_unconverged":
                        continue
                    qf = sd["relax_unconverged"]["q_target"]
                    seg = [
                        s
                        for s in segs
                        if s[0]["origin"] == sd["origin"] and s[0]["q_target"] == qf
                    ]
                    rows["walk_relax_unconverged"].append(
                        base
                        | {
                            "side": view,
                            "other": [v for v in sides if v["origin"] != sd["origin"]],
                            "trial": trial_stats(seg[-1]) if seg else None,
                        }
                    )
            elif st == "unresolved:hysteresis":
                rows["hysteresis"].append(base | {"sides": sides})
            elif st == "unresolved:no_resolved_maximum_bracket":
                r = p.get("refine") or {}
                ends = [
                    s[-1]
                    for s in segments(r.get("history") or [])
                    if s[0]["origin"] in ("left", "right")
                ]
                walk = [sd for sd in p["sides"] if sd["status"] == "maximum_bracketed"]
                rows["no_resolved_maximum_bracket"].append(
                    base
                    | {
                        "walk_pair_rf": [pt["reaction_force"] for pt in walk[0]["points"][-2:]]
                        if walk
                        else None,
                        "refine_end_rf": [f["reaction_force"] for f in ends],
                        "refine_end_fproj": [f["fmax_projected"] for f in ends],
                    }
                )
            elif st == "unresolved:transverse_not_converged" and not p.get("sides"):
                rows["end_transverse_not_converged"].append(
                    base | {"trial": trial_stats(segs[-1]) if segs else None}
                )
    return {
        "run": str(run_dir),
        "status": dict(status),
        "spent": dict(spent),
        "walk_trial_converged_steps": {
            "n": len(walk_steps),
            "p50": q(walk_steps, 0.5),
            "p90": q(walk_steps, 0.9),
            "p99": q(walk_steps, 0.99),
        },
        "rows": rows,
    }


def main():
    out, runs = Path(sys.argv[1]), sys.argv[2:]
    results = [scan(r) for r in runs]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=1) + "\n", encoding="utf-8")
    for r in results:
        print("==", r["run"])
        print("   status", r["status"])
        print("   spent", r["spent"])
        print("   converged walk trials", r["walk_trial_converged_steps"])
        for kind, v in r["rows"].items():
            print("  ", kind, len(v))


if __name__ == "__main__":
    main()
