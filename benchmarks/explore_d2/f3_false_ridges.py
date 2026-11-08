"""Offline count of F3 false ridges in finished b05 D2 runs (TS_HESSIAN_FREE_DESIGN_RECORD §4 F3).

F3: `read_scan` takes the post-peak drop as energies[k] - min(all later energies), while its
convergence check covers only points[: k + 2]. A drop supplied by a later unconverged point
can therefore trigger `ridge_suspected`. The corrected reading used here: the drop may only
come from the contiguous run of converged points right after the peak. Everything else in
the rule (rise above the scan start, local maximum, reaction-force sign change, noise_eV) is
the current rule unchanged, so the difference is attributable to F3 alone.

Reads journals only; no model, no evaluation.

    python f3_false_ridges.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from explore.d2 import _points  # noqa: E402
from explore.probe import read_scan  # noqa: E402


def ridge_converged_drop(points, noise_eV):
    """The current ridge rule with the post-peak drop restricted to converged points."""
    energies = [p.energy_eV for p in points]
    for k in range(1, len(points) - 1):
        if not all(p.converged for p in points[: k + 2]):
            break
        end = k + 1
        while end + 1 < len(points) and points[end + 1].converged:
            end += 1
        rises = energies[k] - energies[0] > noise_eV
        falls = energies[k] - min(energies[k + 1 : end + 1]) > noise_eV
        flips = points[k - 1].reaction_force * points[k + 1].reaction_force < 0
        if (
            rises
            and falls
            and flips
            and energies[k] >= energies[k - 1]
            and energies[k] >= energies[k + 1]
        ):
            return k
    return None


def analyse(run_dir):
    run_dir = Path(run_dir).resolve()
    noise = json.loads((run_dir / "manifest.json").read_text())["settings"]["noise_eV"]
    settled = {}
    with (run_dir / "journal.jsonl").open() as fh:
        for line in fh:
            if '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev["event"] == "settle":
                settled[ev["task"]] = ev
    rows = Counter()
    false_l4 = Counter()
    true_l4 = Counter()
    examples = []
    for task, ev in settled.items():
        if not task.startswith("L1|"):
            continue
        points = _points(ev["payload"])
        if len(points) < 3:
            continue
        old = next(
            (t for t in read_scan(points, noise) if t["kind"] == "ridge_suspected"), None
        )
        new = ridge_converged_drop(points, noise)
        intent = task[3:]
        l4 = settled.get("L4|" + intent)
        l4_status = l4["status"] if l4 else "not_certified"
        if old is not None and new is None:
            rows["false_ridge_F3"] += 1
            false_l4[l4_status] += 1
            if len(examples) < 5:
                k = old["peak"]
                examples.append(
                    {
                        "intent": intent,
                        "peak": k,
                        "energies_rel_eV": [
                            round(p.energy_eV - points[0].energy_eV, 4) for p in points
                        ],
                        "converged": [p.converged for p in points],
                        "l4": l4_status,
                    }
                )
        elif old is not None:
            rows["ridge_both_rules"] += 1
            true_l4[l4_status] += 1
            if new != old["peak"]:
                rows["ridge_both_rules_peak_moved"] += 1
        elif new is not None:
            rows["ridge_only_new_rule"] += 1
        else:
            rows["no_ridge"] += 1
    return {
        "run": str(run_dir.relative_to(ROOT)),
        "noise_eV": noise,
        "counts": dict(rows),
        "l4_of_false_ridges": dict(false_l4),
        "l4_of_ridges_kept": dict(true_l4),
        "examples_false": examples,
    }


def main():
    out = Path(sys.argv[1])
    results = [analyse(d) for d in sys.argv[2:]]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    for r in results:
        print(r["run"], r["counts"])
        print("   L4 of F3-false ridges:", r["l4_of_false_ridges"])
        print("   L4 of kept ridges:   ", r["l4_of_ridges_kept"])


if __name__ == "__main__":
    main()
