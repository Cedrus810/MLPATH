"""Offline replay of the "same-graph ridge" predicate on finished D2 runs (read-only; no model
evaluation). For docs/HF_GENERAL_RULES_OPTIONS_2026-11-12.md §1.

For every ridge (`ridge_suspected`) HS the run settled or skipped, rebuild `read_scan` on the
run's own L1 points with the run's frozen noise_eV and evaluate three predicates:
  adjacent  current D2-HF-v1.2 `_same_graph_ridge`: edges(peak) == edges(peak + 1)
  basin     edges(scan start) == edges(lowest point of the converged run after the peak)
            (the point that supplies read_scan's "falls"); compares the basins on both sides
            of the ridge instead of two neighbouring scan steps
  window    no edge change anywhere from the scan start to the end of that converged run
Reported per predicate: ridges flagged same-graph, and how many of them gave an HS candidate.

    python same_graph_rule_replay.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import collections
import json
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src")]

from explore.probe import read_scan  # noqa: E402


def edges(p):
    return tuple(sorted(tuple(e) for e in p["edges"]))


def predicates(points, noise):
    pts = [
        SimpleNamespace(
            energy_eV=p["energy_eV"],
            converged=p["converged"],
            reaction_force=p["reaction_force"],
            edges=edges(p),
        )
        for p in points
    ]
    ridge = [t for t in read_scan(pts, noise) if t["kind"] == "ridge_suspected"]
    if not ridge:
        return None
    k = ridge[0]["peak"]
    end = k + 1
    while end + 1 < len(pts) and pts[end + 1].converged:
        end += 1
    low = min(range(k + 1, end + 1), key=lambda i: pts[i].energy_eV)
    return {
        "peak": k,
        "adjacent": pts[k].edges == pts[k + 1].edges,
        "basin": pts[0].edges == pts[low].edges,
        "window": all(pts[i].edges == pts[0].edges for i in range(end + 1)),
    }


def replay(run):
    noise = json.loads((run / "manifest.json").read_text())["settings"]["noise_eV"]
    l1, hs = {}, {}
    with open(run / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"L1|' not in line and '"HS|' not in line:
                continue
            ev = json.loads(line)
            task = ev.get("task", "")
            if task.startswith("L1|") and ev.get("event") == "settle":
                l1[task[3:]] = ev["payload"]["points"]
            elif task.startswith("HS|"):
                if ev.get("event") == "settle":
                    p = ev["payload"]
                    if p.get("class") == "ridge_suspected":
                        hs[task[3:]] = {
                            "status": ev["status"],
                            "flag": bool(p.get("same_graph_ridge")),
                        }
                elif ev.get("event") == "skip" and "same_graph_ridge" in ev.get("reason", ""):
                    hs[task[3:]] = {"status": ev["reason"], "flag": True}
    rows = []
    for key, h in hs.items():
        pr = predicates(l1[key], noise) if key in l1 else None
        rows.append({"key": key} | h | {"pred": pr})
    table = {}
    for name in ("adjacent", "basin", "window"):
        c = collections.Counter()
        for r in rows:
            if r["pred"] is None:
                c["unreplayed"] += 1
                continue
            same = r["pred"][name]
            ran = not r["status"].startswith("HS_")
            cand = r["status"] == "converged_candidate"
            tag = "same" if same else "changed"
            c[f"{tag}:ridges"] += 1
            c[f"{tag}:ran"] += ran
            c[f"{tag}:candidates"] += cand
        table[name] = dict(sorted(c.items()))
    check = sum(
        1 for r in rows if r["pred"] is not None and r["pred"]["adjacent"] != r["flag"]
    )
    return {
        "run": str(run.relative_to(ROOT)),
        "noise_eV": noise,
        "ridge_hs": len(rows),
        "adjacent_vs_recorded_flag_mismatch": check,
        "table": table,
        "rows": rows,
    }


def main():
    out = Path(sys.argv[1])
    res = [replay(Path(p).resolve()) for p in sys.argv[2:]]
    out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    for r in res:
        print(r["run"], r["ridge_hs"], "mismatch", r["adjacent_vs_recorded_flag_mismatch"])
        for k, v in r["table"].items():
            print("  ", k, v)


if __name__ == "__main__":
    main()
