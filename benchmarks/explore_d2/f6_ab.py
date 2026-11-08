"""F6 A/B: b05 D2-v1 fp32 rerun vs D2-v1.1 (scan tangent passed to the refiner), same seeds.

Reads journals only; no model, no evaluation. Pairs L4 certifications by intent and reports
the status transition, the saddles reached (c2 + duplicate_c2), and the L4 cost / wall time.

    python f6_ab.py OUT_JSON OLD_RUN NEW_RUN [OLD_RUN NEW_RUN ...]
"""

import json
import sys
from collections import Counter
from pathlib import Path


def l4(run_dir):
    out = {}
    times = []
    with (Path(run_dir) / "journal.jsonl").open() as fh:
        for line in fh:
            if '"L4|' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle":
                out[ev["task"][3:]] = ev
                times.append(ev["t"])
    return out, times


def _median(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else None


def summary(run_dir):
    rows, times = l4(run_dir)
    sad = [
        (ev["payload"]["saddle_energy_eV"], ev["status"], ev["payload"].get("duplicate_of"))
        for ev in rows.values()
        if ev["status"] in ("c2", "duplicate_c2")
    ]
    reasons = Counter()
    overlaps = []
    for ev in rows.values():
        sr = (ev["payload"].get("diagnostics") or {}).get("saddle_refine") or {}
        if sr.get("reason") == "lost_target_mode":
            overlaps.append((sr["target_overlap"], sr["mode_overlap"], sr["steps"]))
        reasons[ev["status"]] += 1
    labels = Counter(
        (ev["payload"].get("event") or {}).get("label")
        for ev in rows.values()
        if ev["status"] == "c2"
    )
    return rows, {
        "certified": len(rows),
        "forward": sum(ev["spent"] for ev in rows.values()),
        "wall_h": (max(times) - min(times)) / 3600 if times else 0,
        "status": dict(reasons),
        "c2_labels": dict(labels),
        "saddle_energies": Counter(round(e, 3) for e, _, _ in sad).most_common(),
        "lost_target_mode": {
            "n": len(overlaps),
            "target_overlap_none": sum(o[0] is None for o in overlaps),
            "target_overlap_median": _median([o[0] for o in overlaps if o[0] is not None]),
            "steps_median": _median([o[2] for o in overlaps]),
        },
    }


def main():
    out = Path(sys.argv[1])
    pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
    results = []
    for old, new in pairs:
        a, sa = summary(old)
        b, sb = summary(new)
        trans = Counter()
        for intent in set(a) | set(b):
            trans[(a[intent]["status"] if intent in a else "-", b[intent]["status"] if intent in b else "-")] += 1
        old_c2 = {i: a[i]["status"] for i in a if a[i]["status"] in ("c2", "duplicate_c2")}
        results.append(
            {
                "old": old,
                "new": new,
                "old_summary": sa,
                "new_summary": sb,
                "old_c2_now": {i.split("|")[-1]: [s, b[i]["status"] if i in b else "-"] for i, s in old_c2.items()},
                "transitions": sorted(([o, n, c] for (o, n), c in trans.items()), key=lambda r: -r[2]),
            }
        )
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    for r in results:
        print("==", r["old"], "->", r["new"])
        for k in ("certified", "forward", "wall_h", "c2_labels", "saddle_energies", "lost_target_mode"):
            print(f"  {k}: {r['old_summary'][k]}  ->  {r['new_summary'][k]}")
        print("  old C2 now:", r["old_c2_now"])
        for o, n, c in r["transitions"][:25]:
            print(f"   {c:4d}  {o:50s} -> {n}")


if __name__ == "__main__":
    main()
