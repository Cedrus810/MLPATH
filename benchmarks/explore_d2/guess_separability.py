"""Offline test: does the L4 starting guess tell "same saddle" apart from "different saddle"?

Follow-up to repeat_response_census.py (docs/TODO.md item 2). There the L3 release endpoint,
grouped per pose, missed most after-the-fact duplicates: certified events span 4-6 poses.
Here the feature is the certification guess itself (frozen rule M5-D1, `_promotion_points`),
compared across poses by symmetric RMSD over the graph automorphisms.

Ground truth exists only for certified events: members of one event (c2 + its duplicate_c2)
are positive pairs, members of different events negative pairs. All promoted items enter
the sweep: for each RMSD radius r, greedy leader clustering in journal order (an item joins
the first earlier leader within r), "first item per cluster only" gives the forwards skipped
and whether every event keeps at least one member. No radius is chosen here; the sweep only
reports what each radius would cost.

Writes the RMSD matrix, guesses, scan tangents and guess-to-own-pose RMSD per run to
guess_<run>.npz next to OUT_JSON, so follow-up reads need no recomputation.

Reads journals only; no model, no evaluation, no writes to the run directory.

HF runs (rule D2-HF-v1): the items are the HS tasks (same promotion guess), an item's
outcome is its HC status when HC ran (else the HS status), its cost HS + HC forwards, and
events are HC `connected` with the `hf_duplicate` items that point at them.

    python guess_separability.py OUT_JSON PROCS RUN_DIR [RUN_DIR ...]
"""

import json
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from ase import Atoms  # noqa: E402
from explore.d2 import _classify, _points, _promotion_points, _task  # noqa: E402
from prrs.chemistry import automorphisms, chemical_key, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

_W = {}


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue


def _init(numbers, group, active, guesses):
    _W.update(atoms=Atoms(numbers=numbers), group=group, active=active, guesses=guesses)


def _row(i):
    a = _W["atoms"].copy()
    a.positions = _W["guesses"][i]
    out = []
    for j in range(i + 1, len(_W["guesses"])):
        b = _W["atoms"].copy()
        b.positions = _W["guesses"][j]
        out.append(symmetric_rmsd(a, b, _W["group"], _W["active"]))
    return out


def analyse(run_dir, procs, out_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    c = SearchConfig.from_dict(manifest["config"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        atoms = Atoms(numbers=r["numbers"], positions=r["positions"])
        atoms.info.update(r["info"])
        poses[r["pose_id"]] = atoms
    rows = {r["intent_id"]: r for r in _lines(run_dir / "population.jsonl")}
    hf = manifest["rule_version"].startswith("D2-HF")
    head = "HS|" if hf else "L4|"
    settled, l4 = {}, []
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") != "settle" or ev["task"][:3] not in ("L1|", "L2|", "L4|", "HS|", "HC|"):
            continue
        if ev["task"].startswith(head) and ev["task"] not in settled:
            l4.append(ev["task"])
        if ev["task"][:3] in ("HS|", "HC|"):  # search histories are not needed here
            p = ev["payload"]
            ev["payload"] = {"class": p.get("class"), "duplicate_of": p.get("duplicate_of")}
        settled[ev["task"]] = ev
    run = SimpleNamespace(settled=settled, settings=SimpleNamespace(noise_eV=manifest["settings"]["noise_eV"]))

    keys = {
        chemical_key(
            p,
            c.bond_scale,
            c.active_atoms,
            charge_sensitive=c.charge_sensitive,
            charge=c.total_charge,
            multiplicity=c.multiplicity,
            parity_tolerance=c.parity_tolerance,
        )["key"]
        for p in poses.values()
    }
    numbers = {tuple(p.numbers) for p in poses.values()}
    assert len(numbers) == 1, "poses differ in atom order"
    template = next(iter(poses.values()))
    group, _ = automorphisms(template, c.bond_scale, c.active_atoms, c.automorphism_limit)

    guesses, tangents, class_mismatch = [], [], 0
    for task in l4:
        row = rows[task[3:]]
        classified = _classify(run, row)
        if classified["class"] != settled[task]["payload"].get("class"):
            class_mismatch += 1
        points = _points(settled[_task("L1", row)]["payload"])
        _, guess, tangent = _promotion_points(classified, points)
        guesses.append(np.asarray(guess, float))
        tangents.append(np.asarray(tangent, float))

    t0 = time.time()
    with Pool(procs, initializer=_init, initargs=(template.numbers, group, c.active_atoms, guesses)) as pool:
        tri = pool.map(_row, range(len(l4)), chunksize=1)
    n = len(l4)
    d = np.zeros((n, n))
    for i, vals in enumerate(tri):
        d[i, i + 1 :] = vals
        d[i + 1 :, i] = vals
    wall = time.time() - t0

    def outcome(t):
        hc = settled.get("HC|" + t[3:]) if hf else None
        return hc or settled[t]

    status = [outcome(t)["status"] for t in l4]
    spent = np.array(
        [settled[t]["spent"] + (outcome(t)["spent"] if outcome(t) is not settled[t] else 0) for t in l4],
        float,
    )
    pose = [rows[t[3:]]["pose_id"] for t in l4]
    found, dup = ("connected", "hf_duplicate") if hf else ("c2", "duplicate_c2")
    event_of = {}
    for i, t in enumerate(l4):
        if status[i] == found:
            event_of[i] = ("HC|" if hf else "L4|") + t[3:]
    for i, t in enumerate(l4):
        if status[i] == dup:
            event_of[i] = outcome(t)["payload"]["duplicate_of"]
    members = list(event_of)
    pos, neg = [], []
    for a in range(len(members)):
        for b in range(a + 1, len(members)):
            i, j = members[a], members[b]
            (pos if event_of[i] == event_of[j] else neg).append(
                {"d": d[i, j], "same_pose": pose[i] == pose[j]}
            )
    off = d[np.triu_indices(n, 1)]
    to_pose = []
    for i, t in enumerate(l4):
        a, b = template.copy(), poses[pose[i]].copy()
        a.positions = guesses[i]
        to_pose.append(symmetric_rmsd(b, a, group, c.active_atoms))
    np.savez_compressed(
        out_dir / f"guess_{run_dir.name}.npz",
        rmsd=d,
        guesses=np.array(guesses),
        tangents=np.array(tangents),
        to_pose=np.array(to_pose),
        spent=spent,
    )
    items = [
        {
            "task": t,
            "status": status[i],
            "spent": float(spent[i]),
            "class": settled[t]["payload"].get("class"),
            "pose": pose[i],
            "event": event_of.get(i),
            "guess_to_pose_rmsd_A": float(to_pose[i]),
        }
        for i, t in enumerate(l4)
    ]

    sweep = []
    for r in sorted({0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5}):
        leaders, skipped = [], []
        for i in range(n):
            if any(d[i, k] <= r for k in leaders):
                skipped.append(i)
            else:
                leaders.append(i)
        kept = set(leaders)
        events = set(event_of.values())
        sweep.append(
            {
                "radius_A": r,
                "clusters": len(leaders),
                "items_skipped": len(skipped),
                "forwards_skipped": float(spent[skipped].sum()),
                "forwards_skipped_fraction": float(spent[skipped].sum() / spent.sum()),
                "skipped_status": dict(Counter(status[i] for i in skipped)),
                "events_kept": sum(any(event_of[i] == e for i in kept if i in event_of) for e in events),
                "events_total": len(events),
            }
        )

    def q(xs):
        xs = np.asarray(xs, float)
        return (
            None
            if xs.size == 0
            else {k: float(np.quantile(xs, v)) for k, v in (("min", 0), ("p50", 0.5), ("p90", 0.9), ("max", 1))}
        )

    return {
        "run": str(run_dir),
        "rule_version": manifest["rule_version"],
        "pose_chemical_keys": len(keys),
        "automorphisms": len(group),
        "l4_items": n,
        "class_mismatch": class_mismatch,
        "pairwise_wall_s": wall,
        "all_pairs_rmsd_A": q(off),
        "positive_pairs": {"n": len(pos), "rmsd_A": q([p["d"] for p in pos]), "same_pose": sum(p["same_pose"] for p in pos)},
        "negative_pairs": {"n": len(neg), "rmsd_A": q([p["d"] for p in neg])},
        "sweep": sweep,
        "items": items,
    }


def main():
    out_path, procs, run_dirs = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3:]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results = [analyse(r, procs, out_path.parent) for r in run_dirs]
    out_path.write_text(json.dumps(results, indent=1), encoding="utf-8")
    print(json.dumps([{k: v for k, v in r.items() if k != "items"} for r in results], indent=1))


if __name__ == "__main__":
    main()
