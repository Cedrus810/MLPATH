"""Offline census of repeated responses among promoted D2 items (finished runs, Hessian path).

Question: before the expensive L4 certify, how many promoted items are "the same response"
as an earlier one, how many certify forwards a representatives-only policy would skip, and
whether every certified event (c2 and its duplicate_c2 members) keeps a representative.

Proxy for "the same response": the L3 release endpoint. Two promoted items share a group
when they start from the same pose and their release endpoints fall in the same basin under
the existing PRRS basin rule (same chemical key, |dE| <= basin_energy_eV, symmetric RMSD <=
basin_rmsd_A over the graph automorphisms; poses.py qualify_conformers). No new threshold.
Groups are built greedily in journal order (each endpoint is compared with the leaders of
earlier groups), i.e. what an online rule would have seen. Items whose release produced no
endpoint are grouped by (pose, release status) and reported separately: they are not
evidence of a shared response.

Reads manifest / poses / population / journal only; no model, no evaluation, no writes to
the run directory.

    python repeat_response_census.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from ase import Atoms  # noqa: E402
from prrs.chemistry import automorphisms, chemical_key, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:  # torn tail of a live journal
                    continue


def load(run_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    config = SearchConfig.from_dict(manifest["config"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        atoms = Atoms(numbers=r["numbers"], positions=r["positions"])
        atoms.info.update(r["info"])
        poses[r["pose_id"]] = atoms
    pose_of = {r["intent_id"]: r["pose_id"] for r in _lines(run_dir / "population.jsonl")}
    settled, order = {}, []
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") == "settle" and ev["task"][:3] in ("L3|", "L4|"):
            if ev["task"] not in settled:
                order.append(ev["task"])
            settled[ev["task"]] = ev
    return manifest, config, poses, pose_of, settled, order


def census(run_dir):
    manifest, c, poses, pose_of, settled, order = load(run_dir)
    l4 = [t for t in order if t.startswith("L4|")]
    groups, leaders, auts = [], defaultdict(list), {}
    group_of, no_endpoint = {}, Counter()
    for task in l4:
        intent = task[3:]
        pose = pose_of[intent]
        rel = settled.get("L3|" + intent)
        end = rel and rel["payload"].get("endpoint_positions")
        if end is None:
            status = rel["status"] if rel else "no_L3"
            no_endpoint[status] += 1
            key = ("no_endpoint", pose, status)
            if key not in leaders:
                leaders[key] = [len(groups)]
                groups.append({"kind": "no_endpoint", "release_status": status, "members": []})
            gi = leaders[key][0]
        else:
            atoms = poses[pose].copy()
            atoms.positions = np.asarray(end, float)
            energy = rel["payload"]["endpoint_energy_eV"]
            ckey = chemical_key(
                atoms,
                c.bond_scale,
                c.active_atoms,
                charge_sensitive=c.charge_sensitive,
                charge=c.total_charge,
                multiplicity=c.multiplicity,
                parity_tolerance=c.parity_tolerance,
            )["key"]
            if ckey not in auts:
                auts[ckey], _ = automorphisms(atoms, c.bond_scale, c.active_atoms, c.automorphism_limit)
            gi = None
            for cand in leaders[("endpoint", pose, ckey)]:
                g = groups[cand]
                if abs(g["energy_eV"] - energy) > c.basin_energy_eV:
                    continue
                if symmetric_rmsd(atoms, g["atoms"], auts[ckey], c.active_atoms) <= c.basin_rmsd_A:
                    gi = cand
                    break
            if gi is None:
                gi = len(groups)
                leaders[("endpoint", pose, ckey)].append(gi)
                groups.append(
                    {
                        "kind": "endpoint",
                        "atoms": atoms,
                        "energy_eV": energy,
                        "returns_to_pose": ckey
                        == chemical_key(
                            poses[pose],
                            c.bond_scale,
                            c.active_atoms,
                            charge_sensitive=c.charge_sensitive,
                            charge=c.total_charge,
                            multiplicity=c.multiplicity,
                            parity_tolerance=c.parity_tolerance,
                        )["key"],
                        "release_label": rel["payload"].get("label"),
                        "members": [],
                    }
                )
        ev = settled[task]
        groups[gi]["members"].append(
            {
                "task": task,
                "status": ev["status"],
                "spent": ev["spent"],
                "duplicate_of": ev["payload"].get("duplicate_of"),
                "class": ev["payload"].get("class"),
            }
        )
        group_of[task] = gi

    # certified events: each c2 with the duplicate_c2 items that point at it
    events = {t: {t} for t in l4 if settled[t]["status"] == "c2"}
    for t in l4:
        dup = settled[t]["payload"].get("duplicate_of")
        if settled[t]["status"] == "duplicate_c2" and dup in events:
            events[dup].add(t)

    total = sum(settled[t]["spent"] for t in l4)
    out = {
        "run": str(run_dir),
        "rule_version": manifest["rule_version"],
        "basin_rule": {"basin_energy_eV": c.basin_energy_eV, "basin_rmsd_A": c.basin_rmsd_A},
        "l4_items": len(l4),
        "l4_forwards": total,
        "l4_status": dict(Counter(settled[t]["status"] for t in l4)),
        "duplicate_c2_forwards": sum(
            settled[t]["spent"] for t in l4 if settled[t]["status"] == "duplicate_c2"
        ),
        "no_endpoint_items": dict(no_endpoint),
        "policies": {},
        "events": [],
    }
    endpoint_groups = [g for g in groups if g["kind"] == "endpoint"]
    sizes = Counter(len(g["members"]) for g in endpoint_groups)
    out["endpoint_groups"] = {
        "groups": len(endpoint_groups),
        "items": sum(len(g["members"]) for g in endpoint_groups),
        "size_histogram": dict(sorted(sizes.items())),
        "groups_returning_to_pose_graph": sum(g["returns_to_pose"] for g in endpoint_groups),
        "largest": [
            {
                "size": len(g["members"]),
                "returns_to_pose_graph": g["returns_to_pose"],
                "release_label": g["release_label"],
                "l4_status": dict(Counter(m["status"] for m in g["members"])),
                "forwards": sum(m["spent"] for m in g["members"]),
            }
            for g in sorted(endpoint_groups, key=lambda g: -len(g["members"]))[:10]
        ],
    }
    for k in (1, 2):
        kept = {m["task"] for g in endpoint_groups for m in g["members"][:k]}
        kept |= {m["task"] for g in groups if g["kind"] == "no_endpoint" for m in g["members"]}
        skipped = [t for t in l4 if t not in kept]
        out["policies"][f"first_{k}_per_endpoint_group"] = {
            "l4_items_run": len(kept),
            "l4_items_skipped": len(skipped),
            "forwards_skipped": sum(settled[t]["spent"] for t in skipped),
            "forwards_skipped_fraction": sum(settled[t]["spent"] for t in skipped) / total
            if total
            else None,
            "skipped_status": dict(Counter(settled[t]["status"] for t in skipped)),
            "events_kept": sum(bool(members & kept) for members in events.values()),
            "events_total": len(events),
        }
    for c2, members in events.items():
        out["events"].append(
            {
                "c2": c2,
                "members": len(members),
                "groups": sorted({group_of[t] for t in members}),
                "group_kinds": sorted({groups[group_of[t]]["kind"] for t in members}),
                "fingerprint": {
                    k: settled[c2]["payload"].get("event", {}).get("fingerprint", {}).get(k)
                    for k in ("broken", "formed", "aligned_rmsd_A", "energy_change_eV")
                },
            }
        )
    # are after-the-fact duplicates already in their original's endpoint group?
    same = diff = 0
    for c2, members in events.items():
        for t in members - {c2}:
            if group_of[t] == group_of[c2] and groups[group_of[t]]["kind"] == "endpoint":
                same += 1
            else:
                diff += 1
    out["duplicate_c2_in_original_endpoint_group"] = {"same": same, "different": diff}
    return out


def main():
    out_path, run_dirs = Path(sys.argv[1]), sys.argv[2:]
    results = [census(d) for d in run_dirs]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=1), encoding="utf-8")
    for r in results:
        print(
            json.dumps(
                {k: r[k] for k in r if k not in ("events",)} | {"events_n": len(r["events"])},
                indent=1,
            )
        )


if __name__ == "__main__":
    main()
