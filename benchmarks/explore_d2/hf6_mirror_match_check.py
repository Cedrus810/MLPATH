"""Read-only: do HF-6 reference matches change when mirror images count (2026-11-12)?

The HF-6 match (hf_event_validate.match_group + symmetric_rmsd) allows atom permutations and
proper rotations only. From an achiral reactant, the mirror image of a saddle is a saddle of
the same energy; the N1 references keep one enantiomer per class (build.py). Here every
order-1 candidate the HF-6 summaries scored is matched again with its mirror image (z -> -z)
also allowed, using the energies the summaries recorded. No model evaluation.

    python hf6_mirror_match_check.py OUT_JSON
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_event_validate import match_group, validation_config  # noqa: E402  (sets sys.path)

from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

from prrs.chemistry import symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

import hf6_e2e_summary as e2e  # noqa: E402
import hf6_local_summary as loc  # noqa: E402

RES = loc.RES


def ref_atoms(spec, run, numbers):
    if "refs" in spec:
        return [
            (p.parent.name if p.name == "geometry.extxyz" else p.stem, read(p).positions)
            for p in spec["refs"]
        ]
    return [
        (f"event{k}", np.asarray(ev["positions"]))
        for k, ev in enumerate(e2e.hc_connected(spec.get("ref_events")))
    ]


def matcher(run, numbers, refs_pos, ref_names, ref_energy, config):
    """refs_pos order = the summaries' reference order; names taken from the local summary."""
    refs = []
    for (_, pos), name, e in zip(refs_pos, ref_names, ref_energy):
        a = Atoms(numbers=numbers, positions=pos)
        refs.append((name, a, e, match_group(a, run, config)))

    def match(pos, energy):
        best = {}
        for mirror in (False, True):
            q = np.asarray(pos) * (np.array([1, 1, -1]) if mirror else 1)
            a = Atoms(numbers=numbers, positions=q)
            for name, r, er, g in refs:
                d = float(symmetric_rmsd(a, r, g, config.active_atoms))
                if d <= config.basin_rmsd_A and abs(energy - er) <= config.basin_energy_eV:
                    best.setdefault(mirror, name)
        return best.get(False), best.get(True)

    return match


def main():
    out = Path(sys.argv[1])
    report = {"local": {}, "e2e": {}}
    local = json.loads((RES / "2026-11-10/hf6_local_summary_v3.json").read_text())
    for name, spec in loc.SYSTEMS.items():
        run = spec["run"]
        config = validation_config(
            SearchConfig.from_dict(json.loads((run / "manifest.json").read_text())["config"])
        )
        numbers = json.loads((run / "poses.jsonl").read_text().splitlines()[0])["numbers"]
        refs_pos = ref_atoms({"refs": spec["refs"]} if "refs" in spec else {"ref_events": spec["ref_events"]}, run, numbers)
        refE = [r["energy_eV"] for r in local[name]["references"]]
        refN = [r["name"] for r in local[name]["references"]]
        match = matcher(run, numbers, refs_pos, refN, refE, config)
        rows = []
        for row in local[name]["rows"]:
            for arm in ("hf", "dimer"):
                c = row[arm]
                if not c.get("candidate") or c.get("order") != 1:
                    continue
                rec = json.loads((spec["result"] / f"{row['item']:03d}_{arm}.json").read_text())
                old = c.get("match")
                plain, mirror = match(rec["positions"], c["energy_eV"])
                if plain != old:
                    raise RuntimeError(f"{name} {row['item']} {arm}: replay {plain} != {old}")
                if mirror and not plain:
                    rows.append({"item": row["item"], "arm": arm, "mirror_match": mirror})
        report["local"][name] = rows
    for tag in ("sA", "sB"):
        summ = json.loads(
            (RES / "2026-11-11" / ("hf6_e2e_summary_sA_v3.json" if tag == "sA" else "hf6_e2e_summary_sB_v2.json")).read_text()
        )
        report["e2e"][tag] = {}
        for name, spec in e2e.systems(tag).items():
            run = spec["d2"]
            config = validation_config(
                SearchConfig.from_dict(json.loads((run / "manifest.json").read_text())["config"])
            )
            numbers = json.loads((run / "poses.jsonl").read_text().splitlines()[0])["numbers"]
            refs_pos = ref_atoms(spec, run, numbers)
            arms = summ[name]["arms"]
            # same references and fp32 single points as the local v3 summary; e2e names differ
            refE = [r["energy_eV"] for r in local[name]["references"]]
            refN = list(summ[name]["references"])
            match = matcher(run, numbers, refs_pos, refN, refE, config)
            pos = {"d2": [ev["positions"] for ev in e2e.hc_connected(run)], "db": []}
            for f in sorted(spec["db"].glob("[0-9][0-9][0-9]_[0-9][0-9][0-9].json")):
                rec = json.loads(f.read_text())
                if rec["status"] == "converged":
                    pos["db"].append((f.stem, rec["positions"]))
            rows = []
            for arm in ("d2", "db"):
                detail = arms[arm]["candidates_detail"]
                plist = pos[arm] if arm == "d2" else [p for _, p in pos[arm]]
                assert len(plist) == len(detail), (tag, name, arm)
                for k, (p, c) in enumerate(zip(plist, detail)):
                    if c.get("order") != 1:
                        continue
                    plain, mirror = match(p, c["energy_eV"])
                    if plain != c.get("match"):
                        raise RuntimeError(f"{tag} {name} {arm} {k}: replay {plain} != {c.get('match')}")
                    if mirror and not plain:
                        rows.append({"arm": arm, "k": k, "cost": c["cost"], "mirror_match": mirror})
            report["e2e"][tag][name] = rows
    out.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
