"""Pure-geometry size count for discovery L0/L1 (PLAN_EXPLORE_V1 M3 Task 5).

No model is loaded. Counts conformers, poses and intents per system and extrapolates the
L1 cost from the b04 smoke (n_points=18, 8.7 forwards/point, 1.07 s/point). Counting only:
none of these geometries is a proposal.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from explore.intents import distance_intents, pair_intents  # noqa: E402
from explore.poses import (  # noqa: E402
    RULE_VERSION,
    conformer_grid,
    contact_ratio,
    fragments,
    sample_poses,
    single_fragment_poses,
)
from prrs.config import SearchConfig  # noqa: E402

STEPS_PER_TURN, CAP, SEED = 6, 10**6, 0
N_SAMPLES = (8, 32)
N_INTENT_POSES = 16
SHELL_DELTAS = (-0.4, 0.0, 0.4)
LOCALITIES = ("shared_atom", "adjacent", "all")
POINTS_PER_SCAN, FORWARDS_PER_POINT, SECONDS_PER_POINT = 18, 8.7, 1.07
COST_NOTE = (
    "b04 (6 atoms) measured, extrapolated; larger systems likely slower per point "
    "(more BFGS steps, slower forward); order of magnitude only"
)
PAIR_ASSUMPTION = (
    "ASSUMPTION: every pair intent is charged the same 18 points x 1.07 s as a single-distance "
    "scan; a pair intent may cost more"
)
GRID0_NOTE = (
    "Per-pose intent counts come from grid-0 conformer-pair poses only (the first 16 placed "
    "poses); est. poses/intents = placed poses x conformer pairs, i.e. they assume every "
    "conformer pair has the same per-pose count (matters for b09 x210, n1 x6)"
)
SYSTEMS = {
    "b04": dict(path="benchmarks/b04_sn2_chloride/input.extxyz", charge=-1, note=""),
    "b05": dict(path="benchmarks/b05_tertiary_halide/input.extxyz", charge=-1, note=""),
    "n1": dict(
        path=None,
        charge=0,
        note=(
            "MMFF geometry, counting only, not a POLAR-1 minimum; the contact ratio is an "
            "artefact of placing the MMFF fragments 8 A apart, not a complex contact ratio"
        ),
    ),
    "b09": dict(
        path="benchmarks/carbene_ambimodal/probe_b09/input.extxyz",
        charge=0,
        note="biased toward TS-2OMe; counting only, not a blind-search start",
    ),
}


def build_n1():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    parts = []
    for smi in ("C1=CC=CC1", "C=CC=O"):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=0)
        AllChem.MMFFOptimizeMolecule(mol)
        pos = mol.GetConformer().GetPositions()
        parts.append((pos - pos.mean(axis=0), [a.GetAtomicNum() for a in mol.GetAtoms()]))
    (pa, na), (pb, nb) = parts
    pb = pb + np.array([8.0, 0.0, 0.0])
    return Atoms(numbers=na + nb, positions=np.vstack([pa, pb]))


def load(name):
    spec = SYSTEMS[name]
    return build_n1() if spec["path"] is None else read(ROOT / spec["path"])


def _mean(xs):
    return float(np.mean(xs)) if xs else 0.0


def cost(n_intents):
    pts = n_intents * POINTS_PER_SCAN
    return {
        "forwards": pts * FORWARDS_PER_POINT,
        "wall_s": pts * SECONDS_PER_POINT,
        "gpu_hours_one_2080ti": pts * SECONDS_PER_POINT / 3600.0,
    }


def count_system(name, config, timing):
    atoms = load(name)
    frags = fragments(atoms, config)
    rec = {
        "system": name,
        "note": SYSTEMS[name]["note"],
        "charge": SYSTEMS[name]["charge"],
        "n_atoms": len(atoms),
        "n_fragments": len(frags),
        "fragments": [{"key": f.key, "n_atoms": len(f.atoms)} for f in frags],
    }
    t0 = time.time()
    confs, creps = [], []
    for f in frags:
        c, r = conformer_grid(atoms, f, config, STEPS_PER_TURN, CAP, SEED)
        confs.append(c)
        creps.append(r)
    rec["conformers"] = [
        {**r, "fragment": f.key, "generated": sum(c.status == "generated" for c in cs)}
        for f, cs, r in zip(frags, confs, creps)
    ]
    timing[f"{name}:conformers"] = time.time() - t0
    gen = [max(1, sum(c.status == "generated" for c in cs)) for cs in confs]
    pairs_of_conformers = int(np.prod(gen))
    rec["conformer_pairs"] = pairs_of_conformers
    rec["grid_input_conformer_status"] = [cs[0].status for cs in confs]

    rec["ratio"] = None
    configs = []  # (label, shell|None, n_samples|None, poses)
    if len(frags) == 2:
        fa, fb = frags
        ratio = contact_ratio(
            atoms.positions[list(fa.atoms)],
            atoms.numbers[list(fa.atoms)],
            atoms.positions[list(fb.atoms)],
            atoms.numbers[list(fb.atoms)],
        )
        rec["ratio"] = ratio
        for d in SHELL_DELTAS:
            shell = round(ratio + d, 6)
            if shell <= 0:
                continue
            label = (
                ("placement f (8 A artificial)" if name == "n1" else "measured f")
                if d == 0
                else f"f{d:+g} (sensitivity only, not a proposal)"
            )
            for n in N_SAMPLES:
                configs.append((label, shell, n, [(frags, confs[0][0], confs[1][0])]))
    elif len(frags) == 1:
        rec["pose_note"] = "single fragment: single_fragment_poses, one pose per conformer"
    else:
        rec["pose_note"] = f"{len(frags)} fragments: v1 poses skipped (needs 1 or 2)"

    rows = []
    t0 = time.time()
    if len(frags) == 2:
        for label, shell, n, ((fr, ca, cb),) in configs:
            poses = sample_poses(atoms, fr, ca, cb, [shell], n, SEED)
            rows.append(
                _row(atoms, config, label, shell, n, poses, pairs_of_conformers, n == 32)
            )
    elif len(frags) == 1:
        poses = single_fragment_poses(atoms, confs[0])
        rows.append(_row(atoms, config, "single fragment", None, None, poses, 1, True))
    timing[f"{name}:poses_intents"] = time.time() - t0
    rec["rows"] = rows
    return rec


def _row(atoms, config, label, shell, n, poses, conf_pairs, with_intents):
    status = {}
    for p in poses:
        status[p.status] = status.get(p.status, 0) + 1
    placed = [p for p in poses if p.status == "placed"]
    row = {
        "label": label,
        "shell": shell,
        "n_samples": n,
        "poses": len(poses),
        "placed": len(placed),
        "status_counts": status,
        "t_A_mean": _mean([p.t_A for p in placed if getattr(p, "t_A", None) is not None]),
    }
    if not with_intents or not placed:
        return row
    sample = placed[:N_INTENT_POSES]
    per = {"singles": [], "bonded": [], "near": [], "intra": [], "inter": [], "far_classes": []}
    per_pair = {k: [] for k in LOCALITIES}
    per_comb = {k: [] for k in LOCALITIES}
    for p in sample:
        a = atoms.copy()
        a.positions = p.positions
        singles, rep = distance_intents(a, p.pose_id, config, 0, SEED)
        per["singles"].append(len(singles))
        per["bonded"].append(sum(s.layer == "bonded" for s in singles))
        per["near"].append(sum(s.layer == "near" for s in singles))
        per["inter"].append(sum(s.scope == "inter" for s in singles))
        per["intra"].append(sum(s.scope != "inter" for s in singles))
        per["far_classes"].append(rep["far_classes_total"])
        for loc in LOCALITIES:
            pi, prep = pair_intents(singles, a, config, loc, p.pose_id)
            per_pair[loc].append(len(pi))
            per_comb[loc].append(prep["combinations_by_locality"][loc])
    mean = {k: _mean(v) for k, v in per.items()}
    row["intents_poses_used"] = len(sample)
    row["per_pose_mean"] = {
        **mean,
        "pair_intents": {k: _mean(v) for k, v in per_pair.items()},
        "pair_combinations": {k: _mean(v) for k, v in per_comb.items()},
    }
    row["far_classes_total_per_pose_mean"] = mean["far_classes"]
    n_poses_total = len(placed) * conf_pairs
    row["total_poses_est"] = n_poses_total
    tot = {"single": mean["singles"] * n_poses_total}
    for loc in LOCALITIES:
        tot[f"single+pair[{loc}]"] = (mean["singles"] + _mean(per_pair[loc])) * n_poses_total
    row["total_intents_est"] = tot
    row["cost_extrapolated"] = {k: cost(v) for k, v in tot.items()}
    return row


def fmt(x, nd=1):
    return f"{x:,.{nd}f}"


def markdown(results, params):
    L = ["# L0 pure-geometry counts", ""]
    L += [
        f"Rule version: `{RULE_VERSION}`. Pure geometry, no model loaded. Counting only; no geometry here is a proposal.",
        "",
    ]
    L += ["Parameters:", ""]
    L += [f"- {k}: {v}" for k, v in params.items()]
    L += [
        "",
        f"Cost note: {COST_NOTE}.",
        "",
        f"Note: {GRID0_NOTE}.",
        "",
        f"Note: {PAIR_ASSUMPTION}.",
        "",
    ]
    summ = []
    for r in results:
        L += [f"## {r['system']}", ""]
        L += [f"- atoms {r['n_atoms']}, charge {r['charge']}, fragments {r['n_fragments']}"]
        if r["note"]:
            L += [f"- NOTE: {r['note']}"]
        if r.get("ratio") is not None:
            what = (
                "placement f (8 A artificial, not a complex contact ratio)"
                if r["system"] == "n1"
                else "measured contact ratio f of the input"
            )
            L += [f"- {what}: {r['ratio']:.4f}"]
        lb = [c for c in r["conformers"] if c["rotors_not_enumerated"]]
        if lb:
            L += [
                "- LOWER BOUND: conformer counts exclude non-enumerated rotors "
                + ", ".join(
                    f"{c['rotors_not_enumerated']} in {c['fragment'].split('|')[0]}"
                    for c in lb
                    if c["rotors_not_enumerated"]
                )
                + "; true conformer counts are larger (ring conformations are never enumerated in v1)"
            ]
        if r.get("pose_note"):
            L += [f"- {r['pose_note']}"]
        L += ["", "Conformers (steps_per_turn 6, cap 1e6)", ""]
        L += [
            "| fragment | genuine torsions | rotors not enumerated | grid | generated | status counts |",
            "|---|---|---|---|---|---|",
        ]
        for c in r["conformers"]:
            L += [
                f"| {c['fragment']} | {c['genuine_torsions']} | {c['rotors_not_enumerated']} | {c['grid_size']} | {c['generated']} | {c['status_counts']} |"
            ]
        L += ["", f"Conformer pairs (generated x generated): {r['conformer_pairs']}", ""]
        if r["rows"]:
            L += [
                "| shell | n_samples | poses | placed | status | singles/pose (bonded, near; intra, inter) [grid-0 poses] | pair intents/pose (shared_atom, adjacent, all) [grid-0 poses] | far classes/pose [grid-0 poses] | est. poses (placed x conformer pairs) | est. intents single (assumes equal per conformer pair) | single+pair adjacent (same assumption) |",
                "|---|---|---|---|---|---|---|---|---|---|---|",
            ]
        for w in r["rows"]:
            sh = "-" if w["shell"] is None else f"{w['shell']:.3f} ({w['label']})"
            if "per_pose_mean" not in w:
                L += [
                    f"| {sh} | {w['n_samples']} | {w['poses']} | {w['placed']} | {w['status_counts']} | n/a | n/a | n/a | n/a | n/a | n/a |"
                ]
                continue
            m = w["per_pose_mean"]
            pp = m["pair_intents"]
            L += [
                f"| {sh} | {w['n_samples']} | {w['poses']} | {w['placed']} | {w['status_counts']} "
                f"| {fmt(m['singles'])} ({fmt(m['bonded'])}, {fmt(m['near'])}; {fmt(m['intra'])}, {fmt(m['inter'])}) "
                f"| {fmt(pp['shared_atom'], 0)}, {fmt(pp['adjacent'], 0)}, {fmt(pp['all'], 0)} "
                f"| {fmt(w['far_classes_total_per_pose_mean'])} | {fmt(w['total_poses_est'], 0)} "
                f"| {fmt(w['total_intents_est']['single'], 0)} | {fmt(w['total_intents_est']['single+pair[adjacent]'], 0)} |"
            ]
            summ.append((r["system"], w))
        L += [""]
        costed = [w for w in r["rows"] if "cost_extrapolated" in w]
        if costed:
            L += ["Cost extrapolation (" + COST_NOTE + ")", ""]
            L += [
                "| shell | n_samples | rule | intents | forwards | wall hours | GPU-hours (1x 2080 Ti) |",
                "|---|---|---|---|---|---|---|",
            ]
            for w in costed:
                for k, v in w["total_intents_est"].items():
                    c = w["cost_extrapolated"][k]
                    L += [
                        f"| {w['label']} {w['shell']:.3f} | {w['n_samples']} | {k} | {fmt(v, 0)} | {c['forwards']:.3e} | {c['wall_s'] / 3600:.2f} | {c['gpu_hours_one_2080ti']:.2f} |"
                    ]
            L += [f"{PAIR_ASSUMPTION}.", "", f"{GRID0_NOTE}.", ""]
    L += [
        "## Summary (n_samples 32, measured-f and sensitivity shells; per-pose intents from first 16 placed poses)",
        "",
    ]
    L += [
        "| system | note | shell | placed poses | singles/pose | single+pair(adjacent)/pose | est. intents single | est. intents single+adjacent | GPU-h single | GPU-h single+adjacent |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    notes = {r["system"]: r["note"] or "-" for r in results}
    for r in results:
        if any(c["rotors_not_enumerated"] for c in r["conformers"]):
            notes[r["system"]] += "; conformer counts are lower bounds (rotors not enumerated)"
    for s, w in summ:
        if w["n_samples"] != 32:
            continue
        m = w["per_pose_mean"]
        c1 = w["cost_extrapolated"]["single"]["gpu_hours_one_2080ti"]
        c2 = w["cost_extrapolated"]["single+pair[adjacent]"]["gpu_hours_one_2080ti"]
        sh = (
            "-"
            if w["shell"] is None
            else f"{w['shell']:.3f} ({w['label'].replace(' (sensitivity only, not a proposal)', ', sensitivity only')})"
        )
        L += [
            f"| {s} | {notes[s]} | {sh} | {w['placed']} | {fmt(m['singles'])} | {fmt(m['singles'] + m['pair_intents']['adjacent'])} | {fmt(w['total_intents_est']['single'], 0)} | {fmt(w['total_intents_est']['single+pair[adjacent]'], 0)} | {c1:.2f} | {c2:.2f} |"
        ]
    L += ["", f"Cost note: {COST_NOTE}.", "", f"{PAIR_ASSUMPTION}.", "", f"{GRID0_NOTE}.", ""]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--systems", nargs="+", default=list(SYSTEMS), choices=list(SYSTEMS))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    timing, results, ratios = {}, [], {}
    t_start = time.time()
    for name in args.systems:
        config = SearchConfig(
            charge_sensitive=True, total_charge=SYSTEMS[name]["charge"], multiplicity=1
        )
        t0 = time.time()
        rec = count_system(name, config, timing)
        timing[name] = time.time() - t0
        print(f"{name}: {timing[name]:.1f} s", flush=True)
        results.append(rec)
        ratios[name] = {
            "contact_ratio": rec["ratio"],
            "n_fragments": rec["n_fragments"],
            "note": rec["note"],
        }
    params = {
        "steps_per_turn": STEPS_PER_TURN,
        "conformer cap": CAP,
        "seed": SEED,
        "n_samples": list(N_SAMPLES),
        "shells": "measured f of the input, f-0.4, f+0.4 (the +-0.4 are sensitivity only, not proposals)",
        "intent poses per (system, shell, n_samples=32)": N_INTENT_POSES,
        "n_far": 0,
        "pair locality rules": list(LOCALITIES),
        "scan points per intent": POINTS_PER_SCAN,
        "forwards per point": FORWARDS_PER_POINT,
        "seconds per point": SECONDS_PER_POINT,
        "systems": args.systems,
        "grid0_note": GRID0_NOTE,
        "pair_cost_assumption": PAIR_ASSUMPTION,
    }
    (out / "contact_ratios.json").write_text(json.dumps(ratios, indent=2))
    (out / "counts.json").write_text(
        json.dumps(
            {
                "rule_version": RULE_VERSION,
                "params": params,
                "cost_note": COST_NOTE,
                "timing_s": timing,
                "systems": results,
            },
            indent=2,
            default=float,
        )
    )
    (out / "COUNTS.md").write_text(markdown(results, params))
    print(f"total {time.time() - t_start:.1f} s -> {out}")


if __name__ == "__main__":
    main()
