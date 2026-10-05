"""M0 preflight: read-only inventory of the maturation layer's named inputs.

PLAN_GEN_COMPAS_MATURATION_2026-10-03 section 4 (M0): read EXACTLY the named
runs -- no recursion over runs/ -- and check, per run, the PARS contract,
the all-atom paths, the A/B region and its mapped microstates, the potential
and the file hashes. The verdict per run is ``supported`` /
``supported_reference_only`` / ``unsupported`` WITH the reasons; confirm and
explore batches are classified separately. The one file this script writes
is ``--out``.

Region identity follows plan M0: A/B are the mapped microstates and their
bond graphs, never the symmetry-quotiented chemical_key -- for P1 both
sides share one chemical_key (the degenerate self-loop), so the region is
defined by microstate ids plus their encoded graphs, and first-entry
shooting is eligible exactly when the two graphs differ.

Usage (repository root, science environment):

    PYTHONPATH=src python docs/experiments/maturation_preflight.py \
        --out docs/experiments/maturation_preflight/preflight.json \
        [--baseline-json docs/experiments/maturation_preflight/baseline.json]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from ase.io import read as ase_read  # noqa: E402

from prrs.config import SearchConfig  # noqa: E402
from prrs.state import encode  # noqa: E402

DEFAULT_RUNS = {
    "source_run": REPO / "runs" / "p1_committor2",
    "confirm_batch": REPO / "runs" / "p1_first_entry_confirm",
    "explore_batch": REPO / "runs" / "p1_path_committor_explore",
    "refinement_baseline": REPO / "runs" / "pars_ground_truth",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prrs_tree_fingerprint() -> str:
    """Same algorithm as search.py's manifest check / implementation_hash."""
    source_hash = hashlib.sha256()
    for path in sorted((REPO / "src" / "prrs").glob("*.py")):
        source_hash.update(path.name.encode())
        source_hash.update(path.read_bytes())
    return source_hash.hexdigest()


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


# ---- the four named inputs --------------------------------------------------


def check_source_run(run_dir: Path) -> dict:
    """runs/p1_committor2: the run whose channels, TS and microstates feed M2/M3."""
    evidence, caveats, failures = {}, [], []
    config = _read_json(run_dir / "config.json")
    evidence["implementation_sha256"] = config["implementation_sha256"]
    evidence["input_structure_sha256"] = config.get("input_structure_sha256")
    evidence["host"] = config.get("host")
    evidence["versions"] = config.get("versions")
    calculator = config["calculator"]
    evidence["calculator"] = {
        "factory": calculator["factory"],
        "kwargs": calculator["kwargs"],
    }
    for model_path, recorded in calculator.get("local_file_sha256", {}).items():
        model = Path(model_path)
        if not model.exists():
            failures.append(f"model checkpoint missing: {model_path}")
            continue
        actual = sha256_file(model)
        evidence["calculator"]["checkpoint_sha256_ok"] = actual == recorded
        if actual != recorded:
            failures.append(f"model checkpoint hash mismatch: {model_path}")

    network = _read_json(run_dir / "network.json")
    evidence["network"] = {
        "schema_version": network.get("schema_version"),
        "status": network.get("status"),
        "attempts": network.get("attempts"),
        "chemical_nodes": len(network.get("chemical_nodes") or []),
        "reaction_channels": len(network.get("reaction_channels") or []),
        "ts_candidates": len(network.get("ts_candidates") or []),
        "product_block_present": "product" in network,
    }
    if network.get("attempts") and not isinstance(network["attempts"], int):
        failures.append("network.attempts is not a count")
    channels = network.get("reaction_channels") or []
    evidence["channels"] = [
        {
            "id": c["id"],
            "class": c["class"],
            "self_loop": c["self_loop"],
            "event_key": c["event_key"],
            "broken": c["broken"],
            "formed": c["formed"],
            "ends": c["ends"],
        }
        for c in channels
    ]
    if not channels:
        failures.append("no reaction_channels: PARS fail-closed reader would refuse this run")

    node = (network.get("chemical_nodes") or [None])[0]
    microstates = (node or {}).get("microstates") or []
    evidence["microstates"] = [
        {
            "id": m["id"],
            "energy_eV": m["energy_eV"],
            "structure": m["structure"],
            "structure_present": (run_dir / m["structure"]).exists(),
        }
        for m in microstates
    ]
    if not all(m["structure_present"] for m in evidence["microstates"]):
        failures.append("a microstate structure file is missing from chemical/")

    trials_dir = run_dir / "trials"
    trial_dirs = sorted(d for d in trials_dir.iterdir() if d.is_dir()) if trials_dir.is_dir() else []
    with_traj = [d for d in trial_dirs if (d / "response.traj").exists()]
    evidence["trials"] = {
        "dirs": len(trial_dirs),
        "with_response_traj": len(with_traj),
    }
    input_atoms = ase_read(str(run_dir / "input.extxyz"))
    evidence["input"] = {"n_atoms": len(input_atoms), "formula": input_atoms.get_chemical_formula()}
    if with_traj:
        frames = ase_read(str(with_traj[0] / "response.traj"), index=":")
        all_atom = len(frames[0]) == len(input_atoms)
        evidence["trials"]["sample_dir"] = with_traj[0].name
        evidence["trials_sample"] = {
            "n_frames": len(frames),
            "n_atoms": len(frames[0]),
            "all_atom": all_atom,
        }
        if not all_atom:
            failures.append("response trajectory is not all-atom (hydrogens lost?)")
    else:
        failures.append("no trial carries a response trajectory")

    paths_dir = run_dir / "paths"
    evidence["paths"] = (
        sorted(p.name for p in paths_dir.iterdir()) if paths_dir.is_dir() else []
    )
    evidence["ts_files"] = (
        sorted(p.name for p in (run_dir / "ts_candidates").iterdir())
        if (run_dir / "ts_candidates").is_dir()
        else []
    )
    if evidence["network"]["ts_candidates"] == 0:
        failures.append("no ts_candidates: nothing for a transition-region diagnosis to locate")

    if config["implementation_sha256"] != prrs_tree_fingerprint():
        caveats.append(
            "the source run's fingerprint is historical provenance; M3 shooting runs "
            "under the CURRENT tree, whose fingerprint is recorded at batch initialize"
        )
    verdict = "supported" if not failures else "unsupported"
    return {"role": "source_run", "verdict": verdict, "evidence": evidence, "failures": failures, "caveats": caveats}


def check_confirm_batch(run_dir: Path) -> dict:
    """runs/p1_first_entry_confirm: the frozen-protocol confirm batch (N=50 x 4 windows)."""
    evidence, failures = {}, []
    batch = _read_json(run_dir / "batch.json")
    verdict = _read_json(run_dir / "verdict.json")
    evidence["protocol"] = batch["protocol"]
    evidence["protocol_sha256"] = batch["protocol_sha256"]
    evidence["source_run"] = batch["source_run"]
    evidence["source_fingerprint"] = batch["source_fingerprint"]
    evidence["N_per_window"] = batch["N"]
    evidence["first_entry_k"] = batch["k"]
    evidence["L_max_run"] = batch["L_max_run"]
    evidence["seed_offset"] = batch["seed_offset"]
    evidence["thresholds"] = batch["thresholds"]
    evidence["V1_whole_batch"] = verdict.get("V1")
    evidence["windows"] = [
        {
            "scan": w["scan"],
            "A": w["A"],
            "B": w["B"],
            "counts": w["counts"],
            "p_B_given_AB": w["p_B_given_AB"],
            "wilson_95": w["wilson_95"],
        }
        for w in verdict.get("windows", [])
    ]
    trials = run_dir / "trials"
    n_trials = len([d for d in trials.iterdir() if d.is_dir()]) if trials.is_dir() else 0
    evidence["trial_dirs"] = n_trials
    expected = batch["N"] * len(evidence["windows"])
    if n_trials != expected:
        failures.append(f"trial dirs {n_trials} != N x windows ({expected})")
    # the region is the UNORDERED pair {A, B}: a window shot from the other
    # side legitimately swaps the labels (window 8 below does), and that
    # convention is recorded, not treated as a second region
    unordered = {frozenset((w["A"], w["B"])) for w in evidence["windows"]}
    if len(unordered) != 1:
        failures.append(
            f"windows span {len(unordered)} distinct unordered (A, B) regions; expected one"
        )
    else:
        evidence["region_pair_unordered"] = sorted(next(iter(unordered)))
        inverted = [
            w["scan"]
            for w in evidence["windows"]
            if (w["A"], w["B"]) != (evidence["region_pair_unordered"][0],
                                    evidence["region_pair_unordered"][1])
        ]
        evidence["windows_with_swapped_labels"] = inverted
    if batch["source_run"] != "runs/p1_committor2":
        failures.append("confirm batch does not name runs/p1_committor2 as its source")
    return {
        "role": "confirm_batch_reference",
        "verdict": "supported" if not failures else "unsupported",
        "evidence": evidence,
        "failures": failures,
        "caveats": [],
    }


def check_explore_batch(run_dir: Path) -> dict:
    """runs/p1_path_committor_explore: N=16 per frame, self-declared not citable."""
    evidence, caveats = {}, []
    summary = _read_json(run_dir / "summary.json")
    evidence["meaning"] = summary["meaning"]
    evidence["source_run"] = summary["source_run"]
    evidence["saddles"] = [
        {
            "ts": s["ts"],
            "A": s["A"],
            "B": s["B"],
            "frames": [
                {"frame": f["frame"], "A": f["A"], "B": f["B"], "other": f["other"],
                 "p_B_given_AB": f["p_B_given_AB"]}
                for f in s["frames"]
            ],
        }
        for s in summary["saddles"]
    ]
    for ts in summary["saddles"]:
        if not (run_dir / ts["ts"]).is_dir():
            caveats.append(f"{ts['ts']} directory missing on disk")
    caveats.append(
        "exploration batch (N=16): usable to LOCATE the next batch, never citable as "
        "a result -- the file says so itself"
    )
    return {
        "role": "explore_batch_not_citable",
        "verdict": "supported_reference_only",
        "evidence": evidence,
        "failures": [],
        "caveats": caveats,
    }


def check_refinement_baseline(run_dir: Path) -> dict:
    """runs/pars_ground_truth: the 27 refined-candidate records (PARS_REVIEW_2026-10-02)."""
    evidence, failures = {}, []
    records = []
    for sub in sorted(p for p in run_dir.iterdir() if p.is_dir()):
        for record_path in sorted(sub.glob("*.json")):
            record = _read_json(record_path)
            records.append((sub.name, record_path.name, record))
    evidence["n_records"] = len(records)
    by_sub = {}
    barriers = []
    sources = set()
    for sub, name, record in records:
        by_sub.setdefault(sub, 0)
        by_sub[sub] += 1
        if record.get("barrier_eV") is not None:
            barriers.append(record["barrier_eV"])
        source = (record.get("provenance") or {}).get("source") or {}
        sources.add(
            json.dumps(
                {
                    "run_dir": source.get("run_dir"),
                    "channel_id": source.get("channel_id"),
                    "event_key": source.get("event_key"),
                },
                sort_keys=True,
            )
        )
    evidence["records_per_source_run"] = by_sub
    evidence["distinct_sources"] = [json.loads(s) for s in sorted(sources)]
    evidence["backend"] = sorted({r.get("backend") for _, _, r in records})
    status_counts = {}
    for _, _, r in records:
        status_counts[r.get("status")] = status_counts.get(r.get("status"), 0) + 1
    evidence["status_counts"] = status_counts
    if barriers:
        evidence["barrier_eV"] = {
            "n": len(barriers),
            "min": min(barriers),
            "max": max(barriers),
            "mean": sum(barriers) / len(barriers),
        }
    # the published baseline: 27/27 refined to the same ~0.419 eV saddle
    if len(records) != 27:
        failures.append(f"expected 27 refinement records, found {len(records)}")
    if barriers and (max(barriers) - min(barriers)) > 1e-3:
        failures.append(
            "refined barriers span more than 1 meV: the 'one saddle' baseline does not hold"
        )
    # the source runs the candidates point at must exist, with their fingerprints
    fingerprints = {}
    for source in evidence["distinct_sources"]:
        source_run = REPO / source["run_dir"]
        config_path = source_run / "config.json"
        if not config_path.exists():
            failures.append(f"baseline source run missing: {source['run_dir']}")
            continue
        fingerprints[source["run_dir"]] = _read_json(config_path)["implementation_sha256"]
    evidence["source_run_fingerprints"] = fingerprints
    caveats = [
        "refinement records postdate the 2026-10-02 PARS fixes (PARS_REVIEW_2026-10-02.md); "
        "they are the baseline, not a new run of this layer"
    ]
    return {
        "role": "refinement_baseline",
        "verdict": "supported" if not failures else "unsupported",
        "evidence": evidence,
        "failures": failures,
        "caveats": caveats,
    }


def check_region_eligibility(run_dir: Path) -> dict:
    """First-entry eligibility of the P1 region, from mapped microstates.

    Both sides share one chemical_key (the degenerate self-loop), so the
    region is the two microstate structures and their encoded graphs. The
    confirm batch names the pair; here the graphs are actually encoded and
    compared.
    """
    report = {"method": "mapped_microstate_graphs", "same_chemical_key": None, "graphs": {}}
    confirm = DEFAULT_RUNS["confirm_batch"]
    verdict = _read_json(confirm / "verdict.json")
    windows = verdict.get("windows") or []
    if not windows:
        return {"error": "no confirm windows name a region"}
    a_id, b_id = windows[0]["A"], windows[0]["B"]
    report["sides"] = {"A": a_id, "B": b_id}
    network = _read_json(run_dir / "network.json")
    node = network["chemical_nodes"][0]
    states = {m["id"]: m for m in node["microstates"]}
    if a_id not in states or b_id not in states:
        return {"error": f"confirm windows name microstates absent from the source run: {a_id}, {b_id}"}
    report["same_chemical_key"] = a_id.split("/")[0] == b_id.split("/")[0]
    config = SearchConfig.from_dict(_read_json(run_dir / "config.json")["config"])
    for side, state_id in (("A", a_id), ("B", b_id)):
        atoms = ase_read(str(run_dir / states[state_id]["structure"]))
        edges = frozenset(encode(atoms, config.bond_scale, active=config.active_atoms).edges)
        report["graphs"][side] = {
            "n_edges": len(edges),
            "graph_sha256": hashlib.sha256(
                json.dumps(sorted([list(e) for e in edges])).encode()
            ).hexdigest(),
        }
    report["first_entry_eligible"] = (
        report["graphs"]["A"]["graph_sha256"] != report["graphs"]["B"]["graph_sha256"]
    )
    report["meaning"] = (
        "eligible: first-entry shooting can tell the sides apart"
        if report["first_entry_eligible"]
        else "unsupported_same_graph: first-entry cannot define this region"
    )
    return report


# ---- assembly ---------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out",
        default=str(REPO / "docs" / "experiments" / "maturation_preflight" / "preflight.json"),
    )
    parser.add_argument("--baseline-json", default=None, help="pytest/digest baseline to embed")
    args = parser.parse_args()

    source_report = check_source_run(DEFAULT_RUNS["source_run"])
    reports = {
        "runs/p1_committor2": source_report,
        "runs/p1_first_entry_confirm": check_confirm_batch(DEFAULT_RUNS["confirm_batch"]),
        "runs/p1_path_committor_explore": check_explore_batch(DEFAULT_RUNS["explore_batch"]),
        "runs/pars_ground_truth": check_refinement_baseline(DEFAULT_RUNS["refinement_baseline"]),
    }
    selection_status = {
        "pars_selection_present": any(
            (DEFAULT_RUNS["source_run"] / name).exists()
            for name in ("pars_selection.json", "selection")
        ),
        "action_required": (
            "no pars selection directory exists for the P1 source run: M2 must run "
            "'pars select' into a NEW output directory under the fixed 2026-10-02 "
            "implementation and record its hash in the protocol; the original run "
            "directory stays untouched"
        ),
    }
    preflight = {
        "kind": "mlpath-maturation-preflight",
        "m0_of_plan": "PLAN_GEN_COMPAS_MATURATION_2026-10-03.md section 4",
        "generated_utc": utc_now(),
        "script": "docs/experiments/maturation_preflight.py",
        "script_sha256": sha256_file(Path(__file__)),
        "current_prrs_implementation_sha256": prrs_tree_fingerprint(),
        "runs": reports,
        "region": check_region_eligibility(DEFAULT_RUNS["source_run"]),
        "pars_selection": selection_status,
        "notes": [
            "only the four named inputs were read; no recursion over runs/",
            "all checks are read-only; the only file written is --out",
            "confirm and explore batches are classified separately and the explore "
            "batch keeps its not-citable marking",
        ],
    }
    if args.baseline_json:
        preflight["baseline"] = _read_json(Path(args.baseline_json))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(preflight, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    verdicts = {name: report["verdict"] for name, report in reports.items()}
    print(json.dumps({"out": str(out), "verdicts": verdicts}, indent=2))
    return 0 if all(v == "supported" for v in verdicts.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
