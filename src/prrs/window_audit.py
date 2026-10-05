"""Matched response-horizon replay; never certifies a static basin boundary.

The baseline must reproduce every selected endpoint before ANY alternate
horizon runs. Trials use the original microstate/probe/thermal seed, with
only response_steps changed. Unresolved samples cannot certify stability.
"""

import hashlib
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from ase.io import read

from .config import SearchConfig
from .io import atomic_json
from .perturbations import Probe
from .provenance import calculator_identity


def window_semantics(config):
    return {
        "rule": "finite_horizon_response_then_quench",
        "response_steps": config.response_steps,
        "timestep_fs": config.timestep_fs,
        "response_time_fs": config.response_steps * config.timestep_fs,
        "horizon_stability": "not_measured",
        "static_basin_boundary": False,
        "meaning": "amplitude windows depend on the response horizon, thermal seed and quench; neither a static basin boundary nor a committor",
    }


def response_windows(labels):
    ordered = sorted((float(a), label) for a, label in labels.items())
    return [
        [a, b]
        for (a, left), (b, right) in zip(ordered, ordered[1:])
        if left == "returned" and right.startswith("departed:")
    ]


def assess_horizons(by_steps, required_steps):
    required = sorted(set(required_steps))
    if len(required) < 2 or any(type(s) is not int or s <= 0 for s in required):
        raise ValueError("at least two distinct positive response horizons required")
    missing = [s for s in required if str(s) not in by_steps]
    if missing:
        return {
            "status": "incomplete",
            "missing_steps": missing,
            "stable_on_tested_grid": False,
        }
    labels = [by_steps[str(s)]["labels"] for s in required]
    if not labels[0] or any(set(v) != set(labels[0]) for v in labels):
        return {"status": "incomplete", "stable_on_tested_grid": False}
    unresolved = any(
        label.startswith("unresolved:") for row in labels for label in row.values()
    )
    moved = any(row != labels[0] for row in labels[1:])
    return {
        "status": "unresolved"
        if unresolved
        else "sensitive"
        if moved
        else "stable_on_tested_grid",
        "stable_on_tested_grid": not unresolved and not moved,
        "label_changes": moved,
        "static_basin_boundary": False,
    }


def _file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replay_window_audit(run, output, factory, metadata, response_steps, scans=None):
    from .network import Registry
    from .runner import run_trial
    from .search import _stamped_with_charge, versions

    run, output = Path(run).resolve(), Path(output).resolve()
    if output == run or output.is_relative_to(run):
        raise ValueError("audit output must be outside the source run")
    if output.exists():
        raise FileExistsError(f"will not reuse {output}")
    names = ("config.json", "network.json", "attempts.jsonl")
    source_hashes = {name: _file_hash(run / name) for name in names}
    manifest = json.loads((run / "config.json").read_text())
    config = SearchConfig.from_dict(manifest["config"])
    horizons = sorted(set([config.response_steps, *response_steps]))
    assess_horizons({}, horizons)  # validate before any writes/calculations
    network = json.loads((run / "network.json").read_text())
    if network["status"] != "completed":
        raise ValueError("audit requires a completed source run")
    all_scans = network["amplitude_scans"]
    indices = list(range(len(all_scans))) if scans is None else sorted(set(scans))
    if not indices or any(i < 0 or i >= len(all_scans) for i in indices):
        raise ValueError("select at least one valid scan")
    identity = calculator_identity(factory, metadata)
    if not identity["verifiable"]:
        raise ValueError("audit calculator identity is unverifiable")
    # New manifests must agree; legacy manifests are replay-gated, never upgraded.
    if "calculator_identity" in manifest and identity != manifest["calculator_identity"]:
        raise ValueError("audit calculator differs from the source run")
    registry = Registry.from_network(network, run, config)
    states = {s.id: s for n in registry.nodes for s in n.microstates}
    records = {
        row["attempt_id"]: row
        for row in (
            json.loads(line)
            for line in (run / "attempts.jsonl").read_text().splitlines()
            if line.strip()
        )
    }
    selected = []
    for index in indices:
        scan = all_scans[index]
        if not scan["samples"]:
            raise ValueError(f"scan {index} has no samples")
        selected.append(
            {
                "scan": index,
                "source_microstate": scan["source_microstate"],
                "by_steps": {},
                "samples": scan["samples"],
            }
        )
    report = {
        "schema_version": 1,
        "source_run": str(run),
        "source_sha256": source_hashes,
        "source_implementation_sha256": manifest.get("implementation_sha256"),
        "audit_implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "calculator_identity": identity,
        "versions": versions(),
        "semantics": window_semantics(config),
        "horizons_steps": horizons,
        "horizons_fs": [s * config.timestep_fs for s in horizons],
        "baseline_reproduced": False,
        "scope": "selected_fixed_amplitude_samples",
        "scan_indices": indices,
        "status": "running",
        "scans": selected,
        "baseline_tolerance_A": 1e-5,
        "static_basin_boundary": False,
    }
    output.mkdir(parents=True)
    atomic_json(output / "report.json", report)
    started = time.perf_counter()
    try:
        # Finish the baseline across ALL scans before replaying alternatives.
        for steps in [
            config.response_steps,
            *[s for s in horizons if s != config.response_steps],
        ]:
            cfg = replace(config, response_steps=steps)
            for entry in selected:
                state = states[entry["source_microstate"]]
                labels = {}
                for sample in entry["samples"]:
                    published = records[sample["attempt"]]
                    probe = Probe(
                        **{
                            k: v
                            for k, v in published["probe"].items()
                            if k in Probe.__dataclass_fields__
                        }
                    )
                    trial_id = f"s{entry['scan']:04d}h{steps:06d}_{sample['attempt']}"
                    outcome = run_trial(
                        _stamped_with_charge(state.structure, cfg),
                        probe,
                        factory,
                        cfg,
                        output / "trials" / trial_id,
                        trial_id,
                    )
                    if steps == config.response_steps:
                        path = (
                            run / published["endpoint"] if published.get("endpoint") else None
                        )
                        expected = read(str(path)) if path and path.is_file() else None
                        if (expected is None) != (outcome.endpoint is None) or (
                            expected is not None
                            and (
                                not np.array_equal(expected.numbers, outcome.endpoint.numbers)
                                or not np.allclose(
                                    expected.positions,
                                    outcome.endpoint.positions,
                                    atol=1e-5,
                                    rtol=0,
                                )
                            )
                        ):
                            raise ValueError(
                                f"baseline endpoint mismatch: scan {entry['scan']} {sample['attempt']}"
                            )
                    label = f"unresolved:{outcome.record['status']}"
                    if outcome.record["status"] == "completed" and outcome.endpoint is not None:
                        _, landed = registry.identify(outcome.endpoint)
                        label = (
                            "unresolved:unidentified_endpoint"
                            if landed is None
                            else "returned"
                            if landed.id == state.id
                            else f"departed:{landed.id}"
                        )
                    labels[str(sample["amplitude"])] = label
                entry["by_steps"][str(steps)] = {
                    "labels": labels,
                    "windows": response_windows(labels),
                }
            if steps == config.response_steps:
                report["baseline_reproduced"] = True
            atomic_json(output / "report.json", report)
        for entry in selected:
            entry["assessment"] = assess_horizons(entry["by_steps"], horizons)
        report["status"] = (
            "stable_on_tested_grid"
            if all(e["assessment"]["stable_on_tested_grid"] for e in selected)
            else "not_certified"
        )
    except Exception as error:
        report.update(status="failed", failure=f"{type(error).__name__}: {error}")
        raise
    finally:
        report["wall_s"] = time.perf_counter() - started
        report["source_unchanged"] = all(
            _file_hash(run / n) == value for n, value in source_hashes.items()
        )
        atomic_json(output / "report.json", report)
    if not report["source_unchanged"]:
        raise ValueError("source run changed during audit")
    return report
