"""Read a bounded b09 evidence snapshot on the writer host; no model imports.

Prints JSON to stdout. Reads only this checkout's named runs/results and the
immediate b09 run names. Does not simulate, alter runs, or certify their PES.
Transfer this script via SSH stdin when NFS may cache newly created files.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/ruigengji/MLPATH")


def main():
    files = {}

    def read(relative, json_lines=False):
        path = ROOT / relative
        before = path.stat()
        data = path.read_bytes()
        after = path.stat()
        files[relative] = {
            "sha256": hashlib.sha256(data).hexdigest(),
            "bytes": len(data),
            "mtime_ns": after.st_mtime_ns,
            "stable_during_read": (before.st_size, before.st_mtime_ns)
            == (after.st_size, after.st_mtime_ns),
        }
        return (
            [json.loads(line) for line in data.splitlines() if line.strip()]
            if json_lines
            else json.loads(data)
        )

    source_hash = hashlib.sha256()
    for path in sorted((ROOT / "src/prrs").glob("*.py")):
        source_hash.update(path.name.encode())
        source_hash.update(path.read_bytes())
    current_hash = source_hash.hexdigest()
    runs = {}
    for name in ("b09_probe1", "b09_compress_coverage", "b09_coarse_smoke", "b09_smoke0"):
        relative = f"runs/{name}"
        network = read(f"{relative}/network.json")
        manifest = read(f"{relative}/config.json")
        attempts = read(f"{relative}/attempts.jsonl", json_lines=True)
        config = manifest.get("config", {})
        ids = [record.get("attempt_id") for record in attempts]
        ledger = network.get("direction_ledger")
        runs[name] = {
            "persisted_status": network.get("status"),
            "termination": network.get("termination"),
            "attempt_count": len(attempts),
            "network_attempt_count": network.get("attempts"),
            "duplicate_attempt_ids": [key for key, count in Counter(ids).items() if count > 1],
            "outcome_counts_from_attempts": dict(Counter(r.get("status") for r in attempts)),
            "classification_counts_from_attempts": dict(
                Counter(str(r.get("classification")) for r in attempts)
            ),
            "phase_counts_from_attempts": dict(Counter(str(r.get("stage")) for r in attempts)),
            "reaction_count": len(network.get("reactions", [])),
            "reaction_channel_count": len(network.get("reaction_channels", [])),
            "ts_candidate_count": len(network.get("ts_candidates", [])),
            "ts_connected_channel_count": len(network.get("ts_connected_channels", [])),
            "chemical_node_count": len(network.get("chemical_nodes", [])),
            "direction_ledger": ledger,
            "manifest_host": manifest.get("host"),
            "run_implementation_sha256": manifest.get("implementation_sha256"),
            "run_source_matches_writer_current": manifest.get("implementation_sha256")
            == current_hash,
            "config": {
                key: config.get(key)
                for key in (
                    "families",
                    "direction_mode",
                    "geometry_amplitudes_A",
                    "max_trials",
                    "quench_steps",
                    "temperature_K",
                    "response_steps",
                    "timestep_fs",
                    "committor_shots",
                    "max_depth",
                )
            },
        }

    result_base = "benchmarks/carbene_ambimodal/results/2026-10-10"
    replays = read(f"{result_base}/replay_quench_coverage/summary.json")
    replay_files = [
        read(f"{result_base}/replay_quench_coverage/{row['label']}.json") for row in replays
    ]
    analysis = read(f"{result_base}/b09_compress_coverage_analysis.json")
    proc = subprocess.run(
        ["ps", "-p", "953421,953745,956356,957756,958049", "-o", "pid=,stat=,etime=,comm="],
        capture_output=True,
        text=True,
    )
    snapshot = {
        "snapshot_utc": datetime.now(timezone.utc).isoformat(),
        "writer_host": platform.node(),
        "workspace": str(ROOT),
        "evidence_scope": "read-only writer-host snapshot; no fresh molecular calculation, physical certification or whole-project run inventory",
        "reader_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        if Path(__file__).is_file()
        else None,
        "writer_current_prrs_sha256": current_hash,
        "immediate_b09_run_names": sorted(
            p.name for p in (ROOT / "runs").iterdir() if p.is_dir() and p.name.startswith("b09")
        ),
        "runs": runs,
        "coverage_replays": {
            "count": len(replays),
            "summary_matches_per_trial_files": replays == replay_files,
            "converged": sum(bool(row.get("converged")) for row in replays),
            "minimum_confirmed": sum(bool(row.get("minimum")) for row in replays),
            "classification_counts": dict(
                Counter(str(r.get("classification")) for r in replays)
            ),
            "unconverged": [row for row in replays if not row.get("converged")],
            "converged_with_bond_changes": [
                row
                for row in replays
                if row.get("converged") and (row.get("formed") or row.get("broken"))
            ],
        },
        "p001_replay": read(f"{result_base}/replay_quench_smoke/summary.json"),
        "stored_analysis": {
            "trial_count": len(analysis.get("trials", [])),
            "ts_candidate_count": len(analysis.get("ts_candidates", [])),
            "ts_candidates": analysis.get("ts_candidates"),
        },
        "historical_pid_check": {
            "pids": [953421, 953745, 956356, 957756, 958049],
            "ps_exit_code": proc.returncode,
            "ps_stdout": proc.stdout,
            "ps_stderr": proc.stderr,
            "caveat": "checks recorded PIDs only, not all possible b09 jobs; persisted running status is not process liveness",
        },
        "file_provenance": files,
        "all_files_stable_during_individual_reads": all(
            entry["stable_during_read"] for entry in files.values()
        ),
    }
    print(json.dumps(snapshot, indent=2))


if __name__ == "__main__":
    main()
