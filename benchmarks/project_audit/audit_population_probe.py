"""Check the scope of the current PARS audit population, without a PES or MD."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pars/src"))

import numpy as np  # noqa: E402
from pars.config import ParsConfig  # noqa: E402
from pars.frontier.archive import ParetoArchive  # noqa: E402
from pars.schema import PathCandidate, Provenance  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    output = parser.parse_args().out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = [
        "pars/src/pars/frontier/archive.py",
        "pars/src/pars/frontier/cover.py",
        "pars/src/pars/frontier/margin.py",
        "pars/src/pars/frontier/objectives.py",
        "benchmarks/project_audit/audit_population_probe.py",
    ]
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    positions = np.repeat(np.array([[[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]]), 3, axis=0)
    candidates = [
        PathCandidate(
            candidate_id=cid,
            channel_id="synthetic",
            trace_id=cid,
            numbers=np.ones(3, dtype=int),
            positions=positions.copy(),
            energies=np.zeros(3),
            uncertainty=None,
            signatures=(frozenset({(0, 1)}),) * 3,
            provenance=Provenance(
                run_dir="synthetic-not-a-run",
                attempt_id=cid,
                channel_id="synthetic",
                event_key=None,
                reactant_state_id=None,
                product_state_id=None,
                frame_start=0,
                frame_stop=3,
            ),
        )
        for cid in ("a", "b", "c")
    ]
    archive = ParetoArchive("synthetic", candidates, ParsConfig(audit_size=5))
    unselected = sorted(set(archive.candidates) - set(archive.representatives))
    assert len(archive.frontier()) == 3
    assert len(archive.representatives) == 1
    assert len(unselected) == 2 and archive.audit_sample() == []
    report = {
        "kind": "selection-audit-population-probe",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": hashes,
        "limits": "Identical synthetic candidates; no physical diversity or missed saddle is asserted. Tests the archive population definition only.",
        "frontier": archive.frontier(),
        "representatives": archive.representatives,
        "unselected": unselected,
        "audit_sample": archive.audit_sample(),
        "selection": archive.selection_dict(),
    }
    assert hashes == {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in ("frontier", "representatives", "unselected", "audit_sample")
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
