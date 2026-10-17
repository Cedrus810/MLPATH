"""Reproduce an accepted negative normalization scale; no PES or MD."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pars/src"))

from pars.config import ParsConfig  # noqa: E402
from pars.frontier.dominance import dominates  # noqa: E402
from pars.frontier.normalize import fit_normalization  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    output = parser.parse_args().out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = [
        "pars/src/pars/config.py",
        "pars/src/pars/frontier/normalize.py",
        "pars/src/pars/frontier/dominance.py",
        "benchmarks/project_audit/normalization_scale_probe.py",
    ]
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    config_input = {"objective_scales": [["energy_excursion", -1.0]]}
    config = ParsConfig.from_dict(config_input)
    raw = {f"A{i}": {"energy_excursion": 0.0} for i in range(8)}
    raw["B"] = {"energy_excursion": 1.0}
    norm = fit_normalization(raw, config)
    normalized = {cid: norm.normalize(vec) for cid, vec in raw.items()}
    assert norm.scale["energy_excursion"] == -1.0
    assert not dominates(raw["B"], raw["A0"])
    assert dominates(normalized["B"], normalized["A0"])
    report = {
        "kind": "negative-normalization-scale-audit",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": hashes,
        "limits": "Synthetic objective vectors via real config parser; no production candidate ingestion or physical saddle outcome.",
        "accepted_config": config_input,
        "raw": raw,
        "normalization": norm.to_dict(),
        "normalized": normalized,
        "worse_B_dominates_A_before": False,
        "worse_B_dominates_A_after": True,
    }
    assert hashes == {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
