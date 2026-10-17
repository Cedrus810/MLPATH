"""Descriptor-level CPU audit probes; no PES or molecular acceptance claim.

Asserts current limitations, not repaired behavior. Writes a fresh evidence
directory only. SimpleNamespace inputs isolate the descriptor functions; they
are not production candidates and do not exercise ingestion or refinement.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pars/src"))

import numpy as np  # noqa: E402
from pars.config import ParsConfig  # noqa: E402
from pars.frontier.dominance import dominates  # noqa: E402
from pars.frontier.normalize import fit_normalization  # noqa: E402
from pars.frontier.objectives import trigger_severity  # noqa: E402
from pars.paths.cluster import build_tubes  # noqa: E402
from pars.paths.fingerprint import frechet_path_distance  # noqa: E402
from pars.provenance import implementation_sha256  # noqa: E402


def severity(amplitude, family, unit):
    candidate = SimpleNamespace(
        provenance=SimpleNamespace(
            perturbation_amplitude=amplitude,
            perturbation_family=family,
            perturbation_unit=unit,
        )
    )
    return trigger_severity(candidate)


def normalized_order(stretch_value, stretch_unit):
    raw = {
        "stretch": {"trigger_severity": severity(stretch_value, "stretch", stretch_unit)},
        "kick": {"trigger_severity": severity(0.2, "kick", "eV")},
    }
    fit = fit_normalization(raw, ParsConfig())
    normalized = {key: fit.normalize(value) for key, value in raw.items()}
    return {
        "raw": raw,
        "normalization": fit.to_dict(),
        "normalized": normalized,
        "stretch_dominates_kick": dominates(normalized["stretch"], normalized["kick"]),
        "kick_dominates_stretch": dominates(normalized["kick"], normalized["stretch"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    out = parser.parse_args().out.resolve()
    out.mkdir(parents=True, exist_ok=False)

    current_units = normalized_order(0.5, "A")
    equivalent_length_representation = normalized_order(0.05, "nm")
    assert current_units["kick_dominates_stretch"]
    assert equivalent_length_representation["stretch_dominates_kick"]
    assert severity(0.2, "stretch", "A") == severity(0.2, "kick", "eV")

    # Same synthetic recorded bond sequence, shared endpoints, distinct internal
    # geometry in the middle. No bond perception or physical saddle is asserted.
    start = np.array([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]])
    end = np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 0]])
    middle_a = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.5, 1.0, 0]])
    middle_b = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.5, 2.0, 0]])
    signatures = (
        frozenset({(0, 1), (1, 2)}),
        frozenset({(0, 1)}),
        frozenset({(0, 1), (0, 2)}),
    )
    candidates = [
        SimpleNamespace(
            candidate_id=name,
            signatures=signatures,
            positions=np.array([start, middle, end]),
        )
        for name, middle in (("path-a", middle_a), ("path-b", middle_b))
    ]
    distance = frechet_path_distance(candidates[0].signatures, candidates[1].signatures)
    tubes = build_tubes(candidates, "synthetic", ParsConfig().tube_distance_threshold)
    middle_pair_distances = [
        float(np.linalg.norm(middle[2] - middle[0])) for middle in (middle_a, middle_b)
    ]
    assert middle_pair_distances[0] != middle_pair_distances[1]
    assert distance == 0.0 and len(tubes) == 1

    report = {
        "workspace": str(ROOT),
        "pars_implementation_sha256": implementation_sha256(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "evidence_scope": "descriptor-only synthetic CPU probes; no PES, ingestion, archive deletion, certified saddle, historical contamination or GPU claim",
        "cross_family_severity": {
            "current_unit_example": current_units,
            "equivalent_length_unit_thought_experiment": equivalent_length_representation,
            "equal_numeric_A_and_eV_are_indistinguishable": True,
            "caveat": "PRRS currently emits fixed units, not nm; changing to nm illustrates missing dimensional comparability, not a reproduced producer unit-change bug. Other objectives are omitted/equal in this probe.",
        },
        "same_signature_distinct_geometry": {
            "path_positions_A": {c.candidate_id: c.positions.tolist() for c in candidates},
            "middle_mapped_pair_0_2_distances_A": middle_pair_distances,
            "signature_distance": distance,
            "tube_count": len(tubes),
            "tube_members": tubes[0].member_ids,
            "caveat": "The descriptor cannot separate these inputs; actual physical multi-saddle loss requires an independently certified reference pool and full selection/refinement comparison.",
        },
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
