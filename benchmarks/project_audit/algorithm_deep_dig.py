"""Synthetic algorithm-stage probes; no calculator, MD, or molecular acceptance.

The assertions pin limitations in the source at audit time. They are evidence
probes, not a regression contract for future fixes. Inputs bypass production
ingestion and do not establish physical basins or saddles.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "pars/src")]

import numpy as np  # noqa: E402
from pars.config import ParsConfig  # noqa: E402
from pars.frontier.margin import assign_tiers  # noqa: E402
from pars.frontier.normalize import fit_normalization  # noqa: E402
from pars.frontier.objectives import evaluate_candidate  # noqa: E402
from pars.paths.cluster import build_tubes  # noqa: E402
from pars.paths.fingerprint import frechet_path_distance, recrossing_count  # noqa: E402
from pars.paths.resample import resample_positions  # noqa: E402
from pars.schema import PathCandidate, Provenance  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.network import ChemicalNode, Microstate, Registry  # noqa: E402


def candidate(cid, signatures, positions=None):
    if positions is None:
        size = max(3, 1 + max((max(pair) for s in signatures for pair in s), default=0))
        frame = np.zeros((size, 3))
        frame[:, 0] = np.arange(size, dtype=float)
        frame[2] = [0.0, 1.0, 0.0]
        positions = np.repeat(frame[None], len(signatures), axis=0)
    return PathCandidate(
        candidate_id=cid,
        channel_id="synthetic",
        trace_id=cid,
        numbers=np.ones(positions.shape[1], dtype=int),
        positions=positions,
        energies=np.zeros(len(positions)),
        uncertainty=None,
        signatures=signatures,
        provenance=Provenance(
            run_dir="synthetic-not-a-run",
            attempt_id=cid,
            channel_id="synthetic",
            event_key=None,
            reactant_state_id=None,
            product_state_id=None,
            frame_start=0,
            frame_stop=len(positions),
        ),
    )


def spectator_dilution():
    cases = []
    for count in (0, 4, 32):
        common = frozenset((3 + 2 * i, 4 + 2 * i) for i in range(count))
        old, new = frozenset({(0, 1)}), frozenset({(1, 2)})
        first = (common | old, common, common | new)
        second = (common | old, common | old | new, common | new)
        distance = frechet_path_distance(first, second)
        tubes = build_tubes([candidate("a", first), candidate("b", second)], "synthetic", 0.3)
        cases.append({"unchanged_bonds": count, "distance": distance, "tubes": len(tubes)})
    assert cases[0]["tubes"] == 2 and cases[-1]["tubes"] == 1
    return cases


def leader_order():
    empty = frozenset()
    signatures = {
        "A": (empty, frozenset({(0, 1)}), empty),
        "B": (empty, frozenset({(0, 1), (1, 2)}), empty),
        "C": (empty, frozenset({(1, 2)}), empty),
    }
    cases = []
    for ordering in ("ABC", "BAC"):
        labels = {f"id{i}": label for i, label in enumerate(ordering)}
        tubes = build_tubes(
            [candidate(cid, signatures[label]) for cid, label in labels.items()],
            "synthetic",
            0.5,
        )
        cases.append(
            {
                "id_order": ordering,
                "members_by_original_label": [
                    [labels[cid] for cid in tube.member_ids] for tube in tubes
                ],
            }
        )
    assert [len(c["members_by_original_label"]) for c in cases] == [2, 1]
    return cases


def pool_normalization():
    baseline = {"A": {"energy_excursion": 0.0}, "B": {"energy_excursion": 0.1}}
    expanded = {
        **baseline,
        **{f"dominated{i}": {"energy_excursion": 100.0 + i} for i in range(20)},
    }
    cases = []
    for name, raw in (("original", baseline), ("extra_dominated", expanded)):
        config = ParsConfig()
        fitted = fit_normalization(raw, config)
        vectors = {cid: fitted.normalize(value) for cid, value in raw.items()}
        tiers, margins = assign_tiers(vectors, config)
        cases.append(
            {
                "pool": name,
                "scale": fitted.scale,
                "A_tier": tiers["A"].value,
                "B_tier": tiers["B"].value,
                "B_margin": margins["B"],
            }
        )
    assert cases[0]["A_tier"] == cases[1]["A_tier"]
    assert cases[0]["B_tier"] != cases[1]["B_tier"]
    return cases


def source_seed_collisions():
    path = ROOT / "src/prrs/search.py"
    tree = ast.parse(path.read_text())
    cls = next(
        n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "ReactionSearch"
    )
    assignments = [
        n
        for n in ast.walk(cls)
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "direction_seed" for t in n.targets)
    ]
    assert len(assignments) == 1
    expression = assignments[0].value
    code = compile(ast.Expression(expression), str(path), "eval")

    def seed(node_id, state_id):
        return eval(  # noqa: S307 - exact repository AST, no user input
            code,
            {"np": np},
            {
                "cfg": SimpleNamespace(seed=42),
                "node": SimpleNamespace(id=node_id, microstates=[None] * 8),
                "state": SimpleNamespace(id=state_id, trials_spent=0),
            },
        )

    nodes = {f"c{i:04d}": seed(f"c{i:04d}", f"c{i:04d}/m0000") for i in range(100)}
    states = {f"c0000/m{i:04d}": seed("c0000", f"c0000/m{i:04d}") for i in range(8)}
    assert len(set(nodes.values())) == 10 and len(set(states.values())) == 1
    return {
        "source_line": assignments[0].lineno,
        "source_expression": ast.unparse(expression),
        "100_node_unique_seeds": len(set(nodes.values())),
        "8_microstate_unique_seeds": len(set(states.values())),
        "example_nodes": {k: nodes[k] for k in ("c0000", "c0001", "c0010", "c0099")},
    }


def response_reward():
    registry = Registry(SearchConfig())
    first = Microstate("c0000/m0000", None, 0.0, "synthetic", {})
    second = Microstate("c0000/m0001", None, 0.0, "synthetic", {})
    node = ChemicalNode("c0000", "synthetic", {}, 0, "synthetic", [], [first, second])
    first.trials_spent = second.trials_spent = 1
    first.outcomes = Counter({"c0000": 1})
    second.outcomes = Counter({"c0000": 1})
    same_node = registry._response_term(node, first)
    first.outcomes = Counter({"failed": 1})
    unique_failure = registry._response_term(node, first)
    assert same_node == 0.0 and unique_failure == 1.0
    return {
        "same_node_or_different_self_loop_channels_reward": same_node,
        "otherwise_unseen_failure_reward": unique_failure,
        "limitation": "note_trial receives target chemical node or status, not event/channel identity",
    }


def external_field_shape():
    template = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    offsets = np.array([[0.0, 0.0, 0.0], [0.5, 1.0, 0.0], [1.0, 0.0, 0.0]])
    positions = template[None] + offsets[:, None]
    path = candidate("moving-in-field", None, positions)
    reported = evaluate_candidate(path, ParsConfig(refinement_remove_rigid_motion=False))[
        "path_curvature"
    ]
    raw = resample_positions(positions, 9).reshape(9, -1)
    length = np.linalg.norm(np.diff(positions.reshape(3, -1), axis=0), axis=1).sum()
    spacing = length / 8
    raw_curvature = float(np.square(np.diff(raw, n=2, axis=0)).sum() / spacing**3)
    assert reported == 0.0 and raw_curvature > 0.0
    return {
        "rigid_motion_removal": False,
        "reported_curvature": reported,
        "raw_curvature": raw_curvature,
    }


def recrossing_semantics():
    a, b, c = (frozenset({pair}) for pair in ((0, 1), (1, 2), (0, 2)))
    cases = {
        "A_B_A_B": recrossing_count((a, b, a, b)),
        "A_B_C_B_C_B": recrossing_count((a, b, c, b, c, b)),
    }
    assert cases == {"A_B_A_B": 1, "A_B_C_B_C_B": 0}
    return {
        "counts": cases,
        "meaning": "source-bond-set returns, not all state boundaries or kinetic basin recrossings",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    output = parser.parse_args().out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = [
        "src/prrs/search.py",
        "src/prrs/network.py",
        "pars/src/pars/frontier/objectives.py",
        "pars/src/pars/frontier/normalize.py",
        "pars/src/pars/frontier/margin.py",
        "pars/src/pars/paths/cluster.py",
        "pars/src/pars/paths/fingerprint.py",
        "pars/src/pars/paths/resample.py",
        "benchmarks/project_audit/algorithm_deep_dig.py",
    ]
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    report = {
        "kind": "algorithm-stage-synthetic-audit",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "numpy": np.__version__,
        "source_sha256": hashes,
        "limits": "No PES/MD, production ingestion, historical-run contamination, real saddles or scientific coverage measured.",
        "spectator_bond_dilution": spectator_dilution(),
        "leader_name_dependence": leader_order(),
        "normalization_pool_dependence": pool_normalization(),
        "direction_seed_collisions": source_seed_collisions(),
        "reservoir_response_reward": response_reward(),
        "external_field_shape": external_field_shape(),
        "recrossing_scope": recrossing_semantics(),
    }
    assert hashes == {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
