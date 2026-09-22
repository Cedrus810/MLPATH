"""PLAN item 6.1: can a Registry be rebuilt from what a run wrote down?

Resume (item 6) is gated on this test existing and passing -- if reconstruction
cannot reproduce the in-memory registry field by field, resume is abandoned rather
than done wrong (HANDOFF §9: acceptance order affects results, and a wrong resume is
worse than none). Compared here, field by field, against the live registry:

node id / key / depth / discovered_by / permutations / trials_spent / evicted /
conformer_transitions / components, and per microstate id / energy / discovered_by /
trials_spent / outcomes / score / free_bonds / torsions.

The one field not serialised is `permutations`; from_network re-enumerates it from the
node's first structure, and the test holds that re-enumeration to be the SAME group in
the SAME order -- deterministic enumeration is what makes that legitimate.
"""

import numpy as np
import pytest
from ase.build import molecule
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

from prrs.config import SearchConfig
from prrs.network import Registry
from prrs.perturbations import step_torsion, torsion_on_bond
from prrs.state import encode
from prrs.chemistry import canonical_labels


def _with_quench(atoms, quench, energy):
    atoms = atoms.copy()
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=np.zeros((len(atoms), 3)))
    atoms.info["quench"] = {"rounds": [{"modes": quench}]}
    return atoms


@pytest.fixture
def scenario(tmp_path):
    """One node, one eviction, some spent trials: everything resume must reproduce."""
    source = molecule("CH3CH2OH")
    graph = encode(source)
    labels = canonical_labels(source.numbers, graph.edges, graph.index)
    torsion = torsion_on_bond(graph, labels, (1, 2))  # the O-H torsion: anti vs gauche
    gauche, _ = step_torsion(source, torsion, np.deg2rad(120.0), in_place=False)
    rotamer, _ = step_torsion(source, torsion, np.deg2rad(240.0), in_place=False)
    hydroxyl = tuple(torsion.indices)
    free = [{"indices": list(hydroxyl), "tier": "free", "coordinate": "torsion"}]
    soft = [{"indices": list(hydroxyl), "tier": "flat_biased", "coordinate": "torsion"}]

    config = SearchConfig(conformer_max_per_node=2)
    live = Registry(config)
    # free on one side, measured-not-free on the other: the coordinate is NOT
    # quotiented (test_curvature's rule), so anti and gauche stay two microstates
    anti = _with_quench(source, free, -100.000)
    first = live.admit(anti, None, "t000000")
    assert first.outcome == "new_chemical_node"
    second = live.admit(_with_quench(gauche, soft, -100.002), first.node, "t000001")
    assert second.outcome == "new_microstate"
    # spend both residents' conformer budgets, then the rotamer evicts the weaker one
    for resident in list(first.node.microstates):
        for _ in range(config.conformer_trials_per_microstate):
            live.note_trial(first.node, resident, "completed")
    third = live.admit(_with_quench(rotamer, soft, -100.050), first.node, "t000002")
    assert third.outcome == "new_microstate"
    assert first.node.evicted and first.node.evicted[0]["id"] == "c0000/m0000"
    live.note_trial(third.node, third.microstate, "completed")

    # serialize exactly the way publish() serialises chemical_nodes
    (tmp_path / "chemical").mkdir()
    network = {"chemical_nodes": []}
    for node in live.nodes:
        microstates = []
        for state in node.microstates:
            relative = f"chemical/{state.id.replace('/', '_')}.extxyz"
            frame = state.structure.copy()
            frame.calc = None
            write(str(tmp_path / relative), frame)
            microstates.append(state.summary(relative))
        network["chemical_nodes"].append(
            {
                "id": node.id,
                "chemical_key": node.key,
                "depth": node.depth,
                "discovered_by": node.discovered_by,
                "key_components": {k: v for k, v in node.components.items() if k != "key"},
                "trials_spent": node.trials_spent,
                "reservoir": {
                    "capacity": config.conformer_max_per_node,
                    "occupied": len(node.microstates),
                    "evicted": node.evicted,
                },
                "microstates": microstates,
                "conformer_transitions": node.conformer_transitions,
            }
        )
    return live, network, config, tmp_path, anti


def test_from_network_reproduces_the_live_registry_field_by_field(scenario):
    live, network, config, output, _ = scenario
    rebuilt = Registry.from_network(network, output, config)

    assert len(rebuilt.nodes) == len(live.nodes)
    for node, original in zip(rebuilt.nodes, live.nodes):
        assert node.id == original.id
        assert node.key == original.key
        assert node.depth == original.depth
        assert node.discovered_by == original.discovered_by
        assert node.trials_spent == original.trials_spent
        assert node.evicted == original.evicted
        assert node.conformer_transitions == original.conformer_transitions
        assert node.components.keys() == original.components.keys()
        for field in original.components:
            if field == "automorphisms":
                continue  # the info dict carries counts; re-enumeration refreshes it
            assert node.components[field] == original.components[field]
        # permutations are not serialised: the re-enumerated group must be the same
        # arrays in the same order
        assert len(node.permutations) == len(original.permutations)
        assert all(
            np.array_equal(a, b) for a, b in zip(node.permutations, original.permutations)
        )
        for state, original_state in zip(node.microstates, original.microstates):
            assert state.id == original_state.id
            assert state.energy_eV == pytest.approx(original_state.energy_eV)
            assert state.discovered_by == original_state.discovered_by
            assert state.trials_spent == original_state.trials_spent
            assert dict(state.outcomes) == dict(original_state.outcomes)
            assert state.score == pytest.approx(original_state.score)
            assert tuple(state.free_bonds) == tuple(original_state.free_bonds)
            assert set(state.torsions) == set(original_state.torsions)
            for key in state.torsions:
                assert state.torsions[key] == pytest.approx(original_state.torsions[key])


def test_rebuilt_registry_behaves_like_the_live_one(scenario):
    """Reconstruction is not a museum piece: it must keep admitting correctly."""
    live, network, config, output, anti = scenario
    rebuilt = Registry.from_network(network, output, config)
    live_node = live.nodes[0]

    # the anti conformer is the EVICTED one; the resident is gauche. Re-admitting the
    # resident's own geometry must match, through the re-enumerated group.
    resident = live_node.microstates[0]
    again = _with_quench(resident.structure, [], resident.energy_eV)
    admission = rebuilt.admit(again, rebuilt.nodes[0], "t000009")
    assert admission.outcome == "existing_microstate"

    # a new conformer takes the next id in the numbering: two residents plus one
    # evicted before it, so the next insert is m0003 -- not m0002, which would collide
    # with history and silently renumber the reservoir
    source = molecule("CH3CH2OH")
    graph = encode(source)
    labels = canonical_labels(source.numbers, graph.edges, graph.index)
    torsion = torsion_on_bond(graph, labels, (1, 2))
    far, _ = step_torsion(source, torsion, np.deg2rad(60.0), in_place=False)
    free = [{"indices": list(torsion.indices), "tier": "free", "coordinate": "torsion"}]
    # spend the rebuilt residents' budgets so one is evictable, then admit
    for resident in list(rebuilt.nodes[0].microstates):
        for _ in range(config.conformer_trials_per_microstate):
            rebuilt.note_trial(rebuilt.nodes[0], resident, "completed")
    admission = rebuilt.admit(_with_quench(far, free, -100.080), rebuilt.nodes[0], "t000010")
    assert admission.outcome == "new_microstate"
    assert admission.microstate.id == "c0000/m0003"


def test_missing_structure_file_is_refused(scenario):
    _, network, config, output, _ = scenario
    (output / "chemical/c0000_m0002.extxyz").unlink()
    with pytest.raises(ValueError, match="missing"):
        Registry.from_network(network, output, config)
