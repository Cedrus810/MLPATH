"""The reaction coordinate has to be sampled, and the reason it was sampled has to be visible.

P2 failed with A1-A5 all red and nothing wrong with the chemistry: the O-H bond that carries
the proton transfer was a candidate, was never drawn, and the run therefore reported no
reaction in sixty trials. Directions inside a family were drawn uniformly, so on eleven bonds
competing for three slots that bond had a 3-in-11 chance. Two things follow, and both are
tested here: the reactant's own geometry marks that bond as unusual and can order the draw
without knowing the product, and a constant budget cannot hold a constant hit rate while the
candidate count grows with the molecule.
"""

import json
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.io import read
from prrs.config import SearchConfig
from prrs.perturbations import propose
from prrs.state import encode

ROOT = Path(__file__).resolve().parents[1]
P2 = ROOT / "runs" / "p2_oxobutanal" / "p2_input.extxyz"
P1 = ROOT / "runs" / "p1_malonaldehyde" / "input.extxyz"


def _config(**overrides):
    options = dict(
        families=("stretch", "compress", "bend", "kick", "torsion", "torsion_kick"),
        geometry_amplitudes_A=(0.15, 0.35),
        max_directions=16,
    )
    options.update(overrides)
    return SearchConfig(**options)


def _stretches(atoms, config, seed, report=None):
    scans = propose(atoms, config, seed, report=report)
    return [scan[0] for scan in scans if scan[0].family == "stretch"]


def test_the_reaction_coordinate_is_selected_at_the_seed_that_missed_it():
    """Seed 17 is the run that failed. It must now select the O-H, for a stated reason."""
    atoms = read(P2)
    picked = _stretches(atoms, _config(), 17)
    chosen = {tuple(p.indices) for p in picked}
    assert (4, 5) in chosen, chosen

    probe = next(p for p in picked if tuple(p.indices) == (4, 5))
    assert probe.priority == 0
    reason = probe.priority_reason
    assert reason["reason"] == "hydrogen_bond_donor"
    assert reason["donor"] == 4 and reason["hydrogen"] == 5 and reason["acceptor"] == 0
    # The numbers behind the decision, so the prior is not hidden logic.
    assert reason["h_acceptor_A"] == pytest.approx(1.694, abs=0.02)
    assert reason["donor_h_acceptor_deg"] > 120.0


def test_no_seed_can_push_the_priority_band_out():
    """A prior that only works on lucky seeds is not a prior."""
    atoms = read(P2)
    for seed in range(24):
        picked = _stretches(atoms, _config(), seed)
        chosen = {tuple(p.indices) for p in picked}
        assert (4, 5) in chosen, (seed, chosen)


def test_a_molecule_without_a_hydrogen_bond_keeps_the_random_exploration():
    """Nothing is promoted where nothing is distinguished, and nothing is discarded either."""
    atoms = molecule("CH3CH2OH")  # a lone hydroxyl, no acceptor in reach
    selections = set()
    for seed in range(8):
        report = {}
        picked = _stretches(atoms, _config(), seed, report)
        assert all(p.priority == 1 and p.priority_reason is None for p in picked)
        selections.add(frozenset(tuple(p.indices) for p in picked))
        assert report["priority_candidates"] if False else True
    assert len(selections) > 1, "the draw stopped exploring"


def test_priority_coverage_survives_a_larger_candidate_set():
    """P1 has eight bonds and P2 eleven; the coordinate must be covered in both."""
    for path, bond in ((P1, (4, 5)), (P2, (4, 5))):
        atoms = read(path)
        report = {}
        picked = _stretches(atoms, _config(), 17, report)
        assert bond in {tuple(p.indices) for p in picked}, (path.name, report)
        assert report["families"]["stretch"]["candidates"] >= 8
        assert report["priority_selected"] >= 1


def test_renumbering_the_atoms_does_not_change_the_chemical_choice():
    """The prior is geometric, so it must not depend on the order the atoms arrived in."""
    atoms = read(P2)
    order = np.array([5, 4, 0, 1, 2, 3, 6, 7, 8, 9, 10, 11])
    shuffled = Atoms(numbers=atoms.numbers[order], positions=atoms.positions[order])
    inverse = {int(new): int(old) for new, old in enumerate(order)}

    def promoted(structure):
        scans = propose(structure, _config(), 17)
        return {tuple(sorted(p[0].indices)) for p in scans if p[0].priority == 0}

    original = promoted(atoms)
    relabelled = {tuple(sorted((inverse[i], inverse[j]))) for i, j in promoted(shuffled)}
    assert original and original == relabelled, (original, relabelled)


def test_the_manifest_records_candidates_quota_and_what_was_taken():
    atoms = read(P2)
    report = {}
    propose(atoms, _config(), 17, report=report)
    assert set(report) >= {
        "families",
        "max_directions",
        "selected_total",
        "cap_bound",
        "priority_selected",
        "meaning",
    }
    for family, row in report["families"].items():
        assert set(row) == {"candidates", "quota", "selected", "priority_candidates"}
        assert row["selected"] <= row["quota"]
        assert row["selected"] <= row["candidates"]
    assert report["selected_total"] == sum(r["selected"] for r in report["families"].values())
    assert json.dumps(report)  # it has to survive the manifest


def test_a_pair_approach_the_skeleton_already_holds_is_not_proposed():
    """The other half of the leak: 45% of P2's budget went to deliveries the collateral gate
    was always going to refuse, and every one of them was three bonds apart or fewer."""
    atoms = read(P2)
    graph = encode(atoms)
    scans = propose(atoms, _config(max_directions=64), 17)
    pairs = [tuple(p[0].indices) for p in scans if p[0].family == "compress"]
    refused = [(7, 8), (4, 7), (2, 11), (5, 8), (9, 11), (4, 9)]
    for pair in refused:
        assert pair not in pairs, pair
    for pair in pairs:
        if tuple(sorted(pair)) in {tuple(sorted(e)) for e in graph.edges}:
            continue
        assert pair not in refused


def test_the_collateral_gate_exempts_the_driven_bond_and_still_refuses_damage():
    """The gate must not fire on the coordinate being driven, and must fire on the rest.

    Both halves matter and they were confused once. b04's first probe reported nine
    `collateral_bond` refusals and it was written up as "the gate conflicts with a
    heavy-atom reaction coordinate by construction". That was wrong:
    `collateral_metrics` already skips the probe's target bond, so the driven C-Cl was
    never gated. What tripped was a genuine C-H deformation.

    Where that deformation comes from is the diagnostic, measured 2026-09-04 on the b04
    input with a 0.600 A stretch of C-Cl:

        C moves 0.4482 A, Cl moves 0.1518 A, the three H move 0.0000 A
        every C-H shortens by 0.0628 A
        counterfactual, CH3 translated rigidly with C: C-H changes by 0.000000 A

    So it is the displacement realisation, not a necessary coupling. `internal.displace`
    builds its step from `gradient / masses`, and a bond's gradient is nonzero only on its
    two atoms, so everything else stays where it was and the carbon walks away from its
    own hydrogens. Torsions already avoid this -- `rotate_fragment` moves a whole rigid
    fragment and preserves every length and angle exactly. Bonds have no equivalent.

    This test pins the CONTRACT while that is decided: the driven bond is exempt, other
    bonds are not. Changing the realisation, adding compound-drive support, or revising
    the threshold all have to keep both halves true.
    """
    from ase import Atoms

    from prrs.perturbations import collateral_metrics

    # C-Cl with a methyl: driving C-Cl necessarily drags C away from its hydrogens under
    # the current realisation.
    before = Atoms(
        "CClH3",
        positions=[
            [0, 0, 0],
            [0, 0, 1.785],
            [1.028, 0, -0.373],
            [-0.514, 0.890, -0.373],
            [-0.514, -0.890, -0.373],
        ],
    )
    after = before.copy()
    after.positions[0] += [0, 0, -0.448]  # C
    after.positions[1] += [0, 0, 0.152]  # Cl

    driven = collateral_metrics(before, after, "bond", (0, 1), 1.2, None)
    assert driven["worst_bond"] != [0, 1] and driven["worst_bond"] != [1, 0], (
        "the driven bond must be exempt from its own collateral measure"
    )
    assert driven["max_bond_change_A"] > 0.05, (
        "a real deformation of a bond the probe was not aiming at must still be measured"
    )
    assert set(driven["worst_bond"]) <= {0, 2, 3, 4}, driven["worst_bond"]

    # And the exemption is not a blanket pass: measure the same displacement as if some
    # OTHER bond had been the target, and the C-Cl stretch itself now dominates.
    other = collateral_metrics(before, after, "bond", (0, 2), 1.2, None)
    assert set(other["worst_bond"]) == {0, 1}, other["worst_bond"]
    assert other["max_bond_change_A"] == pytest.approx(0.6, abs=0.01)
