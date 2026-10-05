"""Connect mode, G1 (PLAN_2026-09-28.md): the product is recorded, never routed.

Methanol -> formaldehyde + H2 carries the one property d has to get right: the three
methyl hydrogens are indistinguishable in the reactant, so which of them the product
file names in the H-H bond is an arbitrary choice, and a state in which a different one
did the chemistry must be at distance 0. The double well carries the search wiring.
"""

import json
from dataclasses import replace

import numpy as np
import pytest
from ase.build import molecule

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.network import Registry
from prrs.product import ProductMismatch, product_reference
from prrs.schema import validate
from prrs.search import ReactionSearch

# molecule("CH3OH") numbering: C0 O1 H2 (methyl) H3 (hydroxyl) H4 H5 (methyl)
METHYL = (2, 4, 5)
HYDROXYL = 3


def _methanol():
    return molecule("CH3OH")


def _formaldehyde_plus_h2(leaving=2):
    """Methyl hydrogen `leaving` and the hydroxyl hydrogen leave together as H2."""
    atoms = _methanol()
    positions = atoms.positions.copy()
    positions[leaving] = [5.0, 0.0, 0.0]
    positions[HYDROXYL] = [5.74, 0.0, 0.0]
    atoms.set_positions(positions)
    return atoms


def _swapped(atoms, a, b):
    """The same geometry with the labels of atoms a and b exchanged."""
    out = atoms.copy()
    positions = out.positions.copy()
    positions[[a, b]] = positions[[b, a]]
    out.set_positions(positions)
    return out


def _reference(product, config=None):
    registry = Registry(config or SearchConfig())
    return product_reference(_methanol(), product, registry)


def test_the_product_must_share_the_reactants_atom_numbering():
    reordered = _formaldehyde_plus_h2()
    reordered.numbers[[0, 1]] = reordered.numbers[[1, 0]]
    with pytest.raises(ProductMismatch, match="element sequence"):
        _reference(reordered)
    with pytest.raises(ProductMismatch, match="atoms"):
        _reference(_formaldehyde_plus_h2()[:-1])


def test_distance_is_zero_at_the_product_and_counts_bond_changes_elsewhere():
    reference = _reference(_formaldehyde_plus_h2())
    assert reference.distance(_formaldehyde_plus_h2()) == 0
    # C-H2 and O-H3 break, H2-H3 forms
    assert reference.distance(_methanol()) == 3
    assert reference.group_order == 6  # the three methyl hydrogens


def test_an_equivalent_hydrogen_doing_the_chemistry_is_still_the_product():
    """The minimum over Aut(G_R) is what makes this 0; without it the answer is 4."""
    reference = _reference(_formaldehyde_plus_h2(leaving=2))
    other = _formaldehyde_plus_h2(leaving=4)
    assert reference.distance(other) == 0
    assert Registry(SearchConfig()).key_of(other)["key"] == reference.key


def test_renumbering_the_product_within_the_reactants_symmetry_changes_no_distance():
    first = _reference(_formaldehyde_plus_h2(leaving=2))
    second = _reference(_swapped(_formaldehyde_plus_h2(leaving=2), 2, 5))
    for structure in (_methanol(), _formaldehyde_plus_h2(2), _formaldehyde_plus_h2(4)):
        assert first.distance(structure) == second.distance(structure)
    assert first.key == second.key


def test_a_reactant_group_that_does_not_certify_is_refused_not_replaced():
    config = replace(SearchConfig(), automorphism_limit=1)
    with pytest.raises(ProductMismatch, match="automorphism_limit"):
        _reference(_formaldehyde_plus_h2(), config)


# ---- search wiring (double well: C2 at 1.2 A bonded, at 2.4 A not) -----------------


def _config():
    return SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
        timestep_fs=0.1,
        response_steps=20,
        quench_fmax_eV_A=0.001,
        closed_shell_only=False,
        max_depth=2,
        max_trials=8,
        refinement_steps=2,
    )


def _dissociated():
    atoms = demo_atoms()
    atoms.positions[:, 0] *= 2.0
    return atoms


def test_a_run_given_a_product_records_distances_and_compatibility(tmp_path):
    network = ReactionSearch(double_well_factory, _config()).run(
        demo_atoms(), tmp_path / "run", product=_dissociated()
    )
    assert network["product"]["d_source"] == 1
    assert network["product"]["structure"] == "product.extxyz"
    assert (tmp_path / "run" / "product.extxyz").exists()
    manifest = json.loads((tmp_path / "run" / "config.json").read_text())
    assert len(manifest["product_structure_sha256"]) == 64
    compatible = {n["id"]: n["product_compatible"] for n in network["chemical_nodes"]}
    assert sum(compatible.values()) == 1, "the dissociated well must have been found"
    for node in network["chemical_nodes"]:
        for state in node["microstates"]:
            assert state["d_to_product"] == (0 if node["product_compatible"] else 1)
    assert validate(network) == []


def test_a_run_without_a_product_writes_nothing_new(tmp_path):
    network = ReactionSearch(double_well_factory, _config()).run(demo_atoms(), tmp_path / "run")
    assert "product" not in network
    for node in network["chemical_nodes"]:
        assert "product_compatible" not in node
        assert all("d_to_product" not in state for state in node["microstates"])
    manifest = json.loads((tmp_path / "run" / "config.json").read_text())
    assert "product_structure_sha256" not in manifest
    assert not (tmp_path / "run" / "product.extxyz").exists()


def test_a_mismatched_product_is_refused_before_anything_is_written(tmp_path):
    wrong = _dissociated()
    wrong.numbers[1] = 1
    with pytest.raises(ProductMismatch):
        ReactionSearch(double_well_factory, _config()).run(
            demo_atoms(), tmp_path / "run", product=wrong
        )
    assert not (tmp_path / "run").exists()


def test_the_schema_holds_a_product_run_to_its_fields(tmp_path):
    network = ReactionSearch(double_well_factory, _config()).run(
        demo_atoms(), tmp_path / "run", product=_dissociated()
    )
    del network["chemical_nodes"][0]["microstates"][0]["d_to_product"]
    assert any("d_to_product" in line for line in validate(network))


def _running(tmp_path, product):
    """A run left at status "running", as a kill leaves it."""
    output = tmp_path / "run"
    ReactionSearch(double_well_factory, _config()).run(demo_atoms(), output, product=product)
    network = json.loads((output / "network.json").read_text())
    network["status"] = "running"
    (output / "network.json").write_text(json.dumps(network))
    return output


@pytest.mark.parametrize(
    "given, message",
    [
        (None, "product"),
        ("moved", "product"),
    ],
)
def test_resume_refuses_a_different_product(tmp_path, given, message):
    output = _running(tmp_path, _dissociated())
    product = None
    if given == "moved":
        product = _dissociated()
        product.positions[1, 0] += 0.1
    with pytest.raises(ValueError, match=message):
        ReactionSearch(double_well_factory, _config()).run(
            demo_atoms(), output, resume=True, product=product
        )


def test_resume_refuses_a_product_the_original_run_did_not_have(tmp_path):
    output = tmp_path / "run"
    ReactionSearch(double_well_factory, _config()).run(demo_atoms(), output)
    network = json.loads((output / "network.json").read_text())
    network["status"] = "running"
    (output / "network.json").write_text(json.dumps(network))
    with pytest.raises(ValueError, match="product"):
        ReactionSearch(double_well_factory, _config()).run(
            demo_atoms(), output, resume=True, product=_dissociated()
        )


def test_distances_are_integers_and_never_negative():
    reference = _reference(_formaldehyde_plus_h2())
    rng = np.random.default_rng(0)
    for _ in range(20):
        jiggled = _methanol()
        jiggled.positions += 0.05 * rng.standard_normal(jiggled.positions.shape)
        d = reference.distance(jiggled)
        assert isinstance(d, int) and d >= 0
