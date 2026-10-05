"""Connect mode, G3 / intervention I-b (PLAN_2026-09-28.md): the key difference goes first.

With `connect_direction_priority` on and a product given, the directions that act on the
bonds the product needs changed -- stretch a bond in B-, compress a pair in B+ -- take a
band of their own ahead of the hydrogen-bond priority band. Off (the default) the proposal
is the one that existed before, byte for byte. The draft put these directions in band -1
while the allocation loop only ever visited bands 0 and 1, so the directions it meant to
put first were never selected at all; the first test here is the one that catches that.
"""

from dataclasses import replace

import pytest
from ase.optimize import LBFGS

from prrs.calculators import ProtonRelay, proton_relay_atoms, proton_relay_factory
from prrs.config import SearchConfig
from prrs.network import Registry
from prrs.perturbations import propose
from prrs.product import product_reference
from prrs.search import ReactionSearch

CONFIG = SearchConfig(families=("stretch", "compress"), closed_shell_only=False)
# relay numbering: heavy 0..3, proton k is atom 4 + k; R has proton k on atom k
BROKEN = {frozenset((k, 4 + k)) for k in range(3)}
FORMED = {frozenset((k + 1, 4 + k)) for k in range(3)}


def _minimum(sides):
    atoms = proton_relay_atoms(sides)
    atoms.calc = ProtonRelay()
    LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=2000)
    atoms.calc = None
    return atoms


def _reference():
    return product_reference(_minimum((0, 0, 0)), _minimum((1, 1, 1)), Registry(CONFIG))


def _flat(proposal):
    return [(p[0].family, frozenset(p[0].indices), p[0].priority) for p in proposal]


def _wanted():
    return {("stretch", pair) for pair in BROKEN} | {("compress", pair) for pair in FORMED}


def test_the_key_difference_directions_are_selected_first():
    reactant = _minimum((0, 0, 0))
    difference = _reference().key_difference(reactant)
    assert difference == (BROKEN, FORMED)
    report = {}
    proposal = _flat(propose(reactant, CONFIG, 7, report=report, key_difference=difference))
    head = [entry for entry in proposal if entry[2] == -1]
    assert proposal[: len(head)] == head, "band -1 must be allocated before any other band"
    assert {(family, pair) for family, pair, _ in head} <= _wanted()
    # I-b reorders inside each family's quota and adds no budget: stretch has three
    # candidates and a quota of two, so one B- stretch is left out -- and named as such
    expected = sum(
        min(report["families"][family]["quota"], sum(1 for f, _ in _wanted() if f == family))
        for family in ("stretch", "compress")
    )
    assert len(head) == expected == 5
    left_out = {(f, p) for f, p in _wanted()} - {(f, p) for f, p, _ in head}
    for family, pair in left_out:
        a, b = sorted(pair)
        sign = "+" if family == "stretch" else "-"
        assert f"{family}:{a}-{b}:{sign}" in report["uncovered_ids"]


def test_exhaustive_mode_puts_every_key_difference_direction_first():
    reactant = _minimum((0, 0, 0))
    config = replace(CONFIG, direction_mode="exhaustive")
    difference = _reference().key_difference(reactant)
    proposal = _flat(propose(reactant, config, 7, key_difference=difference))
    head = proposal[: len(_wanted())]
    assert {(family, pair) for family, pair, _ in head} == _wanted()
    assert all(priority == -1 for _, _, priority in head)
    assert all(priority != -1 for _, _, priority in proposal[len(_wanted()) :])


def test_off_the_proposal_is_exactly_the_old_one():
    reactant = _minimum((0, 0, 0))
    old = propose(reactant, CONFIG, 7)
    same = propose(reactant, CONFIG, 7, key_difference=None)
    assert _flat(old) == _flat(same)
    assert all(priority != -1 for _, _, priority in _flat(old))


def test_the_report_names_what_was_put_first_and_why():
    reactant = _minimum((0, 0, 0))
    report = {}
    propose(
        reactant, CONFIG, 7, report=report, key_difference=_reference().key_difference(reactant)
    )
    assert report["connect_priority_selected"] == 5
    assert set(map(tuple, report["connect_key_difference"]["broken"])) == {
        tuple(sorted(p)) for p in BROKEN
    }
    plain = {}
    propose(reactant, CONFIG, 7, report=plain)
    assert "connect_priority_selected" not in plain
    assert "connect_key_difference" not in plain


def test_the_difference_is_taken_from_the_current_state_not_the_reactant():
    middle = _minimum((1, 0, 1))
    broken, formed = _reference().key_difference(middle)
    assert broken == {frozenset((1, 5))}
    assert formed == {frozenset((2, 5))}


def test_the_switch_without_a_product_refuses_to_start(tmp_path):
    config = replace(CONFIG, connect_direction_priority=True)
    with pytest.raises(ValueError, match="product"):
        ReactionSearch(proton_relay_factory, config).run(
            proton_relay_atoms((0, 0, 0)), tmp_path / "run"
        )
    assert not (tmp_path / "run").exists()


def test_a_connect_run_records_the_priority_on_its_scans(tmp_path):
    config = replace(CONFIG, connect_direction_priority=True, max_trials=6, max_depth=2)
    network = ReactionSearch(proton_relay_factory, config).run(
        proton_relay_atoms((0, 0, 0)), tmp_path / "run", product=_minimum((1, 1, 1))
    )
    first = network["amplitude_scans"][0]
    assert first["proposal"]["priority"] == -1
    assert first["proposal"]["priority_reason"]["connect"] in ("B-", "B+")


def test_tied_numberings_of_the_product_all_contribute_to_the_difference():
    """Methanol -> formaldehyde + H2: any of the three methyl hydrogens can be the one
    that leaves, all three at the same distance, so all three C-H bonds are in B- and
    all three H-H pairs in B+. Taking one numbering would put one hydrogen first on an
    arbitrary choice made when the product file was written."""
    from ase.build import molecule

    methanol = molecule("CH3OH")  # C0 O1 H2 H3(hydroxyl) H4 H5
    product = methanol.copy()
    product.positions[2] = [5.0, 0.0, 0.0]
    product.positions[3] = [5.74, 0.0, 0.0]
    reference = product_reference(methanol, product, Registry(SearchConfig()))
    broken, formed = reference.key_difference(methanol)
    assert broken == {frozenset(p) for p in ((0, 2), (0, 4), (0, 5), (1, 3))}
    assert formed == {frozenset(p) for p in ((2, 3), (3, 4), (3, 5))}
