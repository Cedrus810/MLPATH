"""The delivered probe, continued past the frame where the topology switched.

The saddle search and the two-sided descent establish that two minima are joined. They do
not establish that the motion which was delivered arrives at the product, and on
3-oxobutanal round three those two facts came apart: three saddles supported the proton
transfer channel and no trial had observed A reaching B. These tests pin the three things
that decide whether the continuation is honest -- which way it goes, how far each rung
moves, and what a rung's endpoint is allowed to be called.
"""

from dataclasses import replace
import json
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from prrs.calculators import double_well_factory
from prrs.config import SearchConfig
from prrs.runner import PathRecorder, continue_across_boundary


@pytest.fixture
def config():
    return SearchConfig(
        closed_shell_only=False,
        quench_fmax_eV_A=0.001,
        boundary_continuation_step_A=0.05,
        boundary_continuation_max_steps=8,
        tracked_coordinates=({"name": "r", "kind": "bond", "indices": (0, 1)},),
    )


def ridge_frame():
    """A frame on the double well's ridge, just short of the far basin."""
    return Atoms("C2", positions=[[-0.85, 0, 0], [0.85, 0, 0]])  # r = 1.70


def test_the_walk_goes_the_way_the_tangent_points(config):
    """Sign is carried, not discarded. The same frame with opposite tangents must not
    produce the same verdict, or the continuation would be measuring nothing."""
    frame = ridge_frame()
    outward = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    seen = {}

    def name(endpoint, confirmed, checked):
        r = float(endpoint.get_distance(0, 1))
        far = abs(r - 2.4) < abs(r - 1.2)
        return {
            "verdict": "reached_other_state" if far else "returned_to_source",
            "stop": far,
            "r": r,
        }

    for label, tangent in (("forward", outward), ("backward", -outward)):
        _, report = continue_across_boundary(
            frame, tangent, double_well_factory, config, verdict=name
        )
        seen[label] = report
    assert seen["forward"]["outcome"] == "reached_other_state"
    assert all(rung["verdict"] == "returned_to_source" for rung in seen["backward"]["rungs"])
    assert seen["backward"]["outcome"] == "exhausted"


def test_every_rung_respects_the_single_atom_step_cap(tmp_path, config):
    """The cap is on the largest single-atom displacement, so a rung's geometry is
    predictable from the rung number alone; a cap on the norm would let both atoms move
    the full step and the walk would advance twice as fast as it reports."""
    frame = ridge_frame()
    tangent = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    moved = []

    def name(endpoint, confirmed, checked):
        return {"verdict": "returned_to_source"}

    _, report = continue_across_boundary(
        frame, -tangent, double_well_factory, config, verdict=name
    )
    for rung in report["rungs"]:
        moved.append(rung["displacement_A"])
    assert moved == pytest.approx([0.05 * n for n in range(1, len(moved) + 1)])
    # And the recorded geometries really are that far from the seed.
    recorder = PathRecorder(tmp_path / "path", config, "continuation")
    continue_across_boundary(
        frame, -tangent, double_well_factory, config, verdict=name, recorder=recorder
    )
    rows = [
        json.loads(line)
        for line in (recorder.directory / "continuation.jsonl").read_text().splitlines()
    ]
    rows = [row for row in rows if row["phase"] == "continuation"]
    for index, row in enumerate(rows, start=1):
        assert row["tracked"]["r"] == pytest.approx(1.70 - 2 * 0.05 * index, abs=1e-9)
        assert row["potential_eV"] is not None


def test_a_rung_that_leaves_the_domain_is_named_as_such(config):
    """Driving two atoms into each other is a step-size failure, not a chemical result,
    and the report has to say which of the two it was."""
    tight = replace(config, min_distance_ratio=0.9, boundary_continuation_max_steps=40)
    frame = ridge_frame()
    tangent = np.array([[1.0, 0, 0], [-1.0, 0, 0]])  # compress the pair

    def name(endpoint, confirmed, checked):
        return {"verdict": "returned_to_source"}

    _, report = continue_across_boundary(
        frame, tangent, double_well_factory, tight, verdict=name
    )
    assert report["outcome"] == "left_domain"
    assert report["rungs"][-1]["verdict"] == "left_domain"
    assert report["rungs"][-1]["gate"]


def test_a_tangent_with_no_internal_component_is_refused(config):
    """A pure translation moves no coordinate. Continuing along it would quench straight
    back to the seed's own basin at every rung and report that as a finding."""
    frame = ridge_frame()
    translation = np.array([[1.0, 0, 0], [1.0, 0, 0]])
    _, report = continue_across_boundary(frame, translation, double_well_factory, config)
    assert report["outcome"] == "tangent_has_no_internal_component"
    assert report["rungs"] == []


def test_stopping_on_a_saddle_is_not_reported_as_a_basin(config):
    """A rung that quenches onto the ridge means the delivered path runs along the
    dividing surface. It must not be handed to the caller for naming as a state."""
    seen = []

    def name(endpoint, confirmed, checked):
        seen.append(endpoint)
        return {"verdict": "returned_to_source"}

    # Seeded exactly at the barrier top, where the quench has nowhere to go.
    frame = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])
    tangent = np.array([[0.0, 1.0, 0], [0.0, -1.0, 0]])  # perpendicular: bond unchanged
    _, report = continue_across_boundary(
        frame, tangent, double_well_factory, config, verdict=name
    )
    verdicts = {rung["verdict"] for rung in report["rungs"]}
    assert verdicts <= {
        "stopped_on_saddle",
        "not_a_minimum",
        "returned_to_source",
        "quench_failed",
    }
    if "stopped_on_saddle" in verdicts:
        # The caller was not asked to name it, so nothing could be admitted from it.
        assert len(seen) < len(report["rungs"])


class SharedProton(Calculator):
    """A proton shared between two identical heavy atoms; a degenerate transfer.

    The two minima are mirror images, so they have the same chemical key while their bond
    sets differ -- exactly the case a key-only test gets wrong. Built as a double well in
    q = r(0,1) - r(1,2) with the heavy separation and the sum of the two distances held, so
    the proton is confined to the axis and the transfer is the only soft coordinate.
    """

    implemented_properties = ["energy", "forces"]
    span = 2.6
    width = 0.4

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.positions
        legs = {(0, 1): None, (1, 2): None, (0, 2): None}
        for pair in legs:
            vector = positions[pair[1]] - positions[pair[0]]
            legs[pair] = (float(np.linalg.norm(vector)), vector)
        r01, r12, r02 = (legs[(0, 1)][0], legs[(1, 2)][0], legs[(0, 2)][0])
        total, difference = r01 + r12, r01 - r12
        shape = (difference / self.width) ** 2 - 1
        energy = (
            30.0 * (r02 - self.span) ** 2 + 20.0 * (total - self.span) ** 2 + 0.5 * shape**2
        )
        d_total = 40.0 * (total - self.span)
        d_difference = 0.5 * 2 * shape * 2 * difference / self.width**2
        derivatives = {
            (0, 1): d_total + d_difference,
            (1, 2): d_total - d_difference,
            (0, 2): 60.0 * (r02 - self.span),
        }
        forces = np.zeros_like(positions)
        for (i, j), derivative in derivatives.items():
            distance, vector = legs[(i, j)]
            gradient = derivative * vector / distance
            forces[i] -= gradient
            forces[j] += gradient
        self.results = {"energy": float(energy), "forces": -forces}


def shared_proton_atoms(offset):
    """The proton at q = offset along the axis; +-0.4 are the two minima."""
    r01 = (SharedProton.span + offset) / 2
    return Atoms("CHC", positions=[[0.0, 0, 0], [r01, 0, 0], [SharedProton.span, 0, 0]])


def test_a_degenerate_transfer_is_not_a_return_to_the_source(config):
    """The defect this test exists for: the endpoint's chemical key is unchanged because
    the two wells are mirror images, so a key-only test reports "returned to the source"
    about a proton that has demonstrably moved to the other heavy atom."""
    from prrs.search import rung_verdict
    from prrs.chemistry import chemical_key

    left = shared_proton_atoms(-SharedProton.width)
    right = shared_proton_atoms(+SharedProton.width)
    key_left = chemical_key(left, config.bond_scale)
    key_right = chemical_key(right, config.bond_scale)
    assert key_left["key"] == key_right["key"], (
        "the two wells must be symmetry copies for this test"
    )
    answer = rung_verdict(left, right, key_left, key_right, config)
    assert answer["chemical_key_matches_source"] is True  # what misled the old test
    assert answer["broken"] and answer["formed"]  # and what settles it
    assert answer["transition_class"] == "degenerate_reaction"
    assert answer["verdict"] == "reached_degenerate_product"
    assert answer["new_state"] is False  # so it opens no node


def test_the_continuation_carries_a_shared_proton_across(config):
    """End to end on the same system, through the verdict the search actually installs."""
    from prrs.search import rung_verdict
    from prrs.chemistry import chemical_key

    def factory():
        return SharedProton()

    source = shared_proton_atoms(-SharedProton.width)
    source_components = chemical_key(source, config.bond_scale)
    frame = shared_proton_atoms(-0.05)  # just short of the ridge
    tangent = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 0]])  # push the proton across
    seen = []

    def name(endpoint, confirmed, checked):
        answer = rung_verdict(
            source,
            endpoint,
            source_components,
            chemical_key(endpoint, config.bond_scale),
            config,
        )
        seen.append(answer["verdict"])
        return {**answer, "stop": answer["verdict"] != "returned_to_source"}

    _, report = continue_across_boundary(frame, tangent, factory, config, verdict=name)
    assert report["outcome"] == "reached_degenerate_product"
    assert report["rungs"][-1]["broken"] and report["rungs"][-1]["formed"]
    # And the same walk driven the other way stays where it started.
    _, backward = continue_across_boundary(
        shared_proton_atoms(-0.05), -tangent, factory, config, verdict=name
    )
    assert backward["outcome"] == "exhausted"
    assert {rung["verdict"] for rung in backward["rungs"]} == {"returned_to_source"}


@pytest.mark.parametrize("builder", ["shared_proton", "rotor"])
def test_the_test_potentials_forces_match_their_energies(builder):
    """A test potential whose forces are not minus its energy gradient tests the optimizer
    against itself. The sign was wrong in both of these and invisible: the rotor's bonds
    start at equilibrium, so its bond forces were zero, and the shared proton only showed it
    by quenching onto its own barrier top."""
    if builder == "shared_proton":
        atoms, calculator = shared_proton_atoms(-0.12), SharedProton()
    else:
        rotor = pytest.importorskip("test_soft_mode_tiers")
        atoms, calculator = rotor.rotor_atoms(np.pi / 6), rotor.Rotor(amplitude=0.05)
    atoms.calc = calculator
    analytic = atoms.get_forces()
    numeric = np.zeros_like(analytic)
    step = 1e-6
    for atom in range(len(atoms)):
        for axis in range(3):
            for sign in (1, -1):
                shifted = atoms.copy()
                shifted.calc = calculator
                positions = shifted.positions.copy()
                positions[atom, axis] += sign * step
                shifted.set_positions(positions)
                numeric[atom, axis] += -sign * shifted.get_potential_energy() / (2 * step)
    assert np.abs(analytic - numeric).max() < 1e-6


class ProtonWithBumpyProduct(SharedProton):
    """The shared proton, plus a small ridge inside the product basin.

    This is the 3-oxobutanal situation reduced to three atoms. Every rung of every
    continuation in round four quenched onto a rotor saddle 11 meV above the product's
    minimum -- at -65 to -79 cm^-1 while the proton had demonstrably crossed -- so no rung
    ever handed a minimum to the caller and A3 stayed at zero. A saddle sits on the boundary
    of two basins, so descending both of its sides is what decides whether the rung arrived.
    """

    bump = 0.06
    bump_at = 0.4

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.positions
        step = 1e-6

        def q(pos):
            return float(np.linalg.norm(pos[1] - pos[0])) - float(
                np.linalg.norm(pos[2] - pos[1])
            )

        def extra(value):
            # A narrow gaussian ridge sitting on the product side only.
            return self.bump * np.exp(-(((value - self.bump_at) / 0.06) ** 2))

        self.results["energy"] += extra(q(positions))
        gradient = np.zeros_like(positions)
        for atom in range(len(atoms)):
            for axis in range(3):
                shifted = positions.copy()
                shifted[atom, axis] += step
                plus = extra(q(shifted))
                shifted[atom, axis] -= 2 * step
                gradient[atom, axis] = (plus - extra(q(shifted))) / (2 * step)
        self.results["forces"] = self.results["forces"] - gradient


def _namer(source, config):
    """The verdict callback the search installs, built for a given source structure."""
    from prrs.search import rung_verdict
    from prrs.chemistry import chemical_key

    source_components = chemical_key(source, config.bond_scale)

    def name(endpoint, confirmed, checked):
        answer = rung_verdict(
            source,
            endpoint,
            source_components,
            chemical_key(endpoint, config.bond_scale),
            config,
        )
        return {**answer, "stop": answer["verdict"] != "returned_to_source"}

    return name


def _saddle_and_mode(atoms, factory, config):
    """A stationary point and its unstable mode, via the same check the search uses."""
    from prrs.runner import confirm_minimum

    atoms = atoms.copy()
    atoms.calc = factory()
    confirmed, checked = confirm_minimum(atoms, factory, config, 0)
    assert not confirmed and checked["saddle_order"] == 1, (
        f"this test needs an index-one saddle, got confirmed={confirmed} "
        f"order={checked.get('saddle_order')}"
    )
    return atoms, checked


def test_a_saddle_whose_two_descents_agree_means_arrival(config):
    """The round-four blocker as a rule. The bump turns the product minimum into a ridge, so
    a quench that stops there is not at a minimum -- but both of its descents are the same
    product, so the geometry is inside the product's basin and did arrive."""
    from prrs.runner import arrival_through_saddle

    def factory():
        return ProtonWithBumpyProduct()

    source = shared_proton_atoms(-SharedProton.width)
    ridge, checked = _saddle_and_mode(
        shared_proton_atoms(ProtonWithBumpyProduct.bump_at), factory, config
    )
    agreed, sides = arrival_through_saddle(
        ridge, checked["unstable_mode"], factory, config, _namer(source, config)
    )
    assert agreed is not None, f"descents disagreed: {[s.get('verdict') for s in sides]}"
    assert agreed["verdict"] == "reached_degenerate_product"
    assert [side["verdict"] for side in sides] == ["reached_degenerate_product"] * 2
    assert agreed["broken"] and agreed["formed"]


def test_a_saddle_whose_two_descents_disagree_is_left_as_a_ridge(config):
    """The case the rule must NOT smooth over. The double well's barrier top joins two
    different chemical states, so a geometry sitting on it has not arrived anywhere and
    saying it had would be exactly the fabrication the continuation exists to avoid."""
    from prrs.runner import arrival_through_saddle

    top = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])  # r = 1.8
    ridge, checked = _saddle_and_mode(top, double_well_factory, config)
    source = Atoms("C2", positions=[[-0.6, 0, 0], [0.6, 0, 0]])  # the bonded minimum
    agreed, sides = arrival_through_saddle(
        ridge, checked["unstable_mode"], double_well_factory, config, _namer(source, config)
    )
    assert agreed is None
    assert {side["verdict"] for side in sides} == {"returned_to_source", "reached_other_state"}


def test_the_descent_rule_can_be_switched_off(config):
    """With the rule off, a rung that stops on a ridge reports only that it stopped, which
    is what round four measured and why its A3 read zero."""
    from dataclasses import replace

    plain = replace(
        config, boundary_continuation_descend_saddles=False, boundary_continuation_max_steps=2
    )

    def factory():
        return ProtonWithBumpyProduct()

    source = shared_proton_atoms(-SharedProton.width)
    frame = shared_proton_atoms(0.20)
    tangent = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 0]])
    _, report = continue_across_boundary(
        frame, tangent, factory, plain, verdict=_namer(source, plain)
    )
    assert all("arrival" not in rung for rung in report["rungs"])


def test_the_descent_relays_through_intervening_saddles(config):
    """`descend_saddle` is shared with the saddle registry's two-sided descent, so the relay
    it does has to be pinned here too: nine of round four's ten stalled saddles had one side
    reported as "not a minimum" precisely because a single step landed on another ridge."""
    from prrs.runner import descend_saddle

    top = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])
    ridge, checked = _saddle_and_mode(top, double_well_factory, config)
    landed = {}
    for sign in (1, -1):
        endpoint, side, relays = descend_saddle(
            ridge, checked["unstable_mode"], double_well_factory, config, sign
        )
        assert endpoint is not None, f"sign {sign}: {side.get('reason')}"
        assert side.get("saddle_order") == 0
        assert isinstance(relays, list)
        landed[sign] = round(float(endpoint.get_distance(0, 1)), 2)
    assert sorted(landed.values()) == [1.2, 2.4]


def test_a_walk_that_never_arrives_hands_back_nothing(config):
    """The P2 round-five defect. The naming callback is run on BOTH descents of a stalled
    rung, so one side can name the product while the other fails -- and that rung's own
    verdict is "stopped_on_saddle, not evidence of arrival". An earlier version let the
    callback stash the endpoint in an enclosing dict, and the search admitted it: two of the
    three continuations on that run's directed edge came from rungs that had not arrived.

    A callback run speculatively must not be able to record anything, and a walk whose
    outcome is not an arrival must hand back None.
    """
    stashed = []

    def greedy(endpoint, confirmed, checked):
        # A callback that tries to record. It must have no effect on what comes back.
        stashed.append(endpoint)
        return {"verdict": "returned_to_source"}

    frame = ridge_frame()
    tangent = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    endpoint, report = continue_across_boundary(
        frame, -tangent, double_well_factory, config, verdict=greedy
    )
    assert report["outcome"] == "exhausted"
    assert stashed, "the callback must actually have been called for this to mean anything"
    assert endpoint is None

    # And a callback that names an arrival gets the endpoint back, so the None above is not
    # simply the walk never returning anything.
    def honest(endpoint_, confirmed, checked):
        r = float(endpoint_.get_distance(0, 1))
        return {
            "verdict": "reached_other_state"
            if abs(r - 2.4) < abs(r - 1.2)
            else "returned_to_source"
        }

    arrived, forward = continue_across_boundary(
        frame, tangent, double_well_factory, config, verdict=honest
    )
    assert forward["outcome"] == "reached_other_state"
    assert arrived is not None
