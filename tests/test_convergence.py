"""Mode-resolved convergence: a force tolerance does not bound the geometry error.

In the linear regime the offset along a coordinate is |g_q| / k_q, so one force tolerance
cannot serve a bond stretch at tens of eV/A^2 and a torsion at a fifth of an eV/rad^2.
These tests use an analytic chain whose torsional stiffness is set explicitly, so the
three tiers and the one-degree criterion can be checked without a model.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from prrs import internal
from prrs.config import SearchConfig
from prrs.reliability import GuardedCalculator
from prrs.runner import polish_soft_modes, quench, torsion_response

BONDS = ((0, 1), (1, 2), (2, 3))
ANGLES = ((0, 1, 2), (1, 2, 3))
TORSION = (0, 1, 2, 3)


class Chain(Calculator):
    """Four atoms: stiff bonds and angles, and a torsion of chosen barrier V.

    The torsion term is (V/2)(1 - cos 3 phi), so at a minimum its stiffness is
    d2E/dphi2 = 9V/2 exactly, which is what the tests compare against.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        barrier_eV=0.05,
        bond_k=20.0,
        angle_k=5.0,
        bond_r0=1.53,
        angle_t0=np.radians(112.0),
    ):
        super().__init__()
        self.barrier, self.bond_k, self.angle_k = barrier_eV, bond_k, angle_k
        self.bond_r0, self.angle_t0 = bond_r0, angle_t0

    @property
    def torsional_stiffness(self):
        return 4.5 * self.barrier

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.positions
        energy, forces = 0.0, np.zeros_like(positions)
        for pair in BONDS:
            value = internal.coordinate(positions, "bond", pair)
            energy += 0.5 * self.bond_k * (value - self.bond_r0) ** 2
            forces -= (
                self.bond_k
                * (value - self.bond_r0)
                * internal.gradient(positions, "bond", pair)
            )
        for triple in ANGLES:
            value = internal.coordinate(positions, "angle", triple)
            energy += 0.5 * self.angle_k * (value - self.angle_t0) ** 2
            forces -= (
                self.angle_k
                * (value - self.angle_t0)
                * internal.gradient(positions, "angle", triple)
            )
        phi = internal.coordinate(positions, "dihedral", TORSION)
        energy += 0.5 * self.barrier * (1 - np.cos(3 * phi))
        forces -= (
            1.5
            * self.barrier
            * np.sin(3 * phi)
            * internal.gradient(positions, "dihedral", TORSION)
        )
        self.results = {"energy": float(energy), "forces": forces}


def _chain(phi_deg, barrier=0.05):
    """A chain built at a chosen torsion, with bonds and angles already at their minima."""
    r, theta = 1.53, np.radians(112.0)
    a = np.array([0.0, 0.0, 0.0])
    b = a + [r, 0.0, 0.0]
    c = b + [-r * np.cos(theta), r * np.sin(theta), 0.0]
    axis = (c - b) / np.linalg.norm(c - b)
    # d sits at angle theta from the c->b direction, so its component along b->c is
    # -r cos(theta); the reference arm points the same way as a, making phi_deg = 0 syn.
    arm = a - b
    arm -= np.dot(arm, axis) * axis
    arm = arm / np.linalg.norm(arm) * r * np.sin(theta)
    phi = np.radians(phi_deg)
    d = c - axis * (r * np.cos(theta)) + arm * np.cos(phi) + np.cross(axis, arm) * np.sin(phi)
    atoms = Atoms("C4", positions=[a, b, c, d])
    atoms.calc = Chain(barrier)
    return atoms


def test_torsion_response_recovers_the_analytic_stiffness():
    atoms = _chain(0.0)
    torsion = internal_torsion(atoms)
    gradient_, curvature = torsion_response(atoms, atoms.calc, torsion, 0.02)
    assert gradient_ == pytest.approx(0.0, abs=1e-6)
    assert curvature == pytest.approx(atoms.calc.torsional_stiffness, rel=1e-3)


def internal_torsion(atoms):
    from prrs.chemistry import canonical_labels
    from prrs.perturbations import rotatable_torsions
    from prrs.state import encode

    graph = encode(atoms)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    assert len(genuine + rotors) == 1
    return (genuine + rotors)[0]


@pytest.mark.parametrize("barrier,tier", [(0.05, "soft"), (20.0, "stiff"), (0.0, "free")])
def test_the_three_tiers_are_classified_by_curvature(barrier, tier):
    atoms = _chain(0.0, barrier)
    report = polish_soft_modes(atoms, atoms.calc, SearchConfig())
    assert [mode["tier"] for mode in report["modes"]] == [tier]


def test_force_tolerance_alone_leaves_a_soft_torsion_off_and_the_polish_fixes_it():
    """The whole point, in one comparison.

    Both structures satisfy the same force tolerance. Only the polished one is actually
    at a stationary geometry along the soft coordinate.
    """
    from ase.optimize import FIRE

    config = SearchConfig(quench_fmax_eV_A=0.02, quench_steps=500)
    torsion_offsets = {}
    for label in ("fmax_only", "polished"):
        atoms = _chain(28.0)  # well inside the 60 degree basin
        if label == "fmax_only":
            assert FIRE(atoms, logfile=None).run(fmax=config.quench_fmax_eV_A, steps=500)
        else:
            converged, diagnostics = quench(atoms, atoms.calc, config)
            assert converged, diagnostics
        assert np.linalg.norm(atoms.get_forces(), axis=1).max() <= config.quench_fmax_eV_A
        torsion = internal_torsion(atoms)
        gradient_, curvature = torsion_response(atoms, atoms.calc, torsion, 0.02)
        torsion_offsets[label] = abs(gradient_ / curvature)

    assert torsion_offsets["polished"] < config.torsion_tolerance_rad
    assert torsion_offsets["fmax_only"] > torsion_offsets["polished"]
    # The point is that the loose one is loose in degrees, not just relatively worse.
    assert np.degrees(torsion_offsets["fmax_only"]) > 1.0


def test_a_free_coordinate_is_recorded_rather_than_demanded():
    """A flat coordinate cannot be converged and must not block convergence."""
    atoms = _chain(17.0, barrier=0.0)
    converged, diagnostics = quench(atoms, atoms.calc, SearchConfig())
    assert converged
    modes = diagnostics["rounds"][0]["modes"]
    assert modes[0]["tier"] == "free"
    assert "carries no information" in modes[0]["note"]


def test_polish_leaves_the_stiff_coordinates_where_they_were():
    """Soft-mode polishing must not be a second global optimization."""
    config = SearchConfig(quench_fmax_eV_A=0.02)
    atoms = _chain(28.0)
    assert quench(atoms, atoms.calc, config)[0]
    for pair in BONDS:
        assert internal.coordinate(atoms.positions, "bond", pair) == pytest.approx(
            atoms.calc.bond_r0, abs=0.01
        )
    for triple in ANGLES:
        assert internal.coordinate(atoms.positions, "angle", triple) == pytest.approx(
            atoms.calc.angle_t0, abs=0.01
        )


def _round(movements, converged=True):
    """A polish report shaped like the real one, for the pure contraction tests."""
    return {
        "converged": converged,
        "moved": bool(movements),
        "modes": [],
        "movements": [
            {"indices": list(indices), "offset_rad": float(offset), "kind": "newton"}
            for indices, offset in movements
        ],
    }


def test_polish_contraction_reads_the_coordinate_that_moved_in_every_round():
    """Which coordinate carries a failure is not "the only one that moved".

    The offsets are the measured ones from the P3 benzoylacetone descent: a methyl torsion
    moving one way by a shrinking amount over seven rounds, ratio 0.834. A phenyl torsion
    moved in the first round only. Keying on "exactly one coordinate moved" -- which an
    earlier version of this did -- lets that one transient appearance discard the
    coordinate that is actually carrying the failure, and the whole verdict with it.
    """
    from prrs.runner import polish_contraction

    methyl, phenyl = (14, 0, 1, 3), (3, 4, 6, 11)
    offsets = [0.13449, 0.12360, 0.09873, 0.07535, 0.06764, 0.05465, 0.04534]
    rounds = [_round([(methyl, offsets[0]), (phenyl, 0.04407)])]
    rounds += [_round([(methyl, value)]) for value in offsets[1:]]

    report = polish_contraction(rounds, SearchConfig())
    assert report["contracting"] and report["verdict"] == "contracting"
    assert report["coordinate"] == list(methyl)
    assert report["transient"] == [list(phenyl)]
    assert report["offsets"] == offsets
    assert report["decay_per_round"] == pytest.approx(0.834, abs=0.005)
    # Short by how much, as a number: 0.0453 rad decaying at 0.834 needs about five more
    # rounds to pass one degree, against a budget of six.
    assert report["extra_rounds_to_tolerance"] == pytest.approx(5.25, abs=0.1)


def test_polish_contraction_calls_an_alternating_coordinate_a_limit_cycle():
    """The case no round count fixes has to stay distinguishable from the one it does."""
    from prrs.runner import polish_contraction

    torsion = (0, 1, 2, 3)
    rounds = [_round([(torsion, value)]) for value in (0.09, -0.08, 0.075, -0.07)]
    report = polish_contraction(rounds, SearchConfig())
    assert not report["contracting"]
    assert report["verdict"] == "limit_cycle"
    # Same-sign but not shrinking is also not a contraction, and is named separately.
    rounds = [_round([(torsion, value)]) for value in (0.09, 0.09, 0.09)]
    assert polish_contraction(rounds, SearchConfig())["verdict"] == "same_sign_not_shrinking"


class _PolishThatKeepsMoving:
    """A polish whose own criterion passes while it keeps displacing the structure.

    The surface-level fight this stands in for -- FIRE settling every coordinate, the
    torsion polish then re-optimizing one of them at frozen Cartesian rest, each undoing
    part of the other -- is a property of a coupled potential and is measured, on the P3
    benzoylacetone descents, at 0.834 per round over seven rounds. What is under test here
    is the control flow that reads that sequence, so the sequence is supplied directly
    rather than rebuilt out of a synthetic potential that happens to fight the same way.

    `sign` flips the offsets into a limit cycle, which must NOT be admitted.
    """

    def __init__(self, sign=1, first=0.12, decay=0.85):
        self.calls, self.sign, self.first, self.decay = 0, sign, first, decay

    def __call__(self, atoms, guard, config):
        offset = self.first * self.decay**self.calls * (self.sign**self.calls)
        self.calls += 1
        # A real displacement, so the force tolerance the caller was standing on is
        # genuinely invalidated and a re-minimize has something to do.
        atoms.positions[-1] += [0.0, 0.0, 0.02]
        return _round([((0, 1, 2, 3), offset)])


def test_quench_reminimizes_after_a_trailing_polish_that_moved(monkeypatch):
    """The geometry a failed quench hands back has to satisfy the tolerance it claims.

    The success path returns from inside the loop precisely when the polish left the
    structure alone, so it never had this hole. The failure path ran one last polish and
    returned whatever that displaced the structure to, with no minimization after it:
    measured on the P3 descents, fmax = 2.1e-2 against a 2e-3 criterion, on a geometry a
    caller then read a Hessian off and a chemical key off. Both readings want a structure
    that is where the force tolerance says it is.
    """
    from prrs import runner

    monkeypatch.setattr(runner, "polish_soft_modes", _PolishThatKeepsMoving())
    config = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.002, soft_polish_rounds=6)
    atoms = _chain(30.0)
    guard = GuardedCalculator(Chain(0.05), config, None)
    atoms.calc = guard

    converged, diagnostics = runner.quench(atoms, guard, config)

    assert not converged
    assert diagnostics["reason"] == "soft_mode_moved_on_the_last_round"
    assert diagnostics["force_converged"] is True
    assert np.abs(atoms.get_forces()).max() <= config.quench_fmax_eV_A
    assert diagnostics["polish_contraction"]["contracting"]


def test_a_contracting_soft_coordinate_no_longer_discards_the_descent_endpoint(monkeypatch):
    """The defect that cost P3 a two-sided connection six times out of six.

    `follow_unstable_mode` returned None whenever the quench did not reach its own fixed
    point, so an endpoint whose only unfinished business was a near-free rotor never
    reached the minimum check at all. The gate is `confirm_minimum`, not the polish; what
    this changes is that the endpoint gets to the gate.

    The admission is narrow and the second half of this test is the half that matters: an
    alternating coordinate is a limit cycle, no round count fixes it, and it is still
    refused.
    """
    from prrs import runner
    from prrs.calculators import double_well_factory

    config = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.002, soft_polish_rounds=6)
    saddle = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])  # the well's barrier top
    saddle.calc = double_well_factory()
    mode = np.array([[-1.0, 0, 0], [1.0, 0, 0]])

    monkeypatch.setattr(runner, "polish_soft_modes", _PolishThatKeepsMoving())
    endpoint, diagnostics = runner.follow_unstable_mode(
        saddle, mode, double_well_factory, config, 1
    )
    assert endpoint is not None, diagnostics
    assert diagnostics["admitted_with_unsettled_soft_modes"] is True
    assert diagnostics["reason"] == "soft_mode_moved_on_the_last_round"
    # Admitted, not converged: the record must not start claiming the quench finished.
    assert diagnostics["converged"] is False
    assert np.abs(endpoint.get_forces()).max() <= config.quench_fmax_eV_A

    monkeypatch.setattr(runner, "polish_soft_modes", _PolishThatKeepsMoving(sign=-1))
    refused, diagnostics = runner.follow_unstable_mode(
        saddle, mode, double_well_factory, config, 1
    )
    assert refused is None
    assert diagnostics["polish_contraction"]["verdict"] == "limit_cycle"
    assert "admitted_with_unsettled_soft_modes" not in diagnostics


def test_an_unreadable_direction_makes_the_minimum_check_unverified_not_refuted():
    """A direction the gate refused was not probed, and a probe that never ran cannot vote.

    The defect: `_probe_minimum` dropped gate-refused directions and then computed
    `confirmed = bool(measured) and not negative`. Four refusals and one non-negative
    reading came back `confirmed = True` -- the verdict rested on whichever directions
    happened to survive. And the caller reads `not confirmed`, looks at `saddle_order`,
    and writes `not_a_minimum`, so an unmeasurable structure was being recorded as a
    refuted one and its geometry dropped.

    Three outcomes now. `complete` says every direction was read; `confirmed` adds that
    none was negative; `verdict` names which of the three happened. Incomplete is neither
    a pass nor a refutation, and the caller keeps the candidate.
    """
    import prrs.runner as runner
    from ase import Atoms
    from prrs.calculators import double_well_factory
    from prrs.config import SearchConfig
    from prrs.reliability import GateRejected

    config = SearchConfig(minimum_check="probe", minimum_check_probes=3)
    candidate = Atoms("C2", positions=[[-0.6, 0, 0], [0.6, 0, 0]])
    candidate.calc = double_well_factory()

    real = runner._curvature_along
    calls = {"n": 0}

    def refuse_all_but_one(atoms, guard, direction, step):
        calls["n"] += 1
        if calls["n"] > 1:
            raise GateRejected("overlap", "synthetic refusal", atoms)
        return real(atoms, guard, direction, step)

    runner._curvature_along = refuse_all_but_one
    try:
        confirmed, report = runner.confirm_minimum(candidate, double_well_factory, config, 0)
    finally:
        runner._curvature_along = real

    assert report["method"] == "probe"
    assert report["unmeasured_directions"] > 0, "the synthetic refusals did not land"
    assert report["complete"] is False
    assert confirmed is False, "an unread direction must not be counted as a pass"
    assert report["verdict"] == "unverified", report["verdict"]
    assert report["unmeasured_reasons"], "the reason a direction went unread must survive"

    # And with every direction readable the verdict is a real one again.
    confirmed, report = runner.confirm_minimum(candidate, double_well_factory, config, 0)
    assert report["complete"] is True
    assert report["verdict"] in ("minimum", "not_a_minimum")
    assert confirmed is (report["verdict"] == "minimum")
