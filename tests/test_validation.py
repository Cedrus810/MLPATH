"""A second implementation of the same curvature, because one implementation proves nothing.

The search's curvature numbers are self-consistent: one routine projects the rigid-body
modes out, mass-weights and diagonalizes, and until now nothing else computed the same
quantity. These tests put ASE's vibrational machinery on the other side of the comparison.

Two checks that are not interchangeable. Feeding our Hessian to ASE's spectral code shares
the differentiation and so isolates the mass weighting, the ordering and the unit
conversion. Letting ASE build its own Hessian tests the differentiation too -- but only if
its step differs from ours, because at the same step central differences are the same
computation, and agreement then confirms the column assembly rather than the limit.
"""

import numpy as np
import pytest
from ase import Atoms
from prrs.calculators import double_well_factory
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, hessian_spectrum
from prrs import validation as v

CONFIG = dict(quench_fmax_eV_A=0.002)


def _pair(distance):
    atoms = Atoms("C2", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
    atoms.calc = double_well_factory()
    return atoms


def _chain(phi_deg):
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain as build

    return build(phi_deg), Chain


def _internal(atoms, factory, config, source="fd"):
    eigenvalues, _, _, _ = curvature_spectrum(atoms, factory, config, source=source)
    energies, trivial, _ = v.ase_energies_from_our_hessian(atoms, factory, config, source)
    ours = v.internal_eigenvalues(eigenvalues, trivial)
    theirs = v.internal_eigenvalues(v.eigenvalues_from_energies(energies), trivial)
    return ours, theirs, trivial


def test_the_same_hessian_gives_the_same_spectrum_through_ase():
    """Isolates the convention: mass weighting, mode ordering, unit conversion."""
    config = SearchConfig(**CONFIG)
    for atoms, factory in (
        (_pair(1.2), double_well_factory),
        (_pair(1.8), double_well_factory),
        _chain(0.0),
        _chain(60.0),
    ):
        ours, theirs, trivial = _internal(atoms, factory, config)
        assert len(ours) == 3 * len(atoms) - trivial
        assert np.allclose(ours, theirs, atol=1e-7), (ours, theirs)
        # And the sign structure survives, which is the part that classifies the point.
        assert np.sign(ours).tolist() == np.sign(theirs).tolist()


def test_a_diatomic_leaves_five_rigid_body_modes_not_six():
    config = SearchConfig(**CONFIG)
    _, _, trivial = _internal(_pair(1.2), double_well_factory, config)
    assert trivial == 5
    atoms, factory = _chain(0.0)
    _, _, trivial = _internal(atoms, factory, config)
    assert trivial == 6


def test_ases_own_hessian_at_our_step_confirms_the_column_assembly():
    """Same algorithm and same step, so agreement must be near-exact rather than close.

    This is what says our 3N Hessian-vector products are assembled into the matrix the way
    everyone else assembles theirs. It says nothing about the truncation, which is shared.
    """
    config = SearchConfig(**CONFIG)
    atoms, factory = _chain(60.0)
    ours, _, trivial = _internal(atoms, factory, config)
    energies, trivial_ase, provenance = v.ase_vibrations_energies(atoms, factory, config)
    theirs = v.internal_eigenvalues(v.eigenvalues_from_energies(energies), trivial_ase)
    assert provenance["delta_A"] == config.minimum_check_step_A
    assert provenance["force_evaluations"] == 6 * len(atoms)
    assert np.allclose(ours, theirs, atol=1e-6), (ours, theirs)


def test_disagreement_at_another_step_is_truncation_and_scales_like_one():
    """The sharp version: two central differences at different steps must differ by

        C |h1^2 - h2^2|

    for one constant C. A step-dependent disagreement that does not scale this way is a
    defect; one that does is the limit being approached from two places. Measured against
    our own step, so it also bounds the truncation error our default step carries.
    """
    config = SearchConfig(**CONFIG)
    atoms, factory = _chain(60.0)
    ours, _, _ = _internal(atoms, factory, config)
    reference_step = config.minimum_check_step_A

    implied = []
    for delta in (0.002, 0.05):
        energies, trivial, _ = v.ase_vibrations_energies(atoms, factory, config, delta=delta)
        theirs = v.internal_eigenvalues(v.eigenvalues_from_energies(energies), trivial)
        difference = np.abs(ours - theirs).max()
        assert difference > 1e-5, delta  # it really did change
        implied.append(difference / abs(delta**2 - reference_step**2))
    assert implied[0] == pytest.approx(implied[1], rel=0.15), implied


def test_the_wavenumber_conversion_round_trips():
    """Our icm and ASE's eV are the same quantity in two dresses."""
    from prrs.runner import _wavenumbers
    from ase import units

    config = SearchConfig(**CONFIG)
    atoms, factory = _chain(60.0)
    eigenvalues = hessian_spectrum(atoms, factory, config, source="fd")
    negative = [value for value in eigenvalues if value < -1e-9]
    assert len(negative) == 1
    icm = _wavenumbers(negative)[0]
    recovered = -(((abs(icm) * units.invcm) / v.ENERGY_PER_ROOT_EIGENVALUE) ** 2)
    assert recovered == pytest.approx(negative[0], rel=1e-9)
    assert v.eigenvalues_from_energies([abs(icm) * units.invcm * 1j])[0] == pytest.approx(
        negative[0], rel=1e-9
    )


def test_every_numeric_config_field_has_a_bound():
    """No numeric field may reach a run unvalidated.

    The bound lists drifted behind the dataclass once already: ten integers added for the
    min-mode walk, the response descent and the saddle search were never added to any
    list, and min_mode_refresh = 0 reached `step % config.min_mode_refresh` inside
    follow_min_mode as a plain ZeroDivisionError, an hour into a run. Adding those ten
    names fixes the ten; this test is what stops the eleventh. A new numeric field either
    goes in a bucket or is named in the exemption set with the branch that covers it.
    """
    from dataclasses import fields
    from prrs import config as cfg

    bucketed = (
        set(cfg.INTEGER_POSITIVE)
        | set(cfg.INTEGER_NONNEGATIVE)
        | set(cfg.FLOAT_POSITIVE)
        | set(cfg.FLOAT_NONNEGATIVE)
    )
    numeric = {
        f.name
        for f in fields(cfg.SearchConfig)
        if isinstance(f.default, (int, float)) and not isinstance(f.default, bool)
    }
    # A scalar field typed `int | None` defaults to None, so the default's type does not
    # reveal it. Match on the annotation with the None arm stripped, and require it to be
    # exactly int or float: `tuple[int, ...] | None` is a container with its own branch,
    # not a number that wants a bound.
    for field in fields(cfg.SearchConfig):
        annotation = str(field.type).replace("| None", "").replace("Optional", "").strip()
        if annotation in ("int", "float"):
            numeric.add(field.name)
    unbounded = numeric - bucketed - set(cfg._UNBUCKETED_NUMERIC)
    assert not unbounded, f"numeric config fields with no bound: {sorted(unbounded)}"

    # The exemption set must not outlive the fields it exempts, and a bucket must not
    # name a field the dataclass no longer has.
    names = {f.name for f in fields(cfg.SearchConfig)}
    assert set(cfg._UNBUCKETED_NUMERIC) <= names
    assert bucketed <= names
    # Positive and nonnegative are exclusive; naming a field in both hides which applies.
    assert not set(cfg.INTEGER_POSITIVE) & set(cfg.INTEGER_NONNEGATIVE)
    assert not set(cfg.FLOAT_POSITIVE) & set(cfg.FLOAT_NONNEGATIVE)


def test_every_rationale_entry_names_a_live_field():
    """The rationale table must not outlive the defaults it justifies.

    The reasons the non-obvious defaults are what they are used to sit as comment blocks
    inside the dataclass, where nothing could check them and where they broke the field
    list into fragments. They are the record of a measurement or a failure, so they moved
    out of the body rather than being cut. Out of the body they can go stale silently --
    a renamed field leaves its paragraph pointing at nothing -- which is what this pins.
    """
    from dataclasses import fields
    from prrs import config as cfg

    names = {f.name for f in fields(cfg.SearchConfig)}
    orphaned = set(cfg.RATIONALE) - names
    assert not orphaned, (
        f"RATIONALE entries for fields that no longer exist: {sorted(orphaned)}"
    )
    empty = [k for k, v in cfg.RATIONALE.items() if not v.strip()]
    assert not empty, f"empty rationale entries: {empty}"


def test_loop_and_budget_parameters_are_refused_at_construction():
    """The failure this replaces was a ZeroDivisionError deep inside follow_min_mode."""
    import pytest as _pytest
    from prrs.config import SearchConfig

    with _pytest.raises(ValueError):
        SearchConfig(min_mode_refresh=0)
    for name in (
        "min_mode_max_steps",
        "response_descent_max_steps",
        "saddle_search_seeds",
        "saddle_search_per_run",
        "irc_relay_hops",
        "min_pair_hops",
        "min_directions_per_family",
        "free_alignment_samples",
    ):
        with _pytest.raises(ValueError):
            SearchConfig(**{name: 0})
        with _pytest.raises(ValueError):
            SearchConfig(**{name: -3})
        with _pytest.raises(ValueError):
            SearchConfig(**{name: 2.5})
    # Zero backtracks means take the step or give up, which is a real configuration.
    assert SearchConfig(min_mode_backtracks=0).min_mode_backtracks == 0
    with _pytest.raises(ValueError):
        SearchConfig(min_mode_backtracks=-1)


def test_fractions_are_bounded_above():
    """A quota above the whole, an overlap above a unit dot product, a relative tolerance
    that accepts anything: each silently disables the check it parameterizes."""
    import pytest as _pytest
    from prrs.config import SearchConfig

    for name in ("direction_quota_fraction", "min_mode_overlap"):
        assert SearchConfig(**{name: 1.0})
        with _pytest.raises(ValueError):
            SearchConfig(**{name: 1.5})
    for name in ("response_descent_tol", "min_distance_ratio"):
        assert SearchConfig(**{name: 0.99})
        with _pytest.raises(ValueError):
            SearchConfig(**{name: 1.0})


def test_the_declared_charge_reaches_the_structure_or_the_run_refuses():
    """Two statements about the same charge may not disagree, and may not both be absent.

    `SearchConfig.total_charge` reaches `chemical_key` through network.Registry, so the
    identity layer was charge-aware; the potential hears about the charge only through
    `atoms.info["charge"]`, and no line in the package ever wrote that key. A run could
    declare an anion and evaluate the neutral species with nothing raising.

    It is not a bookkeeping nicety. Measured on Cl- + CH3Cl with MACE-POLAR-1-M,
    declaring -1 puts -0.9729 e on the nucleophilic chlorine and gives a dipole of 2.96;
    with the tag absent both chlorines sit near -0.03 and the dipole is 0.37. Two
    different physical objects, and the difference is silent.

    A charge-blind potential is unaffected -- MACE-OFF24 returns bit-identical energies
    with the tag absent, 0, or -1 -- so stamping is safe everywhere.
    """
    import pytest as _pytest
    from ase import Atoms
    from prrs.config import SearchConfig
    from prrs.search import _stamped_with_charge

    def water():
        return Atoms("OH2", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])

    neutral = SearchConfig()
    anion = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)

    # An unstamped structure is stamped from the configuration, both ways round.
    assert _stamped_with_charge(water(), neutral).info["charge"] == 0
    assert _stamped_with_charge(water(), anion).info["charge"] == -1
    assert _stamped_with_charge(water(), anion).info["spin"] == 1

    # An agreeing structure passes through unchanged.
    agreeing = water()
    agreeing.info["charge"], agreeing.info["spin"] = -1, 1
    assert _stamped_with_charge(agreeing, anion).info["charge"] == -1

    # Disagreement is refused rather than silently resolved in either direction.
    conflicting = water()
    conflicting.info["charge"] = -1
    with _pytest.raises(ValueError, match="two answers to the same question"):
        _stamped_with_charge(conflicting, neutral)

    spin_conflict = water()
    spin_conflict.info["spin"] = 3
    with _pytest.raises(ValueError, match="two answers to the same question"):
        _stamped_with_charge(spin_conflict, anion)

    # The calculator survives the stamping; losing it would break the caller silently.
    from prrs.calculators import double_well_factory

    carried = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    carried.calc = double_well_factory()
    assert _stamped_with_charge(carried, neutral).calc is carried.calc
