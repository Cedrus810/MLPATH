"""Coupled ring probes preserve other bonds without relaxing the damage gate."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from prrs.perturbations import Probe, apply, delivered_direction
from prrs.state import encode


@pytest.mark.parametrize(
    "family,indices,amplitude",
    [("stretch", (0, 1), 0.3), ("bend", (0, 1, 2), 0.17), ("torsion", (0, 1, 2, 3), 0.25)],
)
def test_ring_delivery_reaches_target_without_straining_other_bonds(family, indices, amplitude):
    source = molecule("C6H6")
    moved = source.copy()
    delivered = apply(moved, Probe(family, indices, 1, amplitude, 0))
    assert delivered["target_error"] < 1e-7
    assert delivered["collateral"]["max_bond_change_A"] < 1e-7
    assert moved.get_center_of_mass() == pytest.approx(source.get_center_of_mass(), abs=1e-9)


def test_ring_score_uses_the_constrained_delivery_tangent():
    source = molecule("C6H6")
    moved = []
    for sign in (1, -1):
        atoms = source.copy()
        apply(atoms, Probe("stretch", (0, 1), sign, 1e-4, 0))
        moved.append(atoms.positions)
    tangent = ((moved[0] - moved[1]) * np.sqrt(source.get_masses())[:, None]).ravel()
    closed = delivered_direction(source, "stretch", (0, 1), graph=encode(source))
    assert abs(float(closed @ tangent)) / np.linalg.norm(tangent) == pytest.approx(
        1.0, abs=1e-7
    )


def test_an_impossible_constraint_set_fails_instead_of_tearing_bonds():
    from prrs.internal import constrained_displace

    atoms = Atoms("C3", positions=[(0, 0, 0), (1, 0, 0), (0.5, np.sqrt(3) / 2, 0)])
    with pytest.raises(ValueError, match="did not converge"):
        constrained_displace(atoms, "bond", (0, 1), 2.0, {(0, 1), (0, 2), (1, 2)})
