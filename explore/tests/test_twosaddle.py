import numpy as np
from prrs.config import SearchConfig
from prrs.io import snapshot
from prrs.reliability import GuardedCalculator
from prrs.runner import confirm_minimum, hessian_spectrum
from scipy.spatial.transform import Rotation
from twosaddle import SADDLE_Y, local_coords, two_saddle_atoms, two_saddle_factory

CONFIG = SearchConfig(closed_shell_only=False)


def _snap(x, y):
    atoms = two_saddle_atoms(x, y)
    return snapshot(atoms, GuardedCalculator(two_saddle_factory(), CONFIG))


def _moved(atoms):
    moved = atoms.copy()
    moved.positions = Rotation.from_rotvec([0.3, -1.1, 0.7]).apply(moved.positions) + [1, 2, 3]
    return moved


def test_energy_is_invariant_under_rigid_motion():
    atoms = two_saddle_atoms(0.3, 0.4, 0.1)
    atoms.calc = two_saddle_factory()
    e0 = atoms.get_potential_energy()
    moved = _moved(atoms)
    moved.calc = two_saddle_factory()
    assert abs(moved.get_potential_energy() - e0) < 1e-10


def test_local_coords_roundtrip_and_rigid_motion_invariance():
    atoms = two_saddle_atoms(0.3, -0.4, 0.1)
    assert np.allclose(local_coords(atoms.positions), (0.3, -0.4, 0.1), atol=1e-12)
    assert np.allclose(local_coords(_moved(atoms).positions), (0.3, -0.4, 0.1), atol=1e-12)


def test_both_saddles_have_exactly_one_negative_mode():
    for y in SADDLE_Y:
        eig = hessian_spectrum(_snap(0.0, y), two_saddle_factory, CONFIG)
        assert np.sum(eig < -CONFIG.minimum_check_eigenvalue_tol) == 1


def test_minima_are_confirmed_and_share_one_bond_graph():
    from prrs.state import encode

    graphs = []
    for x in (-1.0, 1.0):
        ok, _ = confirm_minimum(_snap(x, 0.0), two_saddle_factory, CONFIG, seed=0)
        assert ok
        graphs.append(encode(two_saddle_atoms(x, 0.0), CONFIG.bond_scale).edges)
    assert graphs[0] == graphs[1]


def test_minima_are_symmetry_images_but_saddles_are_distinct():
    from prrs.chemistry import automorphisms, chemical_key, symmetric_rmsd
    from prrs.state import same_basin

    lo, hi = _snap(-1.0, 0.0), _snap(1.0, 0.0)
    keys = [chemical_key(s, CONFIG.bond_scale)["key"] for s in (lo, hi)]
    assert keys[0] == keys[1]
    assert not same_basin(lo, hi, CONFIG)
    perms, info = automorphisms(lo, CONFIG.bond_scale)
    assert info["count"] == 6
    assert symmetric_rmsd(lo, hi, perms) < 1e-9
    s0, s1 = (_snap(0.0, y) for y in SADDLE_Y)
    assert symmetric_rmsd(s0, s1, perms) > CONFIG.basin_rmsd_A
