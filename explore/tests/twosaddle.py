"""Rotation- and translation-invariant analytic test potential with two saddles.

This is a TEST potential, not chemistry. Three O anchors form a rigid triangle (harmonic
edges); the mobile H feels, in the anchor frame (x, y, z),

    E = (x^2 + y^2 - 1)^2 + a y^2 + b y^3 + z_k z^2 / 2 + frame term,

the same field as the frozen bench `benchmarks/pars_two_saddles` (a, b, z_k), but written
in coordinates carried by the anchors, so the energy is invariant under rigid motion.
Minima at (x, y) = (+-1, 0); first-order saddles at x = 0, y in SADDLE_Y, joining them.
Forces are central differences of the energy (step H_STEP).

EDGE = 4.0 keeps H-O beyond the O-H bond cutoff (1.164 A) at every stationary point, so
both minima have the same (empty) bond graph.

Symmetry caveat: with an empty graph the three O are equivalent, `automorphisms` has 6
elements, `symmetric_rmsd(min(-1), min(+1)) == 0` and the two minima share a
`chemical_key`: they are symmetry images. `prrs.state.same_basin` is False for them only
because the label-based `aligned_rmsd` is 0.866 A > basin_rmsd_A, so C2 on this potential
relies on label-based same_basin. The two saddles do differ: symmetric_rmsd ~ 0.33 A vs
basin_rmsd_A 0.15 A.
"""

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

SADDLE_Y = (-0.8051774, 0.7451774)
EDGE = 4.0
H_STEP = 1e-6


def _frame(anchors):
    c = anchors.mean(axis=0)
    e1 = anchors[1] - anchors[0]
    e1 = e1 / np.linalg.norm(e1)
    v = anchors[2] - anchors[0]
    v = v - (v @ e1) * e1
    e2 = v / np.linalg.norm(v)
    return c, e1, e2, np.cross(e1, e2)


def local_coords(positions):
    """(x, y, z) of H (atom 3) in the anchor frame; independent of rigid motion."""
    p = np.asarray(positions, dtype=float)
    c, e1, e2, e3 = _frame(p[:3])
    d = p[3] - c
    return float(d @ e1), float(d @ e2), float(d @ e3)


def two_saddle_atoms(x, y, z=0.0):
    h = EDGE * np.sqrt(3.0) / 2.0
    anchors = np.array(
        [[-EDGE / 2, -h / 3, 0.0], [EDGE / 2, -h / 3, 0.0], [0.0, 2 * h / 3, 0.0]]
    )
    return Atoms("OOOH", positions=np.vstack([anchors, [x, y, z]]))


class TwoSaddleInvariant(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, a=0.8, b=0.08, z_k=2.0, frame_k=10.0, step=H_STEP, **kwargs):
        super().__init__(**kwargs)
        self.a, self.b, self.z_k, self.frame_k, self.step = a, b, z_k, frame_k, step

    def _energy(self, p):
        x, y, z = local_coords(p)
        e = (x * x + y * y - 1.0) ** 2 + self.a * y * y + self.b * y**3 + 0.5 * self.z_k * z * z
        for i, j in ((0, 1), (0, 2), (1, 2)):
            e += 0.5 * self.frame_k * (np.linalg.norm(p[i] - p[j]) - EDGE) ** 2
        return float(e)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        p = self.atoms.positions.copy()
        forces = np.zeros_like(p)
        for i in range(p.shape[0]):
            for k in range(3):
                q = p.copy()
                q[i, k] += self.step
                ep = self._energy(q)
                q[i, k] -= 2 * self.step
                forces[i, k] = -(ep - self._energy(q)) / (2 * self.step)
        self.results = {"energy": self._energy(p), "forces": forces}


def two_saddle_factory():
    return TwoSaddleInvariant()
