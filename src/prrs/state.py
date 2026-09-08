"""Atom-mapped molecular graphs and basin comparisons over a declared active subset.

Everything the search compares -- graphs, coordination, RMSD, basin identity -- is
restricted to the active subset and reported with global atom indices. The MVP sets
the subset to every atom, but keeping the distinction explicit is what allows a
reactive region to later sit inside a large MM or QM/MM environment without
rewriting the analysis layer: solvent reorganisation must not make two visits to the
same basin look like different basins, and an all-pairs distance matrix over
thousands of environment atoms is neither affordable nor meaningful.
"""

from dataclasses import dataclass
import numpy as np
from ase.data import covalent_radii


def resolve_active(atoms, active=None):
    """Normalize an active-subset selection into a sorted array of global indices."""
    count = len(atoms)
    if active is None:
        return np.arange(count)
    index = np.asarray(sorted(set(int(i) for i in active)), dtype=int)
    if index.size < 2:
        raise ValueError("The active subset needs at least two atoms")
    if index[0] < 0 or index[-1] >= count:
        raise ValueError("Active atom indices fall outside the structure")
    return index


def validate_atoms(atoms, active=None):
    if len(atoms) < 2:
        raise ValueError("PRRS requires at least two atoms")
    # Scope gate, not an implementation detail: the graph below uses plain
    # displacements, so a periodic cell would be silently mis-measured, not refused.
    if np.any(atoms.pbc):
        raise ValueError(
            "MVP supports finite nonperiodic molecules only; PBC is not implemented"
        )
    if atoms.constraints:
        raise ValueError("Constraints are not supported in this MVP")
    if not np.isfinite(atoms.positions).all():
        raise ValueError("Non-finite input positions")
    masses = atoms.get_masses()
    if not np.isfinite(masses).all() or np.any(masses <= 0):
        raise ValueError("Atomic masses must be finite and positive")
    if np.any(atoms.numbers <= 0):
        raise ValueError("Dummy atoms are not supported")
    resolve_active(atoms, active)


@dataclass
class Graph:
    """Edges carry global atom indices; weights and distances are subset-local."""

    edges: set[tuple[int, int]]
    weights: np.ndarray
    distances: np.ndarray
    index: np.ndarray

    @property
    def coordination(self):
        return self.weights.sum(axis=1)

    @property
    def components(self):
        position = {int(atom): i for i, atom in enumerate(self.index)}
        parent = list(range(len(self.index)))

        def root(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for i, j in self.edges:
            parent[root(position[i])] = root(position[j])
        return len({root(i) for i in parent})


def encode(atoms, scale=1.2, previous=None, hysteresis=0.08, active=None):
    index = resolve_active(atoms, active)
    positions = atoms.positions[index]
    distances = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=-1)
    radii = covalent_radii[atoms.numbers[index]]
    cutoff = scale * (radii[:, None] + radii[None, :])
    weights = 1.0 / (1.0 + np.exp(np.clip(12 * (distances / cutoff - 1), -60, 60)))
    np.fill_diagonal(weights, 0.0)
    edges = set()
    for i in range(len(index)):
        for j in range(i + 1, len(index)):
            pair = (int(index[i]), int(index[j]))
            factor = 1.0
            if previous is not None:
                factor += hysteresis if pair in previous else -hysteresis
            if distances[i, j] < cutoff[i, j] * factor:
                edges.add(pair)
    return Graph(edges, weights, distances, index)


def aligned_rmsd(a, b, active=None):
    """Kabsch RMSD over the active subset, translation and proper-rotation invariant."""
    if not np.array_equal(a.numbers, b.numbers):
        return float("inf")
    index = resolve_active(a, active)
    x = a.positions[index] - a.positions[index].mean(axis=0)
    y = b.positions[index] - b.positions[index].mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((x @ rotation - y) ** 2, axis=1))))


def same_basin(a, b, config):
    active = config.active_atoms
    if (
        encode(a, config.bond_scale, active=active).edges
        != encode(b, config.bond_scale, active=active).edges
    ):
        return False
    return (
        abs(a.get_potential_energy() - b.get_potential_energy()) <= config.basin_energy_eV
        and aligned_rmsd(a, b, active) <= config.basin_rmsd_A
    )


def classify(source, endpoint, config):
    active = config.active_atoms
    ga = encode(source, config.bond_scale, active=active)
    gb = encode(endpoint, config.bond_scale, active=active)
    broken, formed = ga.edges - gb.edges, gb.edges - ga.edges
    if broken or formed:
        label = "dissociation" if gb.components > ga.components else "reactive"
    else:
        label = "elastic" if same_basin(source, endpoint, config) else "rearrangement"
    fingerprint = {
        "broken": [list(e) for e in sorted(broken)],
        "formed": [list(e) for e in sorted(formed)],
        "active_subset": ga.index.tolist(),
        "delta_coordination": (gb.coordination - ga.coordination).tolist(),
        "active_atoms": sorted({i for pair in broken | formed for i in pair}),
        "aligned_rmsd_A": aligned_rmsd(source, endpoint, active),
        "energy_change_eV": endpoint.get_potential_energy() - source.get_potential_energy(),
        "components_before": ga.components,
        "components_after": gb.components,
    }
    return label, fingerprint
