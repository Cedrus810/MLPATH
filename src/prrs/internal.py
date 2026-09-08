"""Primitive internal coordinates and the perturbations that act along them.

The proposal basis must span the internal-coordinate subspaces a molecule actually
moves in -- stretch, bend, torsion -- not only Cartesian pair directions. Ethanol
makes the point: its lowest nontrivial motion is a torsion, so a search built solely
from pair displacements and pair kicks cannot reach its other basins at any amplitude.

Every perturbation here is generated from the gradient of one internal coordinate.
Because an internal coordinate is invariant under rigid translation and rotation, its
gradient g satisfies

    sum_i g_i = 0        and        sum_i x_i x g_i = 0

so a momentum kick dp_i = alpha * g_i conserves total linear and angular momentum
exactly, and a displacement dx_i = s * g_i / m_i preserves the center of mass exactly.
Both follow from the coordinate, not from a projection applied afterwards, and neither
needs the molecule split into two rigid groups -- a split that is simply undefined for
a dihedral inside a ring.

Gradients are taken by central differences of the coordinate function itself. That
function is analytic and cheap, involves no potential evaluation, and reaches ~1e-9
accuracy, so the expensive resource is untouched and the result is easy to verify.
"""

import numpy as np

KINDS = {"bond": 2, "angle": 3, "dihedral": 4}
UNITS = {"bond": "Angstrom", "angle": "radian", "dihedral": "radian"}
PERIODIC = {"bond": False, "angle": False, "dihedral": True}
# Largest single increment, in Angstrom or radian. A displacement along the gradient of
# an angular coordinate is tangent to the arc, so one long step travels the chord: for a
# 2.09 rad torsion of ethanol's hydroxyl it stretched the O-H bond from 0.96 to 2.35 A
# while hitting the requested dihedral exactly. The named coordinate being exact is no
# evidence that the others survived, so large displacements are walked in increments
# with the gradient refreshed each time. The cost is geometry, never the potential.
MAX_STEP = {"bond": 0.05, "angle": 0.05, "dihedral": 0.05}
_STEP = 1e-6


def atom_index(value, count=None):
    """One atom index, validated rather than coerced.

    `int(value)` accepts everything and complains about nothing: -1 becomes a valid
    numpy index onto the last atom, and 1.9 truncates to 1. Both then perturb an atom
    nobody named, and the coordinate that comes back is a real number for the wrong
    triple, so nothing downstream can tell. numpy integers have to pass -- indices
    routinely arrive from np.where or argmax on a bond graph -- so the test is on the
    value being integral, not on `type(value) is int`.
    """
    if isinstance(value, bool):
        raise ValueError(f"atom index must be an integer, got {value!r}")
    try:
        integral = int(value) == value
    except (TypeError, ValueError):  # non-numeric, or a float NaN/inf
        integral = False
    if not integral:
        raise ValueError(f"atom index must be an integer, got {value!r}")
    index = int(value)
    if index < 0:
        raise ValueError(f"atom index must be nonnegative, got {index}")
    if count is not None and index >= count:
        raise ValueError(f"atom index {index} is out of range for {count} atoms")
    return index


def _check(kind, indices, count=None):
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {sorted(KINDS)}")
    if len(indices) != KINDS[kind] or len(set(indices)) != len(indices):
        raise ValueError(f"{kind} needs {KINDS[kind]} distinct atom indices")
    return tuple(atom_index(i, count) for i in indices)


def coordinate(positions, kind, indices):
    """Bond length in Angstrom, or angle/dihedral in radians (dihedral in (-pi, pi])."""
    positions = np.asarray(positions, dtype=float)
    indices = _check(kind, indices, len(positions))
    if kind == "bond":
        i, j = indices
        return float(np.linalg.norm(positions[j] - positions[i]))
    if kind == "angle":
        i, j, k = indices
        u = positions[i] - positions[j]
        v = positions[k] - positions[j]
        u /= np.linalg.norm(u)
        v /= np.linalg.norm(v)
        return float(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0)))
    i, j, k, l = indices  # noqa: E741 - i-j-k-l is the dihedral's own notation
    b1 = positions[j] - positions[i]
    b2 = positions[k] - positions[j]
    b3 = positions[l] - positions[k]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    axis = b2 / np.linalg.norm(b2)
    return float(np.arctan2(np.dot(np.cross(n1, axis), n2), np.dot(n1, n2)))


def named_values(positions, specs):
    """Named internal coordinates evaluated on one frame, for the path record.

    Reaction progress has to be readable off a trajectory by name. The trial record already
    stores every pair distance, so the information was there, but only as an unlabelled
    upper-triangle vector -- which is why "why did this frame not fall into B" could not be
    answered from the record. A named list is supplied by the caller because which
    coordinate is the reaction coordinate is a statement about the chemistry, not something
    to infer.

    `bond_difference` takes four indices and returns r(i,j) - r(k,l). A transfer coordinate
    is that difference: it is negative on one side of the acceptor, positive on the other,
    and zero where the transferring atom is shared, so its sign is the progress variable.

    A spec that names a nonexistent atom is an error at definition time, not something to
    swallow per frame, so no exception is caught here.
    """
    positions = np.asarray(positions, dtype=float)
    values = {}
    for spec in specs or ():
        name = spec["name"]
        kind = spec["kind"]
        # bond_difference is validated here rather than by _check: its four indices share
        # the transferring atom by construction, so only each bond's own pair has to be
        # distinct. The two coordinate() calls below re-check each pair.
        indices = tuple(atom_index(i, len(positions)) for i in spec["indices"])
        if kind == "bond_difference":
            if len(indices) != 4:
                raise ValueError("bond_difference needs four atom indices: i j k l")
            values[name] = coordinate(positions, "bond", indices[:2]) - coordinate(
                positions, "bond", indices[2:]
            )
        else:
            values[name] = coordinate(positions, kind, indices)
    return values


def _wrap(delta):
    return (delta + np.pi) % (2 * np.pi) - np.pi


def _rigid_basis(positions):
    """Orthonormal basis of the rigid-body subspace in the plain Cartesian metric.

    In that metric <g, translation_a> = 0 is exactly sum_i g_i = 0 and
    <g, rotation_a> = 0 is exactly sum_i x_i x g_i = 0. Removing this subspace from a
    finite-difference gradient therefore makes both conservation laws hold to machine
    precision instead of to differencing precision, and discards differencing noise
    that would otherwise show up as spurious drift or spin.
    """
    count = len(positions)
    centred = positions - positions.mean(axis=0)
    basis = []
    for axis in range(3):
        translation = np.zeros((count, 3))
        translation[:, axis] = 1.0
        for candidate in (translation.ravel(), np.cross(np.eye(3)[axis], centred).ravel()):
            for kept in basis:
                candidate = candidate - np.dot(candidate, kept) * kept
            norm = np.linalg.norm(candidate)
            if norm > 1e-8:
                basis.append(candidate / norm)
    return basis


def gradient(positions, kind, indices):
    """d(coordinate)/dx as an (N, 3) array, with rigid-body components removed."""
    positions = np.asarray(positions, dtype=float)
    indices = _check(kind, indices, len(positions))
    result = np.zeros_like(positions)
    wrap = _wrap if PERIODIC[kind] else (lambda d: d)
    for atom in indices:
        for axis in range(3):
            shifted = positions.copy()
            shifted[atom, axis] += _STEP
            plus = coordinate(shifted, kind, indices)
            shifted[atom, axis] -= 2 * _STEP
            minus = coordinate(shifted, kind, indices)
            result[atom, axis] = wrap(plus - minus) / (2 * _STEP)
    if not np.isfinite(result).all() or np.linalg.norm(result) < 1e-9:
        raise ValueError(f"Degenerate {kind} geometry for indices {indices}")
    flat = result.ravel()
    for row in _rigid_basis(positions):
        flat = flat - np.dot(flat, row) * row
    result = flat.reshape(result.shape)
    if np.linalg.norm(result) < 1e-9:
        raise ValueError(f"{kind} {indices} has no internal component at this geometry")
    return result


def displace(atoms, kind, indices, delta, in_place=True, tolerance=1e-9, max_iterations=None):
    """Move by delta along one internal coordinate, preserving the center of mass.

    Each increment is capped at MAX_STEP so the path follows the arc of the coordinate
    instead of its chord; without that cap a large torsion displacement reaches the
    requested angle by tearing a bond off. Iteration also removes the first-order error
    of a single step, so the requested amplitude is the delivered amplitude, which the
    reported amplitude windows depend on. Iteration is free in the resource that
    matters: it re-evaluates the geometry, never the potential.
    """
    indices = _check(kind, indices, len(atoms))
    masses = atoms.get_masses()
    target = atoms if in_place else atoms.copy()
    start = coordinate(target.positions, kind, indices)
    # Guard the requested destination, not the step: a bond length is a norm and can
    # never come out negative, so checking the result after the fact catches nothing.
    if kind == "bond" and start + delta <= 0:
        raise ValueError(
            f"Bond displacement targets a nonpositive length "
            f"({start:.3f} + {delta:.3f} Angstrom)"
        )
    if kind == "angle" and not 0 < start + delta < np.pi:
        raise ValueError(f"Angle displacement targets {start + delta:.3f} rad, outside (0, pi)")
    cap = MAX_STEP[kind]
    if max_iterations is None:
        max_iterations = int(abs(delta) / cap) + 24
    remaining = float(delta)
    for _ in range(max_iterations):
        if abs(remaining) < tolerance:
            break
        g = gradient(target.positions, kind, indices)
        direction = g / masses[:, None]
        slope = float(np.sum(g * direction))
        if abs(slope) < 1e-12:
            raise ValueError("Internal coordinate does not respond to its own gradient")
        increment = float(np.clip(remaining, -cap, cap))
        target.set_positions(target.positions + direction * (increment / slope))
        achieved = coordinate(target.positions, kind, indices) - start
        remaining = float(delta) - (_wrap(achieved) if PERIODIC[kind] else achieved)
    total = coordinate(target.positions, kind, indices) - start
    return target, float(_wrap(total) if PERIODIC[kind] else total)


def rotate_fragment(positions, axis_a, axis_b, moving, angle, masses=None):
    """Rigidly rotate a set of atoms about the axis through two atoms.

    Exact where it applies. Both axis atoms lie on the axis and do not move, so every
    bond length and every bond angle is preserved to machine precision: an angle whose
    three atoms are all stationary is untouched, and an angle j-k-l with j and k on the
    axis is untouched because rotation about the axis preserves l's angle to it. Only
    torsions about that bond change. The centre of mass does move, so it is restored by
    an overall translation, which changes no internal coordinate.

    This is the right realization for a torsion about a bridge bond. A ring torsion has
    no rigid fragment to rotate and must use the continuation scheme in displace().
    """
    positions = np.array(positions, dtype=float)
    count = len(positions)
    axis_a = atom_index(axis_a, count)
    axis_b = atom_index(axis_b, count)
    moving = np.asarray(sorted({atom_index(i, count) for i in moving}), dtype=int)
    if axis_a in moving or axis_b in moving:
        raise ValueError("Axis atoms must not be part of the rotating fragment")
    origin = positions[axis_a]
    axis = positions[axis_b] - origin
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        raise ValueError("Rotation axis is undefined for coincident atoms")
    axis = axis / norm
    before = None
    if masses is not None:
        masses = np.asarray(masses, dtype=float)
        before = np.average(positions, axis=0, weights=masses)
    offsets = positions[moving] - origin
    along = np.outer(offsets @ axis, axis)
    across = offsets - along
    perpendicular = np.cross(axis, across)
    positions[moving] = origin + along + across * np.cos(angle) + perpendicular * np.sin(angle)
    if before is not None:
        positions -= np.average(positions, axis=0, weights=masses) - before
    return positions


def kick(atoms, kind, indices, energy_eV, sign=1):
    """Inject exactly energy_eV of kinetic energy along one internal coordinate.

    Solves A q^2 + B q = energy_eV for q >= 0, so the injected energy is exact even at
    nonzero initial momentum; the cross term B must not be dropped. Total linear and
    angular momentum are unchanged because the direction is an internal-coordinate
    gradient.
    """
    if energy_eV <= 0 or not np.isfinite(energy_eV) or sign not in (-1, 1):
        raise ValueError("Kick requires positive finite energy and sign +/-1")
    g = sign * gradient(atoms.positions, kind, indices)
    masses = atoms.get_masses()
    momenta = atoms.get_momenta()
    a = 0.5 * float(np.sum(g * g / masses[:, None]))
    b = float(np.sum(momenta * g / masses[:, None]))
    radical = np.sqrt(b * b + 4 * a * energy_eV)
    q = 2 * energy_eV / (radical + b) if b >= 0 else (radical - b) / (2 * a)
    atoms.set_momenta(momenta + q * g)
    return float(q)
