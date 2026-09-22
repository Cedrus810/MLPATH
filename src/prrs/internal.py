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


def translate_fragments(positions, i, j, side_i, side_j, delta, masses):
    """Rigidly translate the two sides of a bond along its axis, by delta in total.

    The bond analogue of rotate_fragment, and exact for the same reason: each side moves
    as a rigid body, so every bond length, bond angle and torsion inside a side is
    preserved to machine precision, and the target distance changes by exactly delta
    because both translations are along the bond axis. The two displacements are split in
    inverse proportion to the sides' masses, which keeps the centre of mass fixed without
    a correcting translation afterwards.

    displace() moves along `gradient / masses`, and a bond gradient is nonzero on its two
    atoms only, so the rest of each side stays where it was and the bonded atom walks out
    of its own substituents: pulling C-Cl by 0.600 Angstrom on the b04 input shortened
    every C-H by 0.0628 Angstrom, which the collateral gate refused -- correctly, since
    the damage was real, and it punched the hole that left the amplitude ladder without a
    bracket. Atoms not on either side (a spectator fragment) are left alone.

    A bond inside a ring has no rigid split; `_bond_sides` returns None there and delivery
    falls back to the continuation scheme in displace().
    """
    positions = np.array(positions, dtype=float)
    count = len(positions)
    i = atom_index(i, count)
    j = atom_index(j, count)
    side_i = np.asarray(sorted({atom_index(a, count) for a in side_i}), dtype=int)
    side_j = np.asarray(sorted({atom_index(a, count) for a in side_j}), dtype=int)
    if i not in side_i or j not in side_j:
        raise ValueError("Each side of a bond must contain its own bond atom")
    if set(side_i.tolist()) & set(side_j.tolist()):
        raise ValueError("The two sides of a bond must not overlap")
    axis = positions[j] - positions[i]
    length = float(np.linalg.norm(axis))
    if length < 1e-12:
        raise ValueError("Bond axis is undefined for coincident atoms")
    # Guard the destination, as displace() does: a length is a norm and never comes out
    # negative, so a check after the fact catches nothing.
    if length + delta <= 0:
        raise ValueError(
            f"Bond displacement targets a nonpositive length "
            f"({length:.3f} + {delta:.3f} Angstrom)"
        )
    axis = axis / length
    masses = np.asarray(masses, dtype=float)
    mass_i = float(masses[side_i].sum())
    mass_j = float(masses[side_j].sum())
    share_j = mass_i / (mass_i + mass_j)
    positions[side_j] += axis * (delta * share_j)
    positions[side_i] -= axis * (delta * (1.0 - share_j))
    return positions


def rotate_about_vertex(positions, i, vertex, k, side_i, side_k, delta, masses=None):
    """Open an angle by rigidly rotating its two side fragments about its vertex.

    The angle analogue of rotate_fragment, exact for the same reason: each side turns as a
    rigid body about an axis through the vertex, so every bond length and bond angle
    inside a side is preserved to machine precision, and the two bonds to the vertex are
    preserved because the axis passes through it. Only the angle itself, and torsions
    across the vertex, change. The sides take delta/2 each -- the symmetric choice, and
    the one measured against the collateral gate -- and the centre of mass is restored
    afterwards by an overall translation, which changes no internal coordinate.

    Without this, `displace` moved an angle along its own gradient over the masses, and an
    angle's gradient is nonzero on its three atoms alone, so the vertex walked out of its
    substituents exactly as a bond's atoms used to before translate_fragments. Measured on
    b05's relaxed source: a 0.35 rad bend at the tertiary carbon stretched C0-C6 by
    0.1325 A, more than twice the 0.05 A gate, and 4 of that run's 15 refusals were this.
    The same four deliveries through this routine reach the requested angle to 1e-16 with
    a worst non-target bond change of 2e-16 A.

    `side_i` and `side_k` contain their own arm atom, exclude the vertex, and must not
    overlap; an angle inside a ring has no such split and falls back to the continuation
    scheme in displace().
    """
    positions = np.array(positions, dtype=float)
    count = len(positions)
    i = atom_index(i, count)
    vertex = atom_index(vertex, count)
    k = atom_index(k, count)
    side_i = np.asarray(sorted({atom_index(a, count) for a in side_i}), dtype=int)
    side_k = np.asarray(sorted({atom_index(b, count) for b in side_k}), dtype=int)
    if vertex in side_i.tolist() or vertex in side_k.tolist():
        raise ValueError("The vertex must not be part of either side")
    if set(side_i.tolist()) & set(side_k.tolist()):
        raise ValueError("The two sides of an angle must not overlap")
    if i not in side_i.tolist() or k not in side_k.tolist():
        raise ValueError("Each side must contain its own arm atom")
    # Guard the requested destination, as displace() does. The axis check below only
    # catches an angle that is ALREADY linear; a 178 degree angle passes it and then folds
    # through pi, so arccos reads back 2*pi - (theta + delta), target_error comes out near
    # 0.63 rad for a 0.35 rad bend, and the amplitude ladder still records the requested
    # amplitude. That is a wrong number in the record, not a refusal.
    start = coordinate(positions, "angle", (i, vertex, k))
    if not 0 < start + delta < np.pi:
        raise ValueError(f"Angle displacement targets {start + delta:.3f} rad, outside (0, pi)")
    # The axis is the normal of the plane the ANGLE lies in -- defined by its own two
    # arms, never by a member of a side. Taking it from the lowest-numbered atom of each
    # side instead makes the rotation plane depend on how the atoms happen to be
    # numbered, which the relabelling test in tests/test_internal.py caught immediately.
    axis = np.cross(positions[i] - positions[vertex], positions[k] - positions[vertex])
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        raise ValueError("Angle rotation axis is undefined for a linear or degenerate angle")
    axis = axis / norm
    before = None
    if masses is not None:
        masses = np.asarray(masses, dtype=float)
        before = np.average(positions, axis=0, weights=masses)

    def turn(members, angle):
        offsets = positions[members] - positions[vertex]
        along = np.outer(offsets @ axis, axis)
        across = offsets - along
        perpendicular = np.cross(axis, across)
        positions[members] = (
            positions[vertex] + along + across * np.cos(angle) + perpendicular * np.sin(angle)
        )

    turn(side_i, -delta / 2.0)
    turn(side_k, +delta / 2.0)
    if before is not None:
        positions -= np.average(positions, axis=0, weights=masses) - before
    return positions


def rotate_rigid(positions, members, axis, angle, masses=None):
    """Rotate one fragment about an axis through its own centre of mass.

    The third exact realization, after `rotate_fragment` (axis through two atoms) and
    `rotate_about_vertex` (axis normal to an angle). This one needs no atom on the axis,
    which is the whole point: the coordinate it moves is a fragment's orientation relative
    to the rest of the system, and between two fragments there is no bond to hang an axis
    on. Every intra-fragment bond length, angle and torsion is preserved to machine
    precision because the fragment turns as a rigid body; nothing outside `members` moves.

    That gap is not a budget: `polish_soft_modes` iterates `rotatable_torsions`, which
    enumerates bonds of the graph, and a water molecule beside an alkene shares no bond
    with it. So its orientation was outside the polish operator's DOMAIN, and b05's three
    elimination products stopped on order-1 stationary points whose unstable mode was
    99.9% a rigid rotation of the water -- a coordinate nothing in the polish could move.

    The overall centre of mass is restored afterwards by a translation, which changes no
    internal coordinate.
    """
    positions = np.array(positions, dtype=float)
    members = np.asarray(sorted({atom_index(a, len(positions)) for a in members}), dtype=int)
    if not len(members):
        raise ValueError("A rigid rotation needs at least one atom to move")
    axis = np.asarray(axis, dtype=float)
    norm = float(np.linalg.norm(axis))
    if norm < 1e-12:
        raise ValueError("Rigid rotation axis is undefined")
    axis = axis / norm
    weights = np.ones(len(members)) if masses is None else np.asarray(masses, float)[members]
    origin = np.average(positions[members], axis=0, weights=weights)
    before = None
    if masses is not None:
        masses = np.asarray(masses, dtype=float)
        before = np.average(positions, axis=0, weights=masses)

    offsets = positions[members] - origin
    along = np.outer(offsets @ axis, axis)
    across = offsets - along
    positions[members] = (
        origin + along + across * np.cos(angle) + np.cross(axis, across) * np.sin(angle)
    )
    if before is not None:
        positions -= np.average(positions, axis=0, weights=masses) - before
    return positions


def translate_rigid(positions, members, direction, distance, masses=None):
    """Translate one fragment rigidly along a unit direction (B'', PLAN item 5.1).

    The fourth exact relative-motion realization, and the counterpart of
    `rotate_rigid`: between two fragments the three relative translations have no
    internal coordinate to hang on, exactly as their rotations had no bond. Every
    intra-fragment distance, angle and torsion is preserved to machine precision
    because the fragment moves as a rigid body; the total centre of mass is restored
    afterwards by translating everything, which changes no internal coordinate and
    leaves the relative displacement between the fragments as the only physical change.

    `polish_soft_modes` footnote [1] is why this is a separate operator with its own
    thresholds and not a number reused from the rotational set: rotations are measured
    in radians and translations in Angstrom, and one threshold mechanism fed both unit
    systems is the silent kind of wrong. The thresholds themselves come from a
    measurement on the b05 endpoints (fragment_translation_*), not converted from the
    radian set.
    """
    positions = np.array(positions, dtype=float)
    members = np.asarray(sorted({atom_index(a, len(positions)) for a in members}), dtype=int)
    if not len(members):
        raise ValueError("A rigid translation needs at least one atom to move")
    direction = np.asarray(direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm < 1e-12:
        raise ValueError("Rigid translation direction is undefined")
    direction = direction / norm
    before = None
    if masses is not None:
        masses = np.asarray(masses, dtype=float)
        before = np.average(positions, axis=0, weights=masses)
    positions[members] += distance * direction
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
