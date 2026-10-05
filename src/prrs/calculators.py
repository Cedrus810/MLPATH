"""Physical potentials and calculator factories. No model is downloaded implicitly."""

from importlib import import_module
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.data import atomic_numbers, covalent_radii


class DoubleWell(Calculator):
    """Synthetic C-C radial double well; software benchmark, NOT carbon chemistry.

    Minima r=1.2 and 2.4 Angstrom, barrier 0.4 eV at r=1.8 Angstrom.
    """

    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        if len(atoms) != 2:
            raise ValueError("DoubleWell requires exactly two atoms")
        vector = atoms.positions[1] - atoms.positions[0]
        r = np.linalg.norm(vector)
        if r < 1e-12:
            raise ValueError("Coincident atoms")
        x = (r - 1.8) / 0.6
        energy = 0.4 * (x * x - 1) ** 2
        derivative = 4 * 0.4 * x * (x * x - 1) / 0.6
        force = derivative * vector / r
        self.results = {"energy": float(energy), "forces": np.array([force, -force])}


def double_well_factory():
    return DoubleWell()


def demo_atoms():
    # Only the reactant is supplied to the search.
    return Atoms("C2", positions=[[-0.6, 0, 0], [0.6, 0, 0]])


# ---- connect-mode testbed (PLAN_2026-09-28.md item T0) ------------------------------

RELAY_HEAVY = ("F", "C", "O", "Cl")
RELAY_SPACING_A = 3.0
RELAY_BOND_SCALE = 1.2


def _relay_well_radius(symbol):
    """Well radius 0.8 x the H-X bond cutoff: bonded at the minimum, never overlapping."""
    return 0.8 * RELAY_BOND_SCALE * (covalent_radii[1] + covalent_radii[atomic_numbers[symbol]])


def _spring(x, grad, i, j, k, r0):
    v = x[j] - x[i]
    r = np.linalg.norm(v)
    g = 2 * k * (r - r0) * v / r
    grad[j] += g
    grad[i] -= g
    return k * (r - r0) ** 2


def _angle(x, grad, a, b, c, k, cos0=None):
    """k (1 + cos abc) with cos0 None (straight, regular at pi), else k (cos - cos0)^2."""
    u, v = x[a] - x[b], x[c] - x[b]
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    cos = float(u @ v) / (nu * nv)
    if cos0 is None:
        energy = k * (1 + cos)
        w = k
    else:
        energy = k * (cos - cos0) ** 2
        w = 2 * k * (cos - cos0)
    da = w * (v / (nu * nv) - cos * u / nu**2)
    dc = w * (u / (nu * nv) - cos * v / nv**2)
    grad[a] += da
    grad[c] += dc
    grad[b] -= da + dc
    return energy


def _dihedral(x, grad, i, j, k, m, kd):
    """kd (1 + cos phi): the dihedral held at 180 degrees.

    Never kd (cos - cos0)^2 at 0 or 180: cos is flat there, so that form has zero
    curvature and would leave the very mode it is meant to hold."""
    a, b, c = x[j] - x[i], x[k] - x[j], x[m] - x[k]
    p, q = np.cross(a, b), np.cross(b, c)
    npp, nq = np.linalg.norm(p), np.linalg.norm(q)
    cos = float(p @ q) / (npp * nq)
    energy = kd * (1 + cos)
    w = kd
    gp = w * (q / (npp * nq) - cos * p / npp**2)
    gq = w * (p / (npp * nq) - cos * q / nq**2)
    da = np.cross(b, gp)
    db = np.cross(gp, a) + np.cross(c, gq)
    dc = np.cross(gq, b)
    grad[i] -= da
    grad[j] += da - db
    grad[k] += db - dc
    grad[m] += dc
    return energy


def _well(x, grad, h, j, depth, radius, sigma):
    """A Gaussian shell of radius `radius` around atom j, felt by atom h."""
    v = x[h] - x[j]
    r = np.linalg.norm(v)
    offset = r - radius
    e = -depth * np.exp(-(offset**2) / (2 * sigma**2))
    g = -e * offset / sigma**2 * v / r
    grad[h] += g
    grad[j] -= g
    return e


class ProtonRelay(Calculator):
    """Synthetic proton wire (Grotthuss chain); software testbed, NOT F/C/O/Cl chemistry.

    Four heavy atoms F-C-O-Cl in a planar trans zigzag, 3.0 A apart and therefore never
    bonded to one another; gap k holds proton H_k (atom 4 + k). Each proton sits in one
    of two wells, on its left atom or its right one, so the 2^3 placements are the
    minima and a hop changes exactly one bond. The reactant has every proton on the
    left, the product every proton on the right: d to the product falls 6 -> 4 -> 2 ->
    0 along any single-hop path, the multi-step shape connect mode needs. Four distinct
    elements make the eight placements eight distinct chemical keys (a chain of one
    element would merge mirror placements), and F/C/O/Cl are all in chemistry.VALENCE.

    What the draft testbed (PLAN_2026-09-28.py) got wrong, found in review, and what
    this one is built against:
      * wells at fixed points in space plus tethers -- not translation invariant (the
        draft's own note records the mainline Hessian provenance gate refusing its
        minima). Every term here is a function of distances and angles only;
      * wells centred ON the heavy nuclei -- H-C distance 0, refused by the overlap
        guard (min_distance_ratio). Wells here are shells at 0.8 x the bond cutoff;
      * heavy atoms 1.5 A apart -- a proton half way was bonded to both, so the
        product-side bond formed inside the reactant well. At 3.0 A the new bond forms
        only past the barrier top (tests/test_relay_testbed.py holds that per gap);
      * decoy and bridge atoms with no force on them, drifting in dynamics while their
        wells stayed put. Every atom here is held by a term.
    Measured on the prototype: a straight chain failed confirm_minimum on 2 of 8 minima
    (trivial_mode_floor_above_tolerance -- the axial rotation of a linear molecule); the
    zigzag passes all eight. The held dihedral is NOT what confirm_minimum needs -- it
    passes without it, tolerating the flat torsion -- it is what makes the only zero
    modes the six rigid ones (tests/test_relay_testbed.py counts them).
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        depths_eV=(1.10, 1.15, 1.20, 1.25),
        sigma_A=0.25,
        k_spring=20.0,
        k_bend=5.0,
        k_axis=3.0,
        k_torsion=2.0,
        backbone_angle_deg=150.0,
    ):
        super().__init__()
        self.depths = tuple(depths_eV)
        self.sigma = sigma_A
        self.k_spring, self.k_bend, self.k_axis, self.k_torsion = (
            k_spring,
            k_bend,
            k_axis,
            k_torsion,
        )
        self.cos_backbone = float(np.cos(np.deg2rad(backbone_angle_deg)))
        self.radii = tuple(_relay_well_radius(s) for s in RELAY_HEAVY)

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        n = len(RELAY_HEAVY)
        if len(atoms) != 2 * n - 1:
            raise ValueError(f"ProtonRelay requires {2 * n - 1} atoms (4 heavy + 3 H)")
        grad = np.zeros_like(atoms.positions)
        energy = self._chain(atoms.positions, grad)
        self.results = {"energy": float(energy), "forces": -grad}

    def _chain(self, x, grad):
        """The wire's own terms on atoms 0..6; accumulates into grad, returns energy."""
        n = len(RELAY_HEAVY)
        energy = 0.0
        for j in range(n - 1):
            energy += _spring(x, grad, j, j + 1, self.k_spring, RELAY_SPACING_A)
        for j in range(1, n - 1):
            energy += _angle(x, grad, j - 1, j, j + 1, self.k_bend, self.cos_backbone)
        energy += _dihedral(x, grad, 0, 1, 2, 3, self.k_torsion)
        for k in range(n - 1):
            h = n + k
            energy += _angle(x, grad, k, h, k + 1, self.k_axis)  # stays on its own gap's axis
            for j in range(n):
                energy += _well(x, grad, h, j, self.depths[j], self.radii[j], self.sigma)
        return energy


def proton_relay_factory():
    return ProtonRelay()


def proton_relay_atoms(sides=(0, 0, 0)):
    """sides[k] = 0 puts proton k on its left atom, 1 on its right one; (0,0,0) is R."""
    n = len(RELAY_HEAVY)
    if len(sides) != n - 1 or set(sides) - {0, 1}:
        raise ValueError(f"sides must be {n - 1} values in (0, 1)")
    tilt = (np.pi - np.deg2rad(150.0)) / 2
    heavy = [np.zeros(3)]
    for j in range(n - 1):
        sign = 1 if j % 2 == 0 else -1
        step = np.array([np.cos(tilt), sign * np.sin(tilt), 0.0])
        heavy.append(heavy[-1] + RELAY_SPACING_A * step)
    positions = [h.copy() for h in heavy]
    for k, side in enumerate(sides):
        axis = (heavy[k + 1] - heavy[k]) / RELAY_SPACING_A
        if side == 0:
            positions.append(heavy[k] + _relay_well_radius(RELAY_HEAVY[k]) * axis)
        else:
            positions.append(heavy[k + 1] - _relay_well_radius(RELAY_HEAVY[k + 1]) * axis)
    return Atoms(list(RELAY_HEAVY) + ["H"] * (n - 1), positions=positions)


RELAY_GATE_ATOM = "Br"
RELAY_GATE_SPACING_A = 3.6


class ProtonRelayGated(ProtonRelay):
    """The wire plus a gatekeeper: T-detour of PLAN_2026-09-28.md T0.

    Atoms 0..6 are the wire, unchanged. Atom 7 is Br, 3.6 A from Cl (atom 3), continuing
    the planar trans zigzag; atom 8 is a fourth proton, the gatekeeper, with two wells:
    on Cl (home, deeper) and on Br. The one coupling term,

        A exp(-(r(H5,C) - r(H5,O))^2 / 2 s^2) exp(-(r(H8,Cl) - r_Cl)^2 / 2 s_g^2),

    raises proton 5's transit through the middle of the C-O gap by A while the gatekeeper
    is home, and not at all once it is out on Br. Reactant and product both have it home,
    so the cheap way from R to P sends it out (d up by 2), lets proton 5 through, and
    brings it back: the intermediate that must RAISE the distance to the product before
    anything can lower it, which connect mode's upper layer has to survive by not pruning
    (plan section 1 item 4).

    Why the home is Cl and not C, measured on the first version: with the gatekeeper home
    on C, it and proton 5 are two hydrogens on one carbon in the reactant, so Aut(G_R)
    swaps them, the minimum over it reads "an H leaves C for Br" as proton 5's move as
    well as the gatekeeper's, and the excursion came out d-neutral (6, 4, 4, 4, 2, 0).
    Cl carries no proton in the reactant, so every hydrogen there sits on a different
    element and the group is trivial.

    Distances and angles only, like the wire; the gatekeeper stays on the Cl-Br axis, Br
    is held by a spring, the backbone angle and a trans dihedral.
    """

    def __init__(self, gate_eV=1.5, gate_width_A=0.25, gatekeeper_width_A=0.3, **wire):
        super().__init__(**wire)
        self.gate, self.gate_width, self.gatekeeper_width = (
            gate_eV,
            gate_width_A,
            gatekeeper_width_A,
        )
        self.gate_radius = _relay_well_radius(RELAY_GATE_ATOM)
        self.home_depth, self.away_depth = 1.20, 1.10

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        Calculator.calculate(self, atoms, properties, system_changes)
        if len(atoms) != 9:
            raise ValueError("ProtonRelayGated requires 9 atoms (wire + Br + gatekeeper)")
        x = atoms.positions
        grad = np.zeros_like(x)
        energy = self._chain(x, grad)
        c, o, cl, h5, br, h8 = 1, 2, 3, 5, 7, 8
        energy += _spring(x, grad, cl, br, self.k_spring, RELAY_GATE_SPACING_A)
        energy += _angle(x, grad, o, cl, br, self.k_bend, self.cos_backbone)
        energy += _dihedral(x, grad, c, o, cl, br, self.k_torsion)
        energy += _angle(x, grad, cl, h8, br, self.k_axis)
        energy += _well(x, grad, h8, cl, self.home_depth, self.radii[cl], self.sigma)
        energy += _well(x, grad, h8, br, self.away_depth, self.gate_radius, self.sigma)

        # the gate: proton 5 mid-transit times the gatekeeper at home
        v5c, v5o, v8 = x[h5] - x[c], x[h5] - x[o], x[h8] - x[cl]
        r5c, r5o, r8 = (np.linalg.norm(v) for v in (v5c, v5o, v8))
        u5c, u5o, u8 = v5c / r5c, v5o / r5o, v8 / r8
        delta, home = r5c - r5o, r8 - self.radii[cl]
        transit = np.exp(-(delta**2) / (2 * self.gate_width**2))
        occupied = np.exp(-(home**2) / (2 * self.gatekeeper_width**2))
        gate = self.gate * transit * occupied
        energy += gate
        d_delta = -gate * delta / self.gate_width**2  # dE/d(delta)
        d_home = -gate * home / self.gatekeeper_width**2  # dE/d(r(H8, Cl))
        grad[h5] += d_delta * (u5c - u5o)
        grad[c] -= d_delta * u5c
        grad[o] += d_delta * u5o
        grad[h8] += d_home * u8
        grad[cl] -= d_home * u8
        self.results = {"energy": float(energy), "forces": -grad}


def proton_relay_gated_factory():
    return ProtonRelayGated()


def proton_relay_gated_atoms(sides=(0, 0, 0), gatekeeper_home=True):
    """The wire's placement plus Br and the gatekeeper, home on Cl or out on Br."""
    atoms = proton_relay_atoms(sides)
    x = atoms.positions
    # the zigzag's next link: mirror of the O -> Cl step across the Cl-O line's normal
    step = (x[3] - x[2]) / RELAY_SPACING_A
    step = np.array([step[0], -step[1], 0.0])
    br = x[3] + RELAY_GATE_SPACING_A * step
    if gatekeeper_home:
        h8 = x[3] + _relay_well_radius(RELAY_HEAVY[3]) * step
    else:
        h8 = br - _relay_well_radius(RELAY_GATE_ATOM) * step
    return Atoms(
        atoms.get_chemical_symbols() + [RELAY_GATE_ATOM, "H"],
        positions=np.vstack([x, br, h8]),
    )


class Committee(Calculator):
    """Mean PES plus maximum per-atom RMS vector force disagreement in eV/A."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, members):
        super().__init__()
        if len(members) < 2:
            raise ValueError("An uncertainty committee requires at least two models")
        if len({id(m) for m in members}) != len(members):
            raise ValueError("Committee members must be separate calculators")
        self.members = members

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        energies, forces = [], []
        for member in self.members:
            energies.append(member.get_potential_energy(atoms))
            forces.append(member.get_forces(atoms))
        force_array = np.asarray(forces)
        mean = force_array.mean(axis=0)
        disagreement = np.sqrt(np.mean(np.sum((force_array - mean) ** 2, axis=2), axis=0))
        self.results = {
            "energy": float(np.mean(energies)),
            "forces": mean,
            "force_uncertainty_eV_A": float(np.max(disagreement)),
        }


def load_factory(spec, kwargs=None):
    """Resolve an explicitly trusted local Python factory module:function."""
    if ":" not in spec:
        raise ValueError("Calculator factory must be module:function")
    module, name = spec.split(":", 1)
    function = getattr(import_module(module), name)
    options = dict(kwargs or {})

    def create():
        calc = function(**options)
        if not isinstance(calc, Calculator):
            raise TypeError("Factory must return an ASE Calculator")
        return calc

    return create


def mace_factory(model_paths, device="cpu", default_dtype="float64"):
    """Local checkpoints only; multiple checkpoints enable committee uncertainty."""
    paths = [model_paths] if isinstance(model_paths, str) else list(model_paths)
    if not paths or any(not Path(p).is_file() for p in paths):
        raise ValueError("Supply existing local MACE checkpoint paths")
    try:
        from mace.calculators import MACECalculator
    except ImportError as exc:
        raise ImportError("Install mlpath-prrs[mace] to load MACE checkpoints") from exc
    members = [
        MACECalculator(model_paths=str(p), device=device, default_dtype=default_dtype)
        for p in paths
    ]
    return members[0] if len(members) == 1 else Committee(members)
