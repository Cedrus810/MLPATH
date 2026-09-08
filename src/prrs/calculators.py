"""Physical potentials and calculator factories. No model is downloaded implicitly."""

from importlib import import_module
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes


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
