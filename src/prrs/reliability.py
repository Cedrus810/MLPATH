"""Fail-closed guards executed for every force evaluation, including FIRE."""
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.data import covalent_radii


class GateRejected(RuntimeError):
    def __init__(self, code, message, atoms):
        super().__init__(message)
        self.code = code
        self.atoms = atoms.copy()


class GuardedCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, physical, config, reference_energy=None):
        super().__init__()
        self.physical = physical
        self.config = config
        self.reference_energy = reference_energy

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        cfg = self.config

        def reject(code, message):
            raise GateRejected(code, message, atoms)

        if not np.isfinite(atoms.positions).all():
            reject("nonfinite", "Non-finite coordinates")
        radii = covalent_radii[atoms.numbers]
        limits = cfg.min_distance_ratio * (radii[:, None] + radii[None, :])
        distances = atoms.get_all_distances()
        np.fill_diagonal(distances, np.inf)
        if np.any(distances < limits):
            reject("overlap", "An interatomic distance is below the overlap guard")
        energy = float(self.physical.get_potential_energy(atoms))
        forces = np.asarray(self.physical.get_forces(atoms))
        if forces.shape != (len(atoms), 3) or not np.isfinite(forces).all() or not np.isfinite(energy):
            reject("nonfinite", "Calculator returned invalid energy or forces")
        if np.max(np.linalg.norm(forces, axis=1)) > cfg.max_force_eV_A:
            reject("force_limit", "Physical force exceeds configured limit")
        if self.reference_energy is not None:
            if (energy - self.reference_energy) / len(atoms) > cfg.max_energy_rise_eV_atom:
                reject("energy_limit", "Physical energy rise exceeds configured limit")
        uncertainty = self.physical.results.get("force_uncertainty_eV_A")
        if uncertainty is not None:
            uncertainty = float(uncertainty)
            if not np.isfinite(uncertainty) or uncertainty < 0:
                reject("invalid_uncertainty", "Invalid uncertainty value")
        if cfg.max_uncertainty_eV_A is not None:
            if uncertainty is None:
                reject("uncertainty_unavailable", "Uncertainty required but calculator supplies none")
            if uncertainty > cfg.max_uncertainty_eV_A:
                reject("ood", "Force disagreement exceeds configured limit")
        self.results = {"energy": energy, "forces": forces.copy(),
                        "force_uncertainty_eV_A": uncertainty}


class PairPulse(Calculator):
    """Conservative central bias U=-sign*alpha*(r-r0), held constant during pulse."""
    implemented_properties = ["energy", "forces"]

    def __init__(self, physical, pair, sign, amplitude, reference_distance):
        super().__init__()
        self.physical = physical
        self.pair = pair
        self.sign = sign
        self.amplitude = amplitude
        self.reference_distance = reference_distance

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        from .perturbations import pair_axis
        super().calculate(atoms, properties, system_changes)
        energy = self.physical.get_potential_energy(atoms)
        forces = self.physical.get_forces(atoms).copy()
        axis, distance = pair_axis(atoms, self.pair)
        i, j = self.pair
        force = self.sign * self.amplitude * axis
        forces[i] -= force
        forces[j] += force
        bias = -self.sign * self.amplitude * (distance - self.reference_distance)
        self.results = {"energy": energy + bias, "forces": forces, "bias_energy": bias}

