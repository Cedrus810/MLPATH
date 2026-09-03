"""Verified openmm-ml/MACE deployment primitives for the OpenMM cross-check backend.

Scope: build a Context that evaluates the same physical Hamiltonian as the ASE path,
and nothing more. Trajectory execution belongs to an Engine built on top of this.

Four properties of openmm-ml 1.6 were measured rather than assumed, and each is
enforced here because a silent default would change the physics:

1. The MACE force is an ``openmm.PythonForce`` whose callback rebuilds the MACE
   neighborhood from the current positions on every evaluation. Nothing about the
   bond graph is frozen, so one Context stays valid across bond breaking and forming
   and there is no reason to rebuild it when the graph changes.
2. ``createSystem`` adds a ``CMMotionRemover`` unless told otherwise. That force
   removes center-of-mass motion and would quietly invalidate the linear-momentum
   conservation checks the perturbation protocol depends on, so it is refused here.
3. The local-``modelPath`` branch of openmm-ml's MACE implementation loads with
   ``map_location`` only and never calls ``.to(device)`` the way its named-model
   branch does; constants inside the scripted submodules stay on the host and every
   CUDA evaluation fails. A narrow, reversible patch is applied during createSystem.
4. Platform selection never falls back. A run that silently drops from CUDA to CPU is
   a different experiment, and the manifest must not have to guess which one ran.
"""
from contextlib import contextmanager
from pathlib import Path
import numpy as np

# openmm-ml reports energies in kJ/mol and forces in kJ/mol/nm.
EV_PER_KJ_PER_MOL = 1.0 / 96.4853
NM_PER_ANGSTROM = 0.1
PRECISIONS = ("single", "double")


def build_topology(numbers):
    """One residue of bare atoms, with no bonds declared.

    Declaring bonds would hand the model a fixed adjacency; MACE derives its own
    neighborhood from geometry, which is exactly what reaction discovery requires.
    """
    from openmm import app
    numbers = np.asarray(numbers, dtype=int)
    if numbers.ndim != 1 or numbers.size < 2 or np.any(numbers <= 0):
        raise ValueError("numbers must be at least two positive atomic numbers")
    topology = app.Topology()
    residue = topology.addResidue("MOL", topology.addChain())
    for number in numbers:
        element = app.Element.getByAtomicNumber(int(number))
        topology.addAtom(element.symbol, element, residue)
    return topology


@contextmanager
def _local_model_moved_to_device():
    """Make openmm-ml's local-modelPath MACE branch honour its device argument.

    Remove this once openmm-ml moves the locally loaded model itself; the paired
    regression test fails loudly when the upstream behaviour changes.
    """
    import torch
    original = torch.load

    def load_and_move(f, *args, **kwargs):
        model = original(f, *args, **kwargs)
        device = kwargs.get("map_location")
        return model.to(device) if device is not None and hasattr(model, "to") else model

    torch.load = load_and_move
    try:
        yield
    finally:
        torch.load = original


def create_mace_context(numbers, model_path, device="cuda", precision="double",
                        timestep_fs=0.25, platform_name=None, charge=0, multiplicity=1):
    """Build (context, system, metadata) for a fixed atom set on one explicit platform.

    charge and multiplicity are baked into the System by openmm-ml, so changing either
    requires a new System; they are returned in the metadata to keep that on record.
    """
    import openmm
    from openmm import unit
    from openmmml import MLPotential

    model_path = Path(model_path)
    if not model_path.is_file():
        raise FileNotFoundError(f"MACE checkpoint not found: {model_path}")
    if precision not in PRECISIONS:
        raise ValueError(f"precision must be one of {PRECISIONS}")
    if platform_name is None:
        platform_name = "CUDA" if device.startswith("cuda") else "CPU"

    potential = MLPotential("mace", modelPath=str(model_path))
    with _local_model_moved_to_device():
        system = potential.createSystem(
            build_topology(numbers), device=device, precision=precision,
            returnEnergyType="energy", charge=charge, multiplicity=multiplicity,
            removeCMMotion=False)
    force_names = [type(system.getForce(i)).__name__ for i in range(system.getNumForces())]
    if "CMMotionRemover" in force_names:
        raise RuntimeError("CMMotionRemover present; it would break momentum conservation")

    # No try/except: an unavailable platform must fail, not silently become CPU.
    platform = openmm.Platform.getPlatformByName(platform_name)
    integrator = openmm.VerletIntegrator(timestep_fs * unit.femtoseconds)
    context = openmm.Context(system, integrator, platform)
    metadata = {"backend": "openmm-ml/mace", "model_path": str(model_path),
                "model_sha256": _file_sha256(model_path), "device": device,
                "precision": precision, "platform": platform.getName(),
                "charge": charge, "multiplicity": multiplicity,
                "return_energy_type": "energy", "forces": force_names,
                "timestep_fs": timestep_fs,
                "neighborhood": "rebuilt from positions on every evaluation"}
    return context, system, metadata


def evaluate(context, positions_A):
    """Single-point energy and forces in eV and eV/Angstrom."""
    from openmm import unit
    positions_A = np.asarray(positions_A, dtype=float)
    context.setPositions(positions_A * NM_PER_ANGSTROM)
    state = context.getState(getEnergy=True, getForces=True)
    energy = (state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
              * EV_PER_KJ_PER_MOL)
    forces = (state.getForces(asNumpy=True)
              .value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
              * EV_PER_KJ_PER_MOL * NM_PER_ANGSTROM)
    return float(energy), np.asarray(forces, dtype=float)


def _file_sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
