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
   CUDA evaluation fails. The move is applied during createSystem by
   ``torch_guard.moved_to_device``, which flags the single process-wide ``torch.load``
   wrapper rather than installing a second one. This module used to wrap ``torch.load``
   itself and restore it unconditionally on exit, which deleted any other wrapper
   installed while the context was open. See ``torch_guard`` for why that is now
   structurally impossible rather than merely avoided.
4. Platform selection never falls back. A run that silently drops from CUDA to CPU is
   a different experiment, and the manifest must not have to guess which one ran.
"""

from pathlib import Path
import numpy as np

from . import torch_guard

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


def create_mace_context(
    numbers,
    model_path,
    device="cuda",
    precision="double",
    timestep_fs=0.25,
    platform_name=None,
    charge=0,
    multiplicity=1,
):
    """Build (context, system, metadata) for a fixed atom set on one explicit platform.

    charge and multiplicity are baked into the System by openmm-ml, so changing either
    requires a new System; they are returned in the metadata to keep that on record.

    Notes
    -----
    [1] Before openmmml, because importing it imports e3nn, and e3nn's constants.pt is read
        at import time: registering `slice` after that point is too late to help. The
        checkpoint is declared loadable as a full pickle -- a MACE model cannot be
        reconstructed weights-only -- and its sha256 goes into the metadata below, so what
        was exempted is on the record rather than implied.
    """
    model_path = Path(model_path)
    if not model_path.is_file():
        raise FileNotFoundError(f"MACE checkpoint not found: {model_path}")
    if precision not in PRECISIONS:
        raise ValueError(f"precision must be one of {PRECISIONS}")

    # [1] before openmmml: e3nn reads constants.pt at import
    torch_guard.install()
    torch_guard.trust(model_path)

    import openmm
    from openmm import unit
    from openmmml import MLPotential

    if platform_name is None:
        platform_name = "CUDA" if device.startswith("cuda") else "CPU"

    potential = MLPotential("mace", modelPath=str(model_path))
    with torch_guard.moved_to_device():
        system = potential.createSystem(
            build_topology(numbers),
            device=device,
            precision=precision,
            returnEnergyType="energy",
            charge=charge,
            multiplicity=multiplicity,
            removeCMMotion=False,
        )
    force_names = [type(system.getForce(i)).__name__ for i in range(system.getNumForces())]
    if "CMMotionRemover" in force_names:
        raise RuntimeError("CMMotionRemover present; it would break momentum conservation")

    # No try/except: an unavailable platform must fail, not silently become CPU.
    platform = openmm.Platform.getPlatformByName(platform_name)
    integrator = openmm.VerletIntegrator(timestep_fs * unit.femtoseconds)
    context = openmm.Context(system, integrator, platform)
    metadata = {
        "backend": "openmm-ml/mace",
        "model_path": str(model_path),
        "model_sha256": _file_sha256(model_path),
        "device": device,
        "precision": precision,
        "platform": platform.getName(),
        "charge": charge,
        "multiplicity": multiplicity,
        "return_energy_type": "energy",
        "forces": force_names,
        "timestep_fs": timestep_fs,
        "neighborhood": "rebuilt from positions on every evaluation",
        "torch_load_guard": torch_guard.state(),
        "torch_load_guard_outermost": torch_guard.outermost(),
    }
    return context, system, metadata


def evaluate(context, positions_A):
    """Single-point energy and forces in eV and eV/Angstrom."""
    from openmm import unit

    positions_A = np.asarray(positions_A, dtype=float)
    context.setPositions(positions_A * NM_PER_ANGSTROM)
    state = context.getState(getEnergy=True, getForces=True)
    energy = (
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) * EV_PER_KJ_PER_MOL
    )
    forces = (
        state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        * EV_PER_KJ_PER_MOL
        * NM_PER_ANGSTROM
    )
    return float(energy), np.asarray(forces, dtype=float)


def _file_sha256(path):
    import hashlib

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
