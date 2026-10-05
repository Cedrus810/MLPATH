"""Optional geomeTRIC local refinement on ASE or an existing OpenMM Context.

The optimizer consumes Bohr/Hartree gradients; the project boundary remains
Angstrom/eV/forces. Convergence is reported separately from stationary-point
character and endpoint connectivity, neither of which this adapter certifies.
"""

from dataclasses import dataclass
import hashlib
from importlib.metadata import version
import logging
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

from .io import atomic_json, structure_hash


@dataclass(frozen=True)
class GeometricConfig:
    transition: bool = False
    maxiter: int = 300
    coordsys: str = "tric"
    convergence_set: str = "GAU"

    def validate(self):
        if not isinstance(self.transition, bool):
            raise ValueError("transition must be boolean")
        if (
            isinstance(self.maxiter, bool)
            or not isinstance(self.maxiter, int)
            or self.maxiter < 1
        ):
            raise ValueError("maxiter must be a positive integer")
        if self.coordsys not in ("tric", "dlc", "cart"):
            raise ValueError("coordsys must be tric, dlc or cart")
        if self.convergence_set not in ("GAU", "GAU_TIGHT", "GAU_VERYTIGHT"):
            raise ValueError("unsupported convergence_set")


def _validate_atoms(atoms):
    if len(atoms) < 2 or np.any(atoms.numbers <= 0):
        raise ValueError("geomeTRIC requires at least two real atoms")
    if np.any(atoms.pbc):
        raise ValueError("geomeTRIC refinement supports finite nonperiodic systems only")
    if atoms.constraints:
        raise ValueError("ASE constraints are not supported by this refinement adapter")
    if not np.isfinite(atoms.positions).all():
        raise ValueError("Non-finite input coordinates")


def ase_evaluator(atoms, calculator):
    """Evaluate a private atom copy, preserving charge/spin arrays and info.

    A caller can supply a GuardedCalculator to retain its chosen reliability
    policy. This low-level adapter validates finite values but sets no OOD policy.
    """
    _validate_atoms(atoms)
    work = atoms.copy()
    work.calc = calculator

    def evaluate(positions_A):
        work.set_positions(positions_A)
        return work.get_potential_energy(), work.get_forces()

    return evaluate


def openmm_evaluator(atoms, context):
    """Reuse the supplied Context, including its platform and MLP callbacks.

    Atom order must match the System. Reject constraints/virtual sites rather
    than optimize coordinates that the Hamiltonian does not freely expose.
    """
    _validate_atoms(atoms)
    system = context.getSystem()
    if system.getNumParticles() != len(atoms):
        raise ValueError("OpenMM particle count differs from input atom count")
    if system.usesPeriodicBoundaryConditions():
        raise ValueError("Periodic OpenMM Systems are not supported")
    if system.getNumConstraints() or any(system.isVirtualSite(i) for i in range(len(atoms))):
        raise ValueError("OpenMM constraints and virtual sites are not supported")
    from .openmm_backend import evaluate

    return lambda positions_A: evaluate(context, positions_A)


def make_engine(atoms, evaluator):
    """Build a geomeTRIC Engine without importing geomeTRIC for ordinary PRRS runs."""
    _validate_atoms(atoms)
    try:
        from geometric.engine import Engine
        from geometric.molecule import Molecule
        from geometric.nifty import bohr2ang, au2ev
    except ImportError as exc:
        raise ImportError(
            'Install geomeTRIC support with pip install -e ".[geometric]"'
        ) from exc

    molecule = Molecule()
    molecule.elem = atoms.get_chemical_symbols()
    molecule.xyzs = [atoms.positions.copy()]
    molecule.build_topology()

    class LocalEngine(Engine):
        def __init__(self):
            super().__init__(molecule)
            self.evaluations = 0

        def calc_new(self, coords, dirname):
            coords = np.asarray(coords, dtype=float)
            if coords.size != 3 * len(atoms) or not np.isfinite(coords).all():
                raise ValueError("Invalid geomeTRIC coordinates")
            energy, forces = evaluator(coords.reshape(-1, 3) * bohr2ang)
            forces = np.asarray(forces, dtype=float)
            energy = float(energy)
            if (
                forces.shape != (len(atoms), 3)
                or not np.isfinite(forces).all()
                or not np.isfinite(energy)
            ):
                raise ValueError("Evaluator returned invalid energy or forces")
            self.evaluations += 1
            return {"energy": energy / au2ev, "gradient": (-forces * bohr2ang / au2ev).ravel()}

    return LocalEngine()


def optimize_geometry(atoms, evaluator, output, config=None, metadata=None):
    """Refine one geometry and write a result, including explicit failed runs.

    Output must be new. ``optimization.extxyz`` contains optimizer iterates,
    not dynamics or every force call (TS Hessians also use displaced geometries).
    ``geometry.extxyz`` is the final optimizer frame, including on failure.
    """
    from dataclasses import asdict

    config = config or GeometricConfig()
    config.validate()
    engine = make_engine(atoms, evaluator)
    import geometric
    from geometric.errors import GeomOptNotConvergedError
    from geometric.internal import CartesianCoordinates, DelocalizedInternalCoordinates
    from geometric.nifty import ang2bohr, au2ev, bohr2ang
    from geometric.optimize import Optimizer
    from geometric.params import OptParams

    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write(str(output / "input.extxyz"), atoms)
    record = {
        "status": "RUNNING",
        "optimizer": "geomeTRIC",
        "geometric_version": geometric.__version__,
        "versions": {name: version(name) for name in ("ase", "numpy", "scipy")},
        "implementation_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("geometric_backend.py", "geometric_cli.py", "openmm_backend.py")
        },
        "config": asdict(config),
        "input_structure_sha256": structure_hash(atoms),
        "backend": dict(metadata or {}),
        "stationary_point_certified": False,
        "connectivity_verified": False,
        "output": str(output),
    }
    atomic_json(output / "result.json", record)
    logger = logging.getLogger("geometric")
    handler = logging.FileHandler(output / "geometric.log")
    old_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    optimizer = None
    try:
        cls = (
            CartesianCoordinates
            if config.coordsys == "cart"
            else DelocalizedInternalCoordinates
        )
        internal = cls(engine.M, build=True, connect=config.coordsys == "dlc", addcart=False)
        params = OptParams(
            transition=config.transition,
            maxiter=config.maxiter,
            convergence_set=config.convergence_set,
            hessian="first" if config.transition else "never",
            frequency=False,
        )
        optimizer = Optimizer(
            atoms.positions.ravel() * ang2bohr,
            engine.M,
            internal,
            engine,
            str(output / "scratch"),
            params,
        )
        optimizer.optimizeGeometry()
        record["status"] = "CONVERGED"
    except GeomOptNotConvergedError as exc:
        record.update(status="NOT_CONVERGED", reason=str(exc))
    except Exception as exc:
        record.update(status="FAILED", reason=f"{type(exc).__name__}: {exc}")
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(old_level)

    frames = []
    if optimizer is not None:
        progress = optimizer.progress
        for positions, energy, gradient in zip(
            progress.xyzs, progress.qm_energies, progress.qm_grads
        ):
            frame = atoms.copy()
            frame.positions = positions
            forces = -np.asarray(gradient).reshape(-1, 3) * au2ev / bohr2ang
            frame.calc = SinglePointCalculator(
                frame, energy=float(energy * au2ev), forces=forces
            )
            frames.append(frame)
        record["iterations"] = optimizer.Iteration
    record["force_evaluations"] = engine.evaluations
    if frames:
        write(str(output / "optimization.extxyz"), frames)
        write(str(output / "geometry.extxyz"), frames[-1])
        record.update(
            geometry="geometry.extxyz",
            trajectory="optimization.extxyz",
            energy_eV=float(frames[-1].get_potential_energy()),
            fmax_eV_A=float(np.linalg.norm(frames[-1].get_forces(), axis=1).max()),
        )
    atomic_json(output / "result.json", record)
    return record
