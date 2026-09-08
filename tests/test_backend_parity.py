"""P1 acceptance: the ASE and OpenMM backends must be the same physical Hamiltonian.

Skipped unless a local MACE checkpoint is present. Set PRRS_MACE_MODEL to override the
default location. The throughput comparison is opt-in because it costs tens of seconds.
"""

import os
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip("openmm")
pytest.importorskip("openmmml")
pytest.importorskip("mace")
torch = pytest.importorskip("torch")

from prrs.openmm_backend import (  # noqa: E402 - must follow the importorskip guards
    build_topology,
    create_mace_context,
    evaluate,
)

MODEL = Path(
    os.environ.get("PRRS_MACE_MODEL", Path.home() / ".cache/mace/MACE-OFF24_medium.model")
)
pytestmark = pytest.mark.skipif(not MODEL.is_file(), reason=f"no MACE checkpoint at {MODEL}")

# Absolute tolerances frozen here, not tuned per run: single precision would need
# looser values and a different decision, so both backends run in double.
ENERGY_TOL_EV = 1e-4
FORCE_TOL_EV_A = 1e-3


def _device():
    return "cuda" if torch.cuda.is_available() else "cpu"


@pytest.fixture(scope="module")
def system():
    """A fixed geometry set spanning equilibrium, strain, a broken bond and disorder."""
    from ase.build import molecule

    mol = molecule("CH3CH2OH")
    numbers, positions = mol.numbers.copy(), mol.positions.copy()
    stretched = positions.copy()
    stretched[0] += np.array([0.8, 0.0, 0.0])
    broken = positions.copy()
    broken[7] += np.array([0.0, 0.0, 2.5])
    rattled = positions + np.random.default_rng(0).normal(scale=0.15, size=positions.shape)
    return numbers, {
        "equilibrium": positions,
        "stretched": stretched,
        "broken": broken,
        "rattled": rattled,
    }


@pytest.fixture(scope="module")
def openmm_context(system):
    numbers, _ = system
    context, _, metadata = create_mace_context(numbers, MODEL, device=_device())
    assert metadata["forces"] == ["PythonForce"], metadata["forces"]
    return context


@pytest.fixture(scope="module")
def ase_calculator():
    from mace.calculators import MACECalculator

    return MACECalculator(model_paths=str(MODEL), device=_device(), default_dtype="float64")


def _ase_eval(numbers, positions, calculator):
    from ase import Atoms

    atoms = Atoms(numbers=numbers, positions=positions)
    atoms.calc = calculator
    return atoms.get_potential_energy(), atoms.get_forces()


def test_backends_agree_on_energy_and_forces(system, openmm_context, ase_calculator):
    numbers, geometries = system
    for label, positions in geometries.items():
        energy_openmm, forces_openmm = evaluate(openmm_context, positions)
        energy_ase, forces_ase = _ase_eval(numbers, positions, ase_calculator)
        assert abs(energy_openmm - energy_ase) < ENERGY_TOL_EV, label
        assert np.abs(forces_openmm - forces_ase).max() < FORCE_TOL_EV_A, label


def test_context_survives_a_bond_graph_change(system, openmm_context):
    """A used Context and an untouched one must agree on a broken-bond geometry.

    This is what licenses reusing one Context across a reaction instead of rebuilding
    it whenever connectivity changes.
    """
    numbers, geometries = system
    for positions in geometries.values():
        evaluate(openmm_context, positions)
    used_energy, used_forces = evaluate(openmm_context, geometries["broken"])
    fresh, _, _ = create_mace_context(numbers, MODEL, device=_device())
    fresh_energy, fresh_forces = evaluate(fresh, geometries["broken"])
    assert used_energy == pytest.approx(fresh_energy, abs=1e-9)
    assert np.abs(used_forces - fresh_forces).max() < 1e-9


def test_forces_match_finite_differences(system, openmm_context):
    numbers, geometries = system
    positions = geometries["equilibrium"]
    analytic = evaluate(openmm_context, positions)[1]
    step = 1e-4
    for atom in (0, 5):
        for axis in range(3):
            plus, minus = positions.copy(), positions.copy()
            plus[atom, axis] += step
            minus[atom, axis] -= step
            numeric = -(
                evaluate(openmm_context, plus)[0] - evaluate(openmm_context, minus)[0]
            ) / (2 * step)
            assert numeric == pytest.approx(analytic[atom, axis], abs=1e-3)


def test_no_bonds_are_declared_to_the_model(system):
    numbers, _ = system
    topology = build_topology(numbers)
    assert topology.getNumAtoms() == len(numbers)
    assert list(topology.bonds()) == []


def test_platform_is_never_silently_downgraded(system):
    numbers, _ = system
    with pytest.raises(Exception):
        create_mace_context(numbers, MODEL, device="cpu", platform_name="NoSuchPlatform")


@pytest.mark.skipif(not torch.cuda.is_available(), reason="needs CUDA")
def test_upstream_still_needs_the_explicit_device_move(system):
    """Pin the openmm-ml workaround so a fixed upstream is noticed instead of assumed.

    If this starts failing, openmm-ml now moves a locally loaded model to the requested
    device itself and the torch_guard.moved_to_device scope around createSystem can be
    deleted. The move is deliberately scoped rather than permanent so that this test can
    still observe the unpatched behaviour.
    """
    from openmm import unit
    import openmm
    from openmmml import MLPotential

    numbers, geometries = system
    potential = MLPotential("mace", modelPath=str(MODEL))
    unpatched = potential.createSystem(
        build_topology(numbers),
        device="cuda",
        precision="double",
        returnEnergyType="energy",
        charge=0,
        multiplicity=1,
        removeCMMotion=False,
    )
    context = openmm.Context(
        unpatched,
        openmm.VerletIntegrator(0.25 * unit.femtoseconds),
        openmm.Platform.getPlatformByName("CUDA"),
    )
    with pytest.raises(Exception):
        evaluate(context, geometries["equilibrium"])


@pytest.mark.skipif(
    not os.environ.get("PRRS_RUN_BENCHMARKS"),
    reason="set PRRS_RUN_BENCHMARKS=1 to time the backends",
)
def test_report_backend_throughput(system, openmm_context, ase_calculator, capsys):
    """Reported, not asserted: throughput is evidence for a decision, not a pass criterion."""
    import time
    from ase import Atoms, units
    from ase.md.verlet import VelocityVerlet

    numbers, geometries = system
    positions = geometries["equilibrium"]
    steps = 100
    integrator = openmm_context.getIntegrator()
    openmm_context.setPositions(positions * 0.1)
    integrator.step(5)
    start = time.perf_counter()
    integrator.step(steps)
    openmm_context.getState()
    openmm_seconds = time.perf_counter() - start

    atoms = Atoms(numbers=numbers, positions=positions)
    atoms.calc = ase_calculator
    dynamics = VelocityVerlet(atoms, timestep=0.25 * units.fs, logfile=None)
    dynamics.run(5)
    start = time.perf_counter()
    dynamics.run(steps)
    ase_seconds = time.perf_counter() - start
    with capsys.disabled():
        print(
            f"\n  openmm {1000 * openmm_seconds / steps:.2f} ms/step | "
            f"ase {1000 * ase_seconds / steps:.2f} ms/step | "
            f"ase/openmm {ase_seconds / openmm_seconds:.2f}x"
        )
