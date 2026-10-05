"""Units, real optimization, and failure semantics on analytic surfaces."""

import json

import numpy as np
import pytest
from ase import Atoms
from ase.io import read, write

pytest.importorskip("geometric")

from prrs.calculators import DoubleWell  # noqa: E402
from prrs.cli import main  # noqa: E402
from prrs.geometric_backend import (  # noqa: E402
    GeometricConfig,
    ase_evaluator,
    make_engine,
    openmm_evaluator,
    optimize_geometry,
)


def dimer(distance):
    return Atoms("C2", positions=[[0, 0, 0], [distance, 0, 0]])


def test_atomic_unit_gradient_matches_energy_finite_differences():
    from geometric.nifty import ang2bohr

    atoms = dimer(1.4)
    atoms.info["charge"] = -1
    original = atoms.positions.copy()
    engine = make_engine(atoms, ase_evaluator(atoms, DoubleWell()))
    coords = atoms.positions.ravel() * ang2bohr
    analytic = engine.calc_new(coords, ".")["gradient"]
    numeric = np.zeros_like(coords)
    for index in range(len(coords)):
        delta = np.zeros_like(coords)
        delta[index] = 1e-5
        numeric[index] = (
            engine.calc_new(coords + delta, ".")["energy"]
            - engine.calc_new(coords - delta, ".")["energy"]
        ) / 2e-5
    np.testing.assert_allclose(analytic, numeric, atol=1e-10, rtol=1e-7)
    np.testing.assert_array_equal(atoms.positions, original)


@pytest.mark.parametrize(
    "transition,initial,target,energy", [(False, 1.4, 1.2, 0), (True, 1.75, 1.8, 0.4)]
)
def test_real_geometric_minimum_and_saddle(tmp_path, transition, initial, target, energy):
    atoms = dimer(initial)
    result = optimize_geometry(
        atoms,
        ase_evaluator(atoms, DoubleWell()),
        tmp_path / "opt",
        GeometricConfig(transition=transition, convergence_set="GAU_TIGHT"),
    )
    assert result["status"] == "CONVERGED", result
    final = read(tmp_path / "opt" / result["geometry"])
    assert final.get_distance(0, 1) == pytest.approx(target, abs=1e-4)
    assert result["energy_eV"] == pytest.approx(energy, abs=1e-7)
    assert result["stationary_point_certified"] is False
    assert result["connectivity_verified"] is False
    assert result["force_evaluations"] >= 1
    if transition:
        assert result["force_evaluations"] >= 12  # central differences for 6 coordinates
    assert json.loads((tmp_path / "opt" / "result.json").read_text()) == result


def test_iteration_budget_cannot_be_success(tmp_path):
    atoms = dimer(1.5)
    result = optimize_geometry(
        atoms, ase_evaluator(atoms, DoubleWell()), tmp_path / "opt", GeometricConfig(maxiter=1)
    )
    assert result["status"] == "NOT_CONVERGED"
    assert (tmp_path / "opt" / "geometry.extxyz").is_file()
    with pytest.raises(FileExistsError):
        optimize_geometry(atoms, ase_evaluator(atoms, DoubleWell()), tmp_path / "opt")


def test_nonfinite_evaluation_is_recorded_as_failure(tmp_path):
    result = optimize_geometry(
        dimer(1.4), lambda x: (np.nan, np.zeros_like(x)), tmp_path / "opt"
    )
    assert result["status"] == "FAILED"
    assert "invalid energy or forces" in result["reason"]
    assert "geometry" not in result


def test_refuse_periodic_constraints_and_bad_settings(tmp_path):
    from ase.constraints import FixAtoms

    atoms = dimer(1.4)
    atoms.pbc = True
    with pytest.raises(ValueError, match="nonperiodic"):
        make_engine(atoms, None)
    atoms.pbc = False
    atoms.set_constraint(FixAtoms(indices=[0]))
    with pytest.raises(ValueError, match="constraints"):
        make_engine(atoms, None)
    with pytest.raises(ValueError, match="positive integer"):
        GeometricConfig(maxiter=0).validate()


def test_openmm_reference_units_and_optimization(tmp_path):
    openmm = pytest.importorskip("openmm")
    atoms = dimer(1.5)
    system = openmm.System()
    for _ in atoms:
        system.addParticle(12.0)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, 0.12, 1000.0)
    system.addForce(bond)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    evaluator = openmm_evaluator(atoms, context)
    energy, forces = evaluator(atoms.positions)
    assert energy == pytest.approx(0.45 / 96.4853)
    assert forces[0, 0] == pytest.approx(3.0 / 96.4853)
    result = optimize_geometry(
        atoms, evaluator, tmp_path / "mm", GeometricConfig(convergence_set="GAU_TIGHT")
    )
    assert result["status"] == "CONVERGED", result
    assert read(tmp_path / "mm" / "geometry.extxyz").get_distance(0, 1) == pytest.approx(
        1.2, abs=1e-4
    )
    assert context.getPlatform().getName() == "Reference"


def test_cli_ase_minimum_and_budget_failure(tmp_path):
    source = tmp_path / "seed.extxyz"
    write(source, dimer(1.4))
    args = [
        "optimize",
        str(source),
        "--calculator",
        "prrs.calculators:double_well_factory",
        "--output",
    ]
    assert main(args + [str(tmp_path / "cli")]) == 0
    assert main(args + [str(tmp_path / "short"), "--maxiter", "1"]) == 2


def test_charge_and_spin_are_preserved_or_conflicts_refused():
    from prrs.geometric_cli import _electronic_state

    atoms = dimer(1.4)
    atoms.info.update(charge=-1, spin=2)
    _electronic_state(atoms, None, None, defaults=True)
    assert atoms.info == {"charge": -1, "spin": 2}
    with pytest.raises(ValueError, match="conflicts"):
        _electronic_state(atoms, 0, None)
    with pytest.raises(ValueError, match="at least one"):
        _electronic_state(dimer(1.4), None, 0)


def test_cli_rejects_a_whole_path_instead_of_optimizing_its_last_endpoint(tmp_path):
    source = tmp_path / "path.extxyz"
    write(source, [dimer(1.2), dimer(1.8), dimer(2.4)])
    output = tmp_path / "wrong_guess"
    assert (
        main(
            [
                "optimize",
                str(source),
                "--transition",
                "--calculator",
                "prrs.calculators:double_well_factory",
                "--output",
                str(output),
            ]
        )
        == 2
    )
    assert not output.exists()


def test_cli_openmm_serialized_system(tmp_path):
    openmm = pytest.importorskip("openmm")
    from openmm import app, unit
    from prrs.openmm_backend import build_topology

    atoms = dimer(1.5)
    source = tmp_path / "seed.pdb"
    topology = build_topology(atoms.numbers)
    for i, atom in enumerate(topology.atoms()):
        atom.name = f"C{i + 1}"
    with source.open("w") as handle:
        app.PDBFile.writeFile(topology, atoms.positions * unit.angstrom, handle)
    system = openmm.System()
    for _ in atoms:
        system.addParticle(12.0)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, 0.12, 1000.0)
    system.addForce(bond)
    xml = tmp_path / "system.xml"
    xml.write_text(openmm.XmlSerializer.serialize(system))
    assert (
        main(
            [
                "optimize",
                str(source),
                "--openmm-xml",
                str(xml),
                "--platform",
                "Reference",
                "--output",
                str(tmp_path / "xml"),
            ]
        )
        == 0
    )
    result = json.loads((tmp_path / "xml" / "result.json").read_text())
    assert result["backend"]["platform"] == "Reference"


def test_cli_openmm_forcefield_water_has_no_hidden_constraints(tmp_path):
    pytest.importorskip("openmm")
    from pathlib import Path
    from ase.build import molecule
    from openmm import app, unit

    atoms = molecule("H2O")
    topology = app.Topology()
    residue = topology.addResidue("HOH", topology.addChain())
    names = ["O", "H1", "H2"]
    members = [
        topology.addAtom(name, app.Element.getByAtomicNumber(int(number)), residue)
        for name, number in zip(names, atoms.numbers)
    ]
    topology.addBond(members[0], members[1])
    topology.addBond(members[0], members[2])
    source = tmp_path / "water.pdb"
    with source.open("w") as handle:
        app.PDBFile.writeFile(topology, atoms.positions * unit.angstrom, handle)
    xml = Path(app.__file__).parent / "data" / "amber14" / "tip3p.xml"
    assert (
        main(
            [
                "optimize",
                str(source),
                "--openmm-xml",
                str(xml),
                "--platform",
                "Reference",
                "--output",
                str(tmp_path / "water"),
            ]
        )
        == 0
    )


def test_openmm_constraints_are_refused():
    openmm = pytest.importorskip("openmm")
    system = openmm.System()
    system.addParticle(12.0)
    system.addParticle(12.0)
    system.addConstraint(0, 1, 0.12)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    with pytest.raises(ValueError, match="constraints"):
        openmm_evaluator(dimer(1.2), context)
