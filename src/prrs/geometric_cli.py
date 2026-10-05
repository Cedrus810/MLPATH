"""CLI wiring for optional local refinement; imports backends only on use."""

import json
from pathlib import Path

from ase.io import read

from .geometric_backend import (
    GeometricConfig,
    ase_evaluator,
    openmm_evaluator,
    optimize_geometry,
)


def add_parser(sub):
    parser = sub.add_parser(
        "optimize", help="Local geomeTRIC minimum or TS-candidate refinement"
    )
    parser.add_argument("input", help="Single geometry; atom order must match the backend")
    parser.add_argument("--output", required=True, help="New output directory")
    backend = parser.add_mutually_exclusive_group(required=True)
    backend.add_argument("--calculator", help="ASE calculator factory module:function")
    backend.add_argument(
        "--mace-model", help="Local checkpoint for the existing MACE/OpenMM backend"
    )
    backend.add_argument(
        "--openmm-xml", help="Serialized System XML or ForceField XML (input must be PDB)"
    )
    parser.add_argument("--calculator-kwargs", help="JSON file of ASE factory arguments")
    parser.add_argument(
        "--transition", action="store_true", help="Search for a TS from the input guess"
    )
    parser.add_argument("--maxiter", type=int, default=300)
    parser.add_argument("--coordsys", choices=("tric", "dlc", "cart"), default="tric")
    parser.add_argument(
        "--convergence-set", choices=("GAU", "GAU_TIGHT", "GAU_VERYTIGHT"), default="GAU"
    )
    parser.add_argument("--device", default="cpu", help="MACE torch device")
    parser.add_argument("--precision", choices=("single", "double"), default="double")
    parser.add_argument("--platform", help="Explicit OpenMM platform; no fallback")
    parser.add_argument("--charge", type=int, help="Total charge; must agree with input info")
    parser.add_argument(
        "--multiplicity", type=int, help="Spin multiplicity; must agree with input info"
    )


def _electronic_state(atoms, charge, multiplicity, defaults=False):
    from numbers import Integral

    for key, explicit, fallback in (("charge", charge, 0), ("spin", multiplicity, 1)):
        stored = atoms.info.get(key)
        if stored is not None and (
            isinstance(stored, bool) or not isinstance(stored, Integral)
        ):
            raise ValueError(f"input {key} must be an integer")
        if explicit is not None and stored is not None and explicit != stored:
            raise ValueError(f"input {key}={stored} conflicts with requested {explicit}")
        value = explicit if explicit is not None else stored
        if value is None and defaults:
            value = fallback
        if key == "spin" and value is not None and value < 1:
            raise ValueError("multiplicity must be at least one")
        if value is not None:
            atoms.info[key] = int(value)


def run(args):
    from .cli import factory_metadata, load_json
    from .calculators import load_factory
    from .openmm_backend import _file_sha256

    config = GeometricConfig(args.transition, args.maxiter, args.coordsys, args.convergence_set)
    config.validate()
    if Path(args.output).exists():
        raise FileExistsError(f"Will not overwrite {args.output}")
    frames = read(args.input, index=":")
    if len(frames) != 1:
        raise ValueError(
            "Input must contain one geometry; extract the endpoint or TS guess explicitly"
        )
    atoms = frames[0]
    if args.openmm_xml and (args.charge is not None or args.multiplicity is not None):
        raise ValueError(
            "OpenMM XML defines its own Hamiltonian; charge/spin overrides are unsupported"
        )
    _electronic_state(atoms, args.charge, args.multiplicity, defaults=bool(args.mace_model))
    if args.calculator_kwargs and not args.calculator:
        raise ValueError("--calculator-kwargs requires --calculator")
    if args.calculator:
        kwargs = load_json(args.calculator_kwargs) if args.calculator_kwargs else {}
        evaluator = ase_evaluator(atoms, load_factory(args.calculator, kwargs)())
        metadata = {"engine": "ase", **factory_metadata(args.calculator, kwargs)}
    elif args.mace_model:
        from .openmm_backend import create_mace_context

        context, system, metadata = create_mace_context(
            atoms.numbers,
            args.mace_model,
            device=args.device,
            precision=args.precision,
            platform_name=args.platform,
            charge=atoms.info["charge"],
            multiplicity=atoms.info["spin"],
        )
        evaluator = openmm_evaluator(atoms, context)
        import openmm

        metadata["openmm_version"] = openmm.__version__
    else:
        import openmm
        from openmm import app, unit

        if Path(args.input).suffix.lower() != ".pdb":
            raise ValueError("--openmm-xml requires a PDB input carrying the topology")
        pdb = app.PDBFile(args.input)
        xml_path = Path(args.openmm_xml).resolve()
        text = xml_path.read_text()
        import xml.etree.ElementTree as ET

        if ET.fromstring(text).tag == "System":
            system = openmm.XmlSerializer.deserialize(text)
        else:
            system = app.ForceField(str(xml_path)).createSystem(
                pdb.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=None,
                rigidWater=False,
                removeCMMotion=False,
            )
        # Use OpenMM's PDB parsing for coordinate/element order on this route.
        import numpy as np
        from ase import Atoms

        atoms = Atoms(
            numbers=[atom.element.atomic_number for atom in pdb.topology.atoms()],
            positions=np.asarray(pdb.positions.value_in_unit(unit.angstrom)),
        )
        integrator = openmm.VerletIntegrator(0.25 * unit.femtoseconds)
        platform = openmm.Platform.getPlatformByName(args.platform or "CPU")
        context = openmm.Context(system, integrator, platform)
        evaluator = openmm_evaluator(atoms, context)
        metadata = {
            "engine": "openmm",
            "platform": platform.getName(),
            "openmm_version": openmm.__version__,
            "xml_path": str(xml_path),
            "xml_sha256": _file_sha256(xml_path),
        }
    metadata["input_file_sha256"] = _file_sha256(args.input)
    result = optimize_geometry(atoms, evaluator, args.output, config, metadata)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["status"] == "CONVERGED" else 2
