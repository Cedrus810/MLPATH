"""ORCA ExtOpt bridge  (launched by scripts/mace_engrad.sh, not directly)

ORCA execs ProgExt with its own PATH, which is not the conda environment, so a
bare `#!/usr/bin/env python` shebang finds a python without numpy. The .sh
wrapper next to this file names the interpreter explicitly.

ORCA ExtOpt bridge: let ORCA's optimizer drive a MACE potential.

This exists to make T[E_m] and T[E_QM] come from the SAME optimizer. Comparing
geomeTRIC-on-MLP against ORCA-on-DFT would fold "different optimizer" into
"different potential energy surface", and those are the two things the
comparison is supposed to separate (benchmarks/README.md rule 2).

Protocol, established by execution on 2026-09-04 rather than read from docs:

    ORCA calls   <this script> <base>_EXT.extinp.tmp
    that file =  <base>_EXT.xyz # xyz filename
                 0              # charge
                 1              # multiplicity
                 1              # NCores
                 1              # do gradient
    we must write <base>_EXT.engrad, energies in Eh, gradient in Eh/Bohr

The two file formats are handled by `opi.external_methods.interface.ExtoptInterface`
rather than by hand. Hand-rolled parsing of a vendor format is the part that rots
silently when the vendor changes it, and ORCA ships the reader itself.

Model is chosen by environment so one ORCA input can be pointed at any
checkpoint without editing anything:

    PRRS_ENGRAD_MODEL   path to the .model file          (required)
    PRRS_ENGRAD_TYPE    MACE | PolarMACE                 (default MACE)
    PRRS_ENGRAD_DEVICE  cuda | cpu                       (default cuda)

The charge and multiplicity ORCA passes are written into atoms.info, because a
charge-aware model reads them from there and a charge-blind one silently
ignores them -- which is exactly the failure b04's preflight gate 3 is about.
"""

import os
import sys
from pathlib import Path

import numpy as np

HARTREE_EV = 27.211386245988
BOHR_A = 0.529177210903


def _evaluate(atoms, charge, mult, gradient):
    """Ask the resident server if one is up, otherwise load the model here.

    A numerical Hessian is 6N gradient calls and loading a checkpoint costs about
    5 s against milliseconds of evaluation, so the server is the path that matters.
    The fallback exists so a single-point call still works with no server running,
    and both routes read the same PRRS_ENGRAD_* variables -- they cannot end up
    disagreeing about which checkpoint was used.
    """
    address = os.environ.get("PRRS_ENGRAD_SERVER")
    if address:
        from multiprocessing.connection import Client

        host, _, port = address.rpartition(":")
        with Client((host or "127.0.0.1", int(port))) as conn:
            conn.send({"type": "engrad",
                       "numbers": [int(z) for z in atoms.numbers],
                       "positions": atoms.positions.tolist(),
                       "charge": charge, "spin": mult, "gradient": gradient})
            reply = conn.recv()
        return reply["energy_eV"], reply["forces_eV_A"]

    model = os.environ["PRRS_ENGRAD_MODEL"]
    from prrs import torch_guard

    torch_guard.install(trusted=[model])
    from mace.calculators import MACECalculator

    probe = atoms.copy()
    probe.info["charge"], probe.info["spin"] = charge, mult
    probe.calc = MACECalculator(
        model_paths=model,
        model_type=os.environ.get("PRRS_ENGRAD_TYPE", "MACE"),
        device=os.environ.get("PRRS_ENGRAD_DEVICE", "cuda"),
        default_dtype="float64",
    )
    return float(probe.get_potential_energy()), (probe.get_forces() if gradient else None)


def main(argv):
    control = Path(argv[1])
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from ase.io import read
    from opi.external_methods.interface import ExtoptInterface

    orca = ExtoptInterface()
    xyz, charge, mult, _ncores, do_gradient, _extra = orca.read_extopt_input(control)
    base = control.name[: -len(".extinp.tmp")]

    atoms = read(control.parent / xyz)
    energy_eV, forces = _evaluate(atoms, charge, mult, bool(do_gradient))

    energy = energy_eV / HARTREE_EV
    # ORCA wants the gradient, ASE returns forces: gradient = -force.
    gradient = (-np.asarray(forces) / HARTREE_EV * BOHR_A
                if do_gradient else np.zeros((len(atoms), 3)))

    orca.write_orca_input(
        control.parent / f"{base}.engrad",
        nat=len(atoms),
        etot=energy,
        grad=list(gradient.reshape(-1)) if do_gradient else None,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
