"""Write PRRS's analytic Cartesian Hessian as an ORCA .hess file.

Why this exists: ORCA's ExtOpt interface cannot produce a Hessian. `Freq`, `NumFreq`
and `Calc_Hess` all fail the same way -- `ERROR (SHARK): Failed to read input file` --
because every one of them routes through property integrals, and an external program
has no wavefunction. So a transition-state search over an external potential gets a
model Hessian with no negative mode, picks a mode to invert by guesswork, and walks
away from the saddle: b04's first T[E_m] attempt started 0.001 Eh/Bohr from a
stationary point and left it, RMS gradient climbing 0.0010 -> 0.0133 over 50 cycles
with every step pinned at the trust radius.

`%geom InHess Read InHessName "x.hess" end` takes the Hessian from a file instead,
and this writes that file from the same analytic second derivative PRRS itself uses
(`runner._analytic_hessian`, whose units are proven by a Hessian-vector product
against finite differences rather than assumed).

This does NOT weaken the independence T[E_m] rests on. The independence that matters
is the OPTIMIZER and the STARTING GEOMETRY -- RS-P-RFO in internal coordinates from a
guess taken from no previous result. A Hessian is a property of the potential energy
surface, the same surface both operators already query for gradients; handing it over
is not handing over PRRS's answer.

    python scripts/prrs_to_orca_hess.py <geometry> <out.hess> \
        --calculator module:function [--charge -1] [--multiplicity 1]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from prrs.calculators import load_factory  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import _analytic_hessian  # noqa: E402

EV_PER_HARTREE = 27.211386245988
BOHR_PER_A = 1.0 / 0.529177210903


def orca_hess_text(atoms, hessian_eV_A2, energy_eV=0.0):
    """ORCA .hess: Hessian in Eh/Bohr^2, coordinates in Bohr, five columns per block."""
    scale = (1.0 / EV_PER_HARTREE) / BOHR_PER_A**2
    h = np.asarray(hessian_eV_A2, dtype=float) * scale
    n = h.shape[0]
    lines = ["$orca_hessian_file", "", "$act_atom", "  0", "", "$act_coord", "  0", ""]
    lines += ["$act_energy", f"{energy_eV / EV_PER_HARTREE:20.10f}", "", "$hessian", f"{n}"]
    for start in range(0, n, 5):
        block = range(start, min(start + 5, n))
        lines.append("".join(f"{c:>18d}" for c in block))
        for row in range(n):
            lines.append(f"{row:>7d}" + "".join(f"{h[row, c]:>18.10E}" for c in block))
    lines += ["", "$atoms", f"{len(atoms)}"]
    for symbol, mass, position in zip(
        atoms.get_chemical_symbols(), atoms.get_masses(), atoms.get_positions()
    ):
        x, y, z = position * BOHR_PER_A
        lines.append(f" {symbol:<3s} {mass:12.6f} {x:18.10f} {y:18.10f} {z:18.10f}")
    lines += ["", "$end", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometry")
    parser.add_argument("output")
    parser.add_argument("--calculator", required=True, help="module:function")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    args = parser.parse_args(argv)

    atoms = read(args.geometry)
    # The charge reaches a charge-aware model through atoms.info, exactly as
    # search._stamped_with_charge puts it there; a charge-blind model ignores it.
    atoms.info["charge"] = args.charge
    atoms.info["spin"] = args.multiplicity - 1
    # SearchConfig refuses total_charge unless the potential is declared to respond to
    # it, which is the gate b04's preflight 3 is about -- so a nonzero charge implies it.
    config = SearchConfig(
        charge_sensitive=bool(args.charge),
        total_charge=args.charge,
        multiplicity=args.multiplicity,
    )
    physical = load_factory(args.calculator)()

    hessian, provenance = _analytic_hessian(atoms, physical, config)
    if hessian is None:
        raise SystemExit(f"no analytic Hessian: {provenance}")
    atoms.calc = physical
    energy = float(atoms.get_potential_energy())

    Path(args.output).write_text(orca_hess_text(atoms, hessian, energy))
    eigenvalues = np.linalg.eigvalsh(hessian)
    print(f"provenance: {provenance}")
    print(f"lowest six eigenvalues (eV/A^2): {np.round(eigenvalues[:6], 6)}")
    print(f"negative modes: {int((eigenvalues < -1e-6).sum())}")
    print(f"written {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
