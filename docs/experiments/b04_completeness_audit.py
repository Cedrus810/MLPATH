"""b04: the first completeness audit on a real MLP surface (PLAN item 3, production).

`tests/test_completeness.py` exercises the full chain on an ANALYTIC potential. This
script is the same chain on POLAR-1, against a minimum b04 actually admitted -- the
first time delta / L / V / R and a membership table are measured on a surface the
project reports results from.

It issues nothing. Decision 4 (2026-09-20) left the implication with the mathematics,
so `audit_basin` stamps every text DRAFT and this script only writes the record. What
it produces is the measurement the ratification argument would have to stand on.

The subspace V (d <= 3) is the three softest non-trivial mass-weighted modes at the
minimum, and that construction travels inside the record: a competing basin, if one
is near, is nearest along the directions the surface resists least. It is a choice,
not a derivation, which is exactly why it is recorded rather than assumed.

R / delta / X / gap_w are DECLARATIONS, written here before the run and not tuned
after it. Cost is exponential in d: at R = 0.5 A and delta = 0.1 A the grid spacing is
2*delta/sqrt(3) = 0.115 A, so the ball holds a few hundred points and each one is a
quench. Raising delta is cheaper and weaker; both facts belong in the record.

Needs the POLAR-1 checkpoint: run on kasuga02, not yayoi.

    python docs/experiments/b04_completeness_audit.py [minimum.extxyz]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.completeness import audit_basin, certificate_text  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import hessian_spectrum  # noqa: E402

DEFAULT_MINIMUM = ROOT / "runs" / "b04_frozen2" / "chemical" / "c0000_m0000.extxyz"
OUT = ROOT / "runs" / "b04_completeness_audit"

# Declared before the run. See the module docstring.
RADIUS_A = 0.5
DELTA_A = 0.1
GAP_W_A = 0.1
X_BELOW_MINIMUM_eV = 0.05  # rule out basins deeper than E(a*) - this
# Below this an eigenvalue is the projected-out rigid null space, not a soft mode.
NULL_EIGENVALUE = 1e-8
# The basin tolerance is NOT declared here: audit_basin reads config.basin_rmsd_A, so
# membership uses the same number the frozen b04 run used. One place decides it.


def load_config():
    text = (ROOT / "benchmarks" / "b04_sn2_chloride" / "config.json").read_text()
    return SearchConfig(
        **json.loads(
            "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
        )
    )


def soft_subspace(atoms, factory, config, dimension=3):
    """The `dimension` softest non-trivial modes, converted to CARTESIAN directions.

    hessian_spectrum returns eigenvalues ascending with the rigid modes already
    projected out (the same projection the frozen b04 results used), so "non-trivial"
    means what it returns, not a re-filtering.

    The eigenvectors are mass-weighted; `cover` displaces positions in Angstrom. Feeding
    one to the other is the silent unit error `polish_soft_modes` footnote [1] names, so
    the conversion x = v / sqrt(m) is done here, explicitly, before normalising.

    hessian_spectrum PROJECTS the rigid modes out, it does not remove their dimensions:
    they come back as eigenvalues at ~1e-15 with arbitrary eigenvectors. Taking the
    lowest three without filtering therefore hands back the null space -- which is what
    the 2026-09-22 first run did, and why its 484 `other_basin` rows said nothing about
    the surface. NULL_EIGENVALUE is the cut, and how many modes it removed is returned
    so the record can state that six were found where six are expected.
    """
    values, vectors = hessian_spectrum(atoms, factory, config, return_vectors=True)
    values = np.asarray(values)
    live = np.flatnonzero(np.abs(values) > NULL_EIGENVALUE)
    if len(live) < dimension:
        raise ValueError(f"only {len(live)} non-null modes; cannot build a d={dimension} V")
    weights = np.repeat(np.sqrt(atoms.get_masses()), 3)
    basis = []
    for i in live[:dimension]:
        cartesian = np.asarray(vectors)[:, i].ravel() / weights
        basis.append(cartesian / np.linalg.norm(cartesian))
    return np.array(basis), [float(values[i]) for i in live[:5]], len(values) - len(live)


def main(argv):
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT_MINIMUM
    atoms = read(path)
    config = load_config()
    OUT.mkdir(parents=True, exist_ok=True)

    basis, spectrum, null_modes = soft_subspace(atoms, polar1, config)
    atoms.calc = polar1()
    energy = float(atoms.get_potential_energy())

    record = audit_basin(
        atoms,
        polar1,
        config,
        basis=basis,
        radius=RADIUS_A,
        delta=DELTA_A,
        energy_ceiling_X=energy - X_BELOW_MINIMUM_eV,
        gap_w=GAP_W_A,
    )
    record["subspace_construction"] = (
        f"three softest NON-NULL modes at the minimum, mass-weighted eigenvectors "
        f"converted to Cartesian; {null_modes} null modes removed at |lambda| <= "
        f"{NULL_EIGENVALUE} (six expected for a non-linear molecule); lowest five "
        f"live eigenvalues {spectrum} eV/(A^2 amu)"
    )
    record["minimum_source"] = str(path.relative_to(ROOT))
    record["minimum_energy_eV"] = energy
    record["basin_rmsd_A"] = config.basin_rmsd_A

    (OUT / "record.json").write_text(json.dumps(record, indent=1, default=str))
    (OUT / "certificate_DRAFT.txt").write_text(certificate_text(record))
    print(certificate_text(record))
    print(f"\nrecord: {OUT / 'record.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
