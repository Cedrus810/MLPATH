"""b05: measure the curvature magnitude of fragment relative TRANSLATION (PLAN item 5.3).

The three thresholds `fragment_translation_*` must come from a measurement on the
structures that motivated them -- the three order-1 elimination endpoints whose
residual imaginary mode is coupled inter-fragment motion -- not converted from the
radian set (`polish_soft_modes` footnote [1]: two unit systems through one threshold
mechanism is the silent kind of error).

What this script measures, per endpoint: displace the fragments rigidly along their
COM-connecting line (both signs, several offsets; the operator is
`internal.translate_rigid`), compute the energy at each, and report the local
second derivative d²E/dΔ² near Δ=0 in eV/A² -- the magnitude the curvature floor and
cap must bracket -- together with the gradient |dE/dΔ| in eV/A, which the tolerance
must be able to see. From those numbers the thresholds are SET, not guessed.

Needs the MACE checkpoint (POLAR-1): run on kasuga02, not yayoi. Usage:

    python docs/experiments/b05_fragment_translation_measurement.py

Writes runs/b05_fragment_translation_measurement/summary.json plus this stdout table.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from prrs import internal, torch_guard  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

ENDPOINTS = ROOT / "runs" / "b05_quench_continuation"
OUT = ROOT / "runs" / "b05_fragment_translation_measurement"
OFFSETS_A = np.array([-0.4, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.4])


def fragments_of(atoms, bond_scale=1.2):
    """Fragment membership by connectivity at the endpoint's own geometry."""
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix
    from ase.data import covalent_radii

    n = len(atoms)
    radii = covalent_radii[atoms.numbers]
    distances = atoms.get_all_distances()
    rows, cols = [], []
    for i in range(n):
        for j in range(i + 1, n):
            if distances[i, j] <= bond_scale * (radii[i] + radii[j]):
                rows += [i, j]
                cols += [j, i]
    adjacency = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    _, labels = connected_components(adjacency)
    groups = {}
    for atom, label in enumerate(labels):
        groups.setdefault(int(label), []).append(atom)
    return [np.asarray(members) for members in groups.values() if len(members) > 1], labels


def main():
    spec = json.loads((ROOT / "benchmarks" / "b05_tertiary_halide" / "model.json").read_text())
    from mace.calculators import MACECalculator

    torch_guard.install(trusted=[spec["path"]])

    def factory():
        return MACECalculator(
            model_paths=spec["path"],
            model_type=spec["mace_model_type"],
            device="cuda",
            default_dtype="float64",
        )

    config = SearchConfig(
        **json.loads(
            "\n".join(
                line
                for line in (ROOT / "benchmarks" / "b05_tertiary_halide" / "config.json")
                .read_text()
                .splitlines()
                if not line.lstrip().startswith("//")
            )
        )
    )

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for trial in ("t000002", "t000058", "t000080"):
        atoms = read(str(ENDPOINTS / f"{trial}_endpoint.extxyz"))
        atoms.calc = factory()
        fragments, _ = fragments_of(atoms, config.bond_scale)
        if len(fragments) < 2:
            print(f"{trial}: fewer than two multi-atom fragments; skipped")
            continue
        # decision 5: the fragment-pair COM line. Biggest fragment pair by mass.
        masses = atoms.get_masses()
        pair = max(
            ((a, b) for a in fragments for b in fragments if not np.array_equal(a, b)),
            key=lambda ab: masses[ab[0]].sum() * masses[ab[1]].sum(),
        )
        left, right = pair
        line = (
            atoms[right].get_positions().mean(axis=0) - atoms[left].get_positions().mean(axis=0)
        )
        energies = []
        for offset in OFFSETS_A:
            moved = atoms.copy()
            moved.set_positions(
                internal.translate_rigid(
                    moved.get_positions(),
                    list(left),
                    line,
                    float(offset),
                    masses,
                )
            )
            moved.calc = factory()
            energies.append(float(moved.get_potential_energy()))
        energies = np.asarray(energies)
        step = 0.05  # the tightest sampled spacing, around Delta=0
        curvature = float((energies[5] - 2 * energies[6] + energies[7]) / step**2)
        gradient = float((energies[7] - energies[5]) / (2 * step))
        rows.append(
            {
                "trial": trial,
                "fragments": [int(len(left)), int(len(right))],
                "curvature_eV_A2": curvature,
                "gradient_eV_A": gradient,
                "energies_eV": energies.tolist(),
                "offsets_A": OFFSETS_A.tolist(),
            }
        )
        print(
            f"{trial}: fragments {len(left)}+{len(right)}  "
            f"d2E/dD2 = {curvature:8.3f} eV/A^2   |dE/dD| = {gradient:8.4f} eV/A"
        )

    (OUT / "summary.json").write_text(
        json.dumps({"rows": rows, "direction": "fragment-pair COM line (decision 5)"}, indent=2),
        encoding="utf-8",
    )
    print(
        "\nSet the three thresholds FROM THESE NUMBERS (PLAN 5.3): the floor below the\n"
        "smallest |curvature| that is still a real well, the cap above the largest, the\n"
        "tolerance at the gradient scale a translation of that size produces. RATIONALE\n"
        "entries cite this run by directory."
    )


if __name__ == "__main__":
    main()
