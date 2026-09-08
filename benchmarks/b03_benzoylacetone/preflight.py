"""b03 preflight: the geometric reachability gate the P3 probes said must exist.

`docs/P3_BENZOYLACETONE_PROBES.md` section 1.1 records the most dangerous error of that
round. A probe built on the E isomer returned a complete, self-consistent, all-green set
of numbers: both endpoints real minima, different `chemical_key`, different `graph_hash`,
a clean single-proton `event_key`, dE = -128 meV. Every one of those is correct. The
reaction is still impossible, because the E isomer's hydroxyl points away from the
carbonyl and r(O...O) never drops below 4.33 A in 300 ETKDG conformers.

The identity layer does not stop it, and the reason matters. It is NOT that E and Z share
a graph -- `locked_bond_parity` records the C=C configuration and the two keys really do
differ (Z `8710064d12ce6a5a`, E `909db28462d72d1b`). It is that **the identity layer does
not answer reachability**. On the E graph the A->B proton transfer is a perfectly legal
cross-node reaction: bond delta, event key and direction quotient all hold. Nothing in
that layer can see "r(O...O) = 4.29 A, so the proton cannot get across".

So this gate is not a second identity check. It answers a different question, and it has
to exist separately for exactly that reason.

    python benchmarks/b03_benzoylacetone/preflight.py <source.extxyz>
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402

# The chelation criterion from the probes: r(O...O) < 3.0 A is "closed". Z spans
# 2.615-4.223 with 151/300 closed; E spans 4.334-4.667 with 0/300.
CLOSED_A = 3.0
# A source structure for this benchmark is the chelated enol, so it must already be closed
# rather than merely capable of closing.
SOURCE_MAX_A = 3.0
CONFORMERS = 300


def oxygen_pair(atoms):
    index = [i for i, z in enumerate(atoms.numbers) if z == 8]
    if len(index) != 2:
        raise ValueError(f"expected exactly two oxygens, found {len(index)}")
    return tuple(index)


def gate_source(atoms):
    """The supplied source must already be chelated."""
    i, j = oxygen_pair(atoms)
    distance = float(atoms.get_distance(i, j))
    return distance <= SOURCE_MAX_A, distance


def gate_isomer(smiles):
    """Can this constitution ever chelate at all? The E/Z question, asked geometrically.

    Conformers rather than the bond graph, because the graph is exactly what cannot answer
    it. rdkit is used here and nowhere in `prrs` -- this is a preflight, not a search path.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 0xF00D
    AllChem.EmbedMultipleConfs(molecule, numConfs=CONFORMERS, params=parameters)
    oxygens = [a.GetIdx() for a in molecule.GetAtoms() if a.GetSymbol() == "O"]
    if len(oxygens) != 2:
        raise ValueError("expected exactly two oxygens")
    distances = []
    for conformer in molecule.GetConformers():
        p = conformer.GetPositions()
        distances.append(float(np.linalg.norm(p[oxygens[0]] - p[oxygens[1]])))
    closed = sum(1 for d in distances if d < CLOSED_A)
    return closed, len(distances), min(distances), max(distances)


def main(argv):
    print(f"b03 preflight  ·  chelation reachability  ·  r(O...O) < {CLOSED_A} A = closed")
    print("=" * 78)
    failed = []

    for label, smiles in (
        ("Z (required)", r"C\C(=O)\C=C(/O)c1ccccc1"),
        ("E (must fail)", r"C\C(=O)\C=C(\O)c1ccccc1"),
    ):
        closed, total, low, high = gate_isomer(smiles)
        wanted = closed > 0 if label.startswith("Z") else closed == 0
        print(
            f"  [{'PASS' if wanted else 'FAIL'}] {label:14s} {closed:3d}/{total} closed   "
            f"r(O...O) {low:.3f} - {high:.3f} A"
        )
        if not wanted:
            failed.append(label)
    print(
        "         The E isomer returning 0/300 is the point: it is a complete, "
        "self-consistent,\n         all-green benchmark that is chemically impossible, "
        "and only geometry says so.\n"
    )

    for path in argv:
        atoms = read(path)
        ok, distance = gate_source(atoms)
        print(
            f"  [{'PASS' if ok else 'FAIL'}] source {Path(path).name:28s} "
            f"r(O...O) = {distance:.4f} A"
            + ("" if ok else f"  > {SOURCE_MAX_A} A, not chelated")
        )
        if not ok:
            failed.append(path)

    print("=" * 78)
    if failed:
        print("blocking:", ", ".join(str(f) for f in failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
