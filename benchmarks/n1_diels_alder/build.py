"""N1 (cyclopentadiene + acrolein Diels-Alder): input complex and TS guesses
(docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §5 n0).

Reference side only: the guesses use the known product and are never shown to the search.

  python build.py OUT_DIR

Atom order everywhere (the D2 runs read input.extxyz): cyclopentadiene C1..C5 (C1, C4 the
diene termini, C5 the CH2), its 6 H, then acrolein C6 (CH2=), C7 (=CH-), C8 (CHO), O9, its
4 H; hydrogens follow their heavy atom in map order. Writes
  fragments_polar1.extxyz   the two fragments, each relaxed on float32 POLAR-1
  input.extxyz              face-to-face pre-reaction complex relaxed on float32 POLAR-1
  guess_<endo|exo>_<s-cis|s-trans>.xyz   MMFF guesses: forming C1-C6 / C4-C7 held at
                            GUESS_FORMING_A, O9=C8-C7=C6 dihedral held at 0 / 180 deg
  build.json                energies, contact ratio of the complex, endo/exo assignment
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from ase import Atoms  # noqa: E402
from ase.io import write  # noqa: E402
from ase.optimize import BFGS  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402

SEED = 20261109
FMAX = 0.01
GUESS_FORMING_A = 2.25
CPD = "[CH:1]1=[CH:2][CH:3]=[CH:4][CH2:5]1"
ACR = "[CH2:6]=[CH:7][CH:8]=[O:9]"
PRODUCT = "[CH:1]12[CH:2]=[CH:3][CH:4]([CH2:5]1)[CH:7]([CH:8]=[O:9])[CH2:6]2"
COVALENT = {1: 0.31, 6: 0.76, 8: 0.66}


BLOCKS = ((1, 2, 3, 4, 5), (6, 7, 8, 9))  # cyclopentadiene, acrolein
# indices in the module order: C1 0, C4 3 (diene termini); C6 11, C7 12 (dienophile)
FORMING = ((0, 11), (3, 12))


def ordered(mol):
    """Indices of mol in the module order: per block (cyclopentadiene, then acrolein) the
    heavy atoms by map number, then their hydrogens grouped by parent map number."""
    idx = []
    maps = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    for block in BLOCKS:
        present = [maps[m] for m in block if m in maps]
        hyd = []
        for i in present:
            hyd += sorted(
                n.GetIdx()
                for n in mol.GetAtomWithIdx(i).GetNeighbors()
                if n.GetAtomicNum() == 1
            )
        idx += present + hyd
    return idx


def embed(smiles, seed=SEED):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise RuntimeError(f"embed failed: {smiles}")
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def mmff_constrained(mol, distances=(), dihedrals=()):
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    for i, j, d in distances:
        ff.MMFFAddDistanceConstraint(i, j, False, d, d, 1.0e4)
    for i, j, k, m, t in dihedrals:
        ff.MMFFAddTorsionConstraint(i, j, k, m, False, t, t, 1.0e3)
    ff.Minimize(maxIts=20000)
    return mol


def by_map(mol):
    return {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}


def atoms_of(mol, conf_id=-1):
    idx = ordered(mol)
    if len(idx) != mol.GetNumAtoms():
        raise RuntimeError("ordering lost atoms")
    pos = mol.GetConformer(conf_id).GetPositions()
    return Atoms(
        numbers=[mol.GetAtomWithIdx(i).GetAtomicNum() for i in idx], positions=pos[idx]
    )


def relax(atoms, calc, charge=0):
    a = atoms.copy()
    a.info.update({"charge": charge, "spin": 1})
    a.calc = calc
    opt = BFGS(a, logfile=None)
    converged = opt.run(fmax=FMAX, steps=2000)
    return a, bool(converged), int(opt.nsteps)


def contact_ratio(atoms, n_first):
    """min inter-fragment distance / sum of covalent radii of that pair (count_l0 definition)."""
    best = None
    for i in range(n_first):
        for j in range(n_first, len(atoms)):
            d = atoms.get_distance(i, j)
            r = d / (COVALENT[atoms.numbers[i]] + COVALENT[atoms.numbers[j]])
            best = r if best is None or r < best else best
    return best


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    rec = {"seed": SEED, "fmax_eV_A": FMAX, "forming_A": GUESS_FORMING_A}

    # fragments: cyclopentadiene, s-trans acrolein (dihedral O9=C8-C7=C6 = 180)
    cpd = embed(CPD)
    acr = embed(ACR)
    m = by_map(acr)
    mmff_constrained(acr, dihedrals=[(m[9], m[8], m[7], m[6], 180.0)])
    frags = []
    for name, mol in (("cyclopentadiene", cpd), ("acrolein_s-trans", acr)):
        a, ok, n = relax(atoms_of(mol), calc)
        rec[name] = {"energy_eV": float(a.get_potential_energy()), "converged": ok, "steps": n}
        frags.append(a)
    write(out / "fragments_polar1.extxyz", [f.copy() for f in frags])  # copies: no shared calc

    # complex: acrolein C=C 3.4 A above the diene plane, centred over C1..C4
    cp, ac = frags[0].positions.copy(), frags[1].positions.copy()
    ring = cp[:5]
    centre = ring[:4].mean(axis=0)
    normal = np.cross(ring[1] - ring[0], ring[3] - ring[0])
    normal /= np.linalg.norm(normal)
    if np.dot(cp[4] - centre, normal) > 0:  # put acrolein on the face away from the CH2
        normal = -normal
    ac = ac - ac[[0, 1]].mean(axis=0) + centre + 3.4 * normal
    cplx = Atoms(
        numbers=list(frags[0].numbers) + list(frags[1].numbers), positions=np.vstack([cp, ac])
    )
    cplx, ok, n = relax(cplx, calc)
    rec["complex"] = {
        "energy_eV": float(cplx.get_potential_energy()),
        "converged": ok,
        "steps": n,
        "contact_ratio": contact_ratio(cplx, len(frags[0])),
        "binding_eV": float(cplx.get_potential_energy())
        - rec["cyclopentadiene"]["energy_eV"]
        - rec["acrolein_s-trans"]["energy_eV"],
    }
    cplx.info.update({"charge": 0, "spin": 1})
    write(out / "input.extxyz", cplx.copy(), format="extxyz")

    # TS guesses from the mapped product: endo / exo x s-cis / s-trans
    guesses = {}
    for iso in EnumerateStereoisomers(Chem.MolFromSmiles(PRODUCT)):
        mol = Chem.AddHs(iso)
        if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
            continue
        AllChem.MMFFOptimizeMolecule(mol)
        mp = by_map(mol)
        pos = mol.GetConformer().GetPositions()
        # endo: the CHO carbon C8 is nearer the C2=C3 bridge than the CH2 bridge C5
        d_alkene = np.linalg.norm(pos[mp[8]] - (pos[mp[2]] + pos[mp[3]]) / 2)
        d_ch2 = np.linalg.norm(pos[mp[8]] - pos[mp[5]])
        kind = "endo" if d_alkene < d_ch2 else "exo"
        for conf, torsion in (("s-cis", 0.0), ("s-trans", 180.0)):
            key = f"{kind}_{conf}"
            if key in guesses:
                continue
            g = Chem.Mol(mol)
            mmff_constrained(
                g,
                distances=[
                    (mp[1], mp[6], GUESS_FORMING_A),
                    (mp[4], mp[7], GUESS_FORMING_A),
                ],
                dihedrals=[(mp[9], mp[8], mp[7], mp[6], torsion)],
            )
            a = atoms_of(g)
            if list(a.numbers) != list(cplx.numbers):
                raise RuntimeError(
                    f"{key}: atom order {list(a.numbers)} != {list(cplx.numbers)}"
                )
            a.info.update({"charge": 0, "spin": 1})
            write(out / f"guess_{key}.xyz", a, format="extxyz")
            guesses[key] = {
                "chirality": Chem.MolToSmiles(iso),
                "d_C8_to_alkene_A": float(d_alkene),
                "d_C8_to_CH2_A": float(d_ch2),
                "forming_A": [float(a.get_distance(i, j)) for i, j in FORMING],
                # module order: C6 11, C7 12, C8 13, O9 14
                "dihedral_O9C8C7C6_deg": float(a.get_dihedral(14, 13, 12, 11)),
            }
    rec["guesses"] = guesses
    (out / "build.json").write_text(json.dumps(rec, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
