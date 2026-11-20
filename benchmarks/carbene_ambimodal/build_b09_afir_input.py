"""Unbiased b09 input for the Non-Hessian AFIR runner (docs/B09_DISTANCE_2026-11-20.md §5 A).

The older build_b09_input.py starts from the complex below TS-2OMe, i.e. it is built from the
answer. Here the two reactants are the SI's own reactant structures (lowest-G dimethoxycarbene
and 1-cyclohexene-1-carboxaldehyde, already re-optimised on POLAR-1 on 2026-10-08,
results/2026-10-08/opt/<name>/geometry.extxyz); no saddle geometry enters the positions.

The POLAR-1 TS-2OMe is used only as an atom-order template, so that classify_b09.py and the
reference RMSDs compare atom i with atom i: each reactant is mapped onto its fragment of the
template by an element-labelled bond-graph isomorphism.

The two relaxed fragments are placed with explore.poses.place along a seeded random direction
and orientation (contact ratio PLACE_SHELL), and the complex is relaxed on float32 POLAR-1. The
runner draws its own random poses around these fragments; the measured contact ratio of the
relaxed complex is printed as the --shell value (same rule as b05: measured contact ratio).

  python build_b09_afir_input.py OUT_DIR

Writes OUT_DIR/input.extxyz and OUT_DIR/build.json (energies, contact ratio, and per-fragment
heavy-atom RMSD of each reactant conformer to its fragment in the four POLAR-1 saddles -- a
diagnostic of which saddle the input conformers resemble; nothing is decided from it).
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.io import read, write
from ase.optimize import BFGS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.poses import place  # noqa: E402
from prrs.state import centred_kabsch_rmsd  # noqa: E402

OPT = HERE / "results/2026-10-08/opt"
TEMPLATE = OPT / "TS-2OMe/geometry.extxyz"
REACTANTS = {
    "C3H6O2": OPT / "Dimethoxycarbene/geometry.extxyz",
    "C7H10O": OPT / "1-Cyclohexene-1-carboxaldehyde_6O/geometry.extxyz",
}
SADDLES = ("TS-2OMe", "TS-2OMe-cis", "TS-2OMe-C2", "TS-2OMe-2")
BOND_SCALE = 1.2  # same connectivity scale as the runner's bond graphs
SEED = 20261120
PLACE_SHELL = 2.5
FMAX = 0.01


def graph(atoms):
    d = atoms.get_all_distances()
    z = atoms.numbers
    g = nx.Graph()
    for i, zi in enumerate(z):
        g.add_node(i, z=int(zi))
    for i in range(len(z)):
        for j in range(i + 1, len(z)):
            if d[i, j] < BOND_SCALE * (covalent_radii[z[i]] + covalent_radii[z[j]]):
                g.add_edge(i, j)
    return g


def formula(atoms, idx):
    return Atoms(numbers=atoms.numbers[idx]).get_chemical_formula(mode="hill")


def heavy_rmsd(x_atoms, y_atoms):
    heavy = [i for i, z in enumerate(x_atoms.numbers) if z > 1]
    x = x_atoms.positions[heavy] - x_atoms.positions[heavy].mean(axis=0)
    best = np.inf
    for sign in (1.0, -1.0):
        y = sign * y_atoms.positions[heavy]
        best = min(best, centred_kabsch_rmsd(x, y - y.mean(axis=0)))
    return float(best)


def relax(atoms, calc):
    a = atoms.copy()
    a.info.update({"charge": 0, "spin": 1})
    a.calc = calc
    opt = BFGS(a, logfile=None)
    ok = opt.run(fmax=FMAX, steps=2000)
    return a, bool(ok), int(opt.nsteps)


def contact_ratio(atoms, block_a, block_b):
    d = atoms.get_all_distances()
    z = atoms.numbers
    return float(
        min(d[i, j] / (covalent_radii[z[i]] + covalent_radii[z[j]]) for i in block_a for j in block_b)
    )


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    tmpl = read(TEMPLATE)
    blocks = [sorted(c) for c in nx.connected_components(graph(tmpl))]
    if len(blocks) != 2:
        raise RuntimeError(f"TS-2OMe template splits into {len(blocks)} fragments, expected 2")
    by_formula = {formula(tmpl, b): b for b in blocks}
    if set(by_formula) != set(REACTANTS):
        raise RuntimeError(f"template fragments {sorted(by_formula)} != {sorted(REACTANTS)}")
    rec = {"seed": SEED, "fmax_eV_A": FMAX, "place_shell": PLACE_SHELL, "template": str(TEMPLATE)}

    frag_pos = {}
    for key, path in REACTANTS.items():
        block = by_formula[key]
        sub = tmpl[block]
        react = read(path)
        gm = nx.algorithms.isomorphism.GraphMatcher(
            graph(sub), graph(react), node_match=lambda a, b: a["z"] == b["z"]
        )
        if not gm.is_isomorphic():
            raise RuntimeError(f"{key}: reactant graph is not isomorphic to the template fragment")
        m = gm.mapping  # template-fragment index -> reactant index
        ordered = Atoms(
            numbers=sub.numbers, positions=react.positions[[m[i] for i in range(len(sub))]]
        )
        relaxed, ok, n = relax(ordered, calc)
        frag_pos[key] = relaxed.positions.copy()
        rmsd = {}
        for s in SADDLES:
            ts = read(OPT / s / "geometry.extxyz")
            rmsd[s] = heavy_rmsd(relaxed, ts[block])
        rec[key] = {
            "source": str(path),
            "template_atoms": block,
            "energy_eV": float(relaxed.get_potential_energy()),
            "converged": ok,
            "steps": n,
            "heavy_rmsd_to_saddle_fragment_A": rmsd,
        }

    # complex in template order: place fragment B (aldehyde) around fragment A (carbene)
    key_a, key_b = "C3H6O2", "C7H10O"
    blk_a, blk_b = by_formula[key_a], by_formula[key_b]
    rng = np.random.default_rng(SEED)
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    quat = rng.normal(size=4)
    quat /= np.linalg.norm(quat)
    new_b, ratio, status, _ = place(
        frag_pos[key_a], tmpl.numbers[blk_a], frag_pos[key_b], tmpl.numbers[blk_b],
        direction, quat, PLACE_SHELL,
    )
    if status != "placed":
        raise RuntimeError(f"placement failed: {status}")
    pos = np.zeros((len(tmpl), 3))
    pos[blk_a] = frag_pos[key_a]
    pos[blk_b] = new_b
    cplx, ok, n = relax(Atoms(numbers=tmpl.numbers, positions=pos), calc)
    e_sum = rec[key_a]["energy_eV"] + rec[key_b]["energy_eV"]
    cr = contact_ratio(cplx, blk_a, blk_b)
    rec["complex"] = {
        "energy_eV": float(cplx.get_potential_energy()),
        "converged": ok,
        "steps": n,
        "binding_eV": float(cplx.get_potential_energy()) - e_sum,
        "placed_contact_ratio": float(ratio),
        "contact_ratio": cr,
        "suggested_shell": round(cr, 1),
        "n_fragments_after_relax": nx.number_connected_components(graph(cplx)),
    }
    cplx.info.update({"charge": 0, "spin": 1})
    write(out / "input.extxyz", cplx.copy(), format="extxyz")
    (out / "build.json").write_text(json.dumps(rec, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
