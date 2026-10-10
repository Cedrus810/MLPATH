"""Where did b09's reaction directions rank, before and after scoring the delivered motion?

The b09 probe (runs/b09_probe1) never selected compress/kick 2-15 or 2-13 -- the carbene C
approaching the two ring alkene carbons -- so it never pushed along the bonds that form.
The direction budget is ranked by overlap with the four softest non-torsional modes. Two
questions, both on the run's own source microstate and its own config:

1. Scoring. The old score used the endpoint gradient g/sqrt(m); a displacement family is
   delivered as a rigid motion of whole fragments (perturbations.apply). Print each
   family's ranking under both scores and which directions each would select.
2. The soft modes themselves. soft_internal_modes projects torsions out and renormalises,
   keeping the spectral order; print how much of each kept mode survived the projection
   (a mode that is nearly all torsion is amplified by the renormalisation).

Usage: python benchmarks/carbene_ambimodal/rank_b09.py <source extxyz> <config json> <network.json>
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs import internal, perturbations  # noqa: E402
from prrs.canonical import canonical_form  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.perturbations import FAMILIES, propose, rotatable_torsions  # noqa: E402
from prrs.runner import _trivial_modes, curvature_spectrum, soft_internal_modes  # noqa: E402
from prrs.state import encode  # noqa: E402

WATCH = ("compress:2-15:-", "compress:2-13:-", "kick:2-15:+", "kick:2-15:-", "kick:2-13:+")


def old_overlap(atoms, family, indices, modes, bond_scale=1.2, active=None, graph=None):
    gradient = internal.gradient(atoms.positions, FAMILIES[family]["kind"], indices)
    vector = (gradient / np.sqrt(atoms.get_masses())[:, None]).reshape(-1)
    vector = vector / np.linalg.norm(vector)
    return max(abs(float(vector @ mode)) for _, mode in modes)


def survival(atoms, factory, config):
    """Fraction of each soft mode left after the torsional projection, in spectral order."""
    eigenvalues, vectors, _, _ = curvature_spectrum(
        atoms, factory, config, source=config.hessian_source
    )
    trivial = _trivial_modes(atoms.positions, atoms.get_masses())
    graph = encode(atoms, config.bond_scale, active=config.active_atoms)
    form = canonical_form(atoms.numbers, graph.edges, graph.index)
    orbits = form.orbits if form.info.get("enumerated") else form.colours
    genuine, rotors = rotatable_torsions(graph, form.colours, orbits)
    columns = []
    for torsion in list(genuine) + list(rotors):
        g = internal.gradient(atoms.positions, "dihedral", torsion.indices)
        v = (g / np.sqrt(atoms.get_masses())[:, None]).reshape(-1)
        columns.append(v / np.linalg.norm(v))
    torsional = np.linalg.qr(np.array(columns).T)[0] if columns else None
    rows = []
    for column in range(vectors.shape[1]):
        if len(rows) >= config.spectral_rank_modes:
            break
        mode = vectors[:, column]
        if float(np.linalg.norm(trivial @ mode)) > 0.5:
            continue
        left = mode - torsional @ (torsional.T @ mode) if torsional is not None else mode
        rows.append((float(eigenvalues[column]), float(np.linalg.norm(left))))
    return rows


def ranking(atoms, config, modes, score):
    perturbations.spectral_overlap = score
    report = {}
    list(propose(atoms, config, config.seed, report=report, soft_modes=modes))
    scores = {}
    for identifier in report["candidate_ids"]:
        family, indices, _ = identifier.split(":")
        indices = tuple(int(i) for i in indices.split("-"))
        scores[identifier] = score(
            atoms, family, indices, modes, config.bond_scale, config.active_atoms
        )
    return report["selected_ids"], scores


def main():
    atoms = read(sys.argv[1])
    raw = json.loads(Path(sys.argv[2]).read_text())
    config = SearchConfig.from_dict(raw)
    network = json.loads(Path(sys.argv[3]).read_text())
    atoms.info.update(charge=config.total_charge, spin=config.multiplicity)
    calc = polar1()
    atoms.calc = calc
    factory = lambda: calc  # noqa: E731
    modes, note = soft_internal_modes(atoms, factory, config)
    print(f"soft modes: {len(modes)} ({note})")
    for eigenvalue, left in survival(atoms, factory, config):
        print(f"  eigenvalue {eigenvalue:+.5f}  non-torsional fraction kept {left:.3f}")

    new_score = perturbations.spectral_overlap
    old_sel, old = ranking(atoms, config, modes, old_overlap)
    new_sel, new = ranking(atoms, config, modes, new_score)
    perturbations.spectral_overlap = new_score
    run_sel = set(network["direction_ledger"]["executed_ids"]) | set(
        network["direction_ledger"]["not_executed_ids"]
    )
    print(f"\nrun selected {len(run_sel)}; old score reproduces it: {set(old_sel) == run_sel}")
    print(f"old selection ∩ run {len(set(old_sel) & run_sel)}/{len(run_sel)}")
    for name, scores, sel in (("old", old, old_sel), ("new", new, new_sel)):
        print(f"\n{name} score: selected {len(sel)}")
        for family in ("compress", "kick", "stretch", "bend"):
            ids = sorted(
                (k for k in scores if k.startswith(family + ":")), key=lambda k: -scores[k]
            )
            picked = [k for k in sel if k.startswith(family + ":")]
            print(f"  {family}: {len(ids)} candidates, {len(picked)} selected")
            for k in ids[:8]:
                print(
                    f"    {ids.index(k) + 1:3d} {k:22s} {scores[k]:.3f}{'  *' if k in sel else ''}"
                )
            for k in WATCH:
                if k.startswith(family + ":") and k in scores:
                    print(
                        f"    watch {k:22s} rank {ids.index(k) + 1}/{len(ids)} score {scores[k]:.3f}"
                        f"{'  selected' if k in sel else ''}"
                    )
    print(json.dumps({"old_selected": old_sel, "new_selected": new_sel}))


if __name__ == "__main__":
    main()
