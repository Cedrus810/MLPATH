"""Does a spectral ranking of proposal directions beat the one in the code?

Idea A of the 2026-09-13 brainstorm: every admitted minimum already pays for an analytic
mass-weighted Hessian (`minimum_check = "hessian"` is the default, and `curvature_spectrum`
diagonalises it), and `confirm_minimum` keeps six eigenvalues and throws every eigenvector
away. The proposal layer then orders directions by a seeded shuffle inside each family.
The proposal is to rank instead by overlap with the softest modes that are NOT rotors --
soft because a soft direction is where large-amplitude motion is cheap, not-rotor because
the softest modes at any minimum are torsions by construction and a torsion reaches a
conformer, not a reaction.

Nothing is changed here. This measures whether the gain exists before anything is built,
because the project's own rule is that a knob is verified before an experiment is designed
around it (docs/HANDOFF.md section 2), and idea A has no failure behind it -- it has a
plausibility argument. The measurement it needs is narrow:

    on a run where the answer is known, where does the direction that actually
    reacted land under the spectral score, versus where it landed in the run?

b04 is the case with a known answer: 19 executed directions, exactly one of which moved
the system (`stretch (0,1) +1`, the C-Cl bond). It came out first in the execution order,
but there is no mechanism that put it there -- with no hydrogen-bond donor in Cl- + CH3Cl
every candidate sits in priority band 1, so the order inside a family is the seeded
shuffle. So the question is not "did the run fail" (it did not) but "was it luck", and
the answer decides whether a spectral ranking converts a lottery into a first pick.

    python docs/experiments/spectral_proposal_probe.py [run_dir]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from prrs import internal  # noqa: E402
from prrs.chemistry import canonical_labels  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.perturbations import (  # noqa: E402
    FAMILIES,
    _angle_directions,
    _bond_directions,
    _dihedral_directions,
    candidate_id,
    rotatable_torsions,
)
from prrs.runner import _trivial_modes, curvature_spectrum, floor_residual  # noqa: E402
from prrs.state import encode  # noqa: E402

SOFT_MODES = 4


def mass_weighted_tangent(atoms, kind, indices):
    """The direction `displace` actually moves in, expressed in mass-weighted coordinates.

    `displace` steps along g/m in Cartesian coordinates -- the internal coordinate's own
    gradient over the masses -- so in mass-weighted coordinates q = sqrt(m) x that is
    g/sqrt(m). Using g or g*sqrt(m) instead is the same error the relay mass-weighting
    made once already (M^(1/2) g where M^-1 g was meant); the overlap it computes would
    be with a direction nothing travels in.
    """
    gradient = internal.gradient(atoms.positions, kind, indices)
    vector = (gradient / np.sqrt(atoms.get_masses())[:, None]).reshape(-1)
    norm = np.linalg.norm(vector)
    return vector / norm if norm > 1e-30 else None


def rotor_subspace(atoms, graph, labels):
    """An orthonormal basis for the torsional directions, mass-weighted.

    Both halves of `rotatable_torsions`: the genuine torsions and the symmetry rotors. A
    genuine torsion is as much a conformer coordinate as a methyl spin is, and the failure
    this projection exists to avoid -- eight rotor barriers at -65 to -80 cm^-1 while the
    proton transfer sits at -3186 -- does not care which half a torsion came from.
    """
    genuine, rotors = rotatable_torsions(graph, labels)
    columns = []
    for torsion in list(genuine) + list(rotors):
        vector = mass_weighted_tangent(atoms, "dihedral", torsion.indices)
        if vector is not None:
            columns.append(vector)
    if not columns:
        return np.zeros((len(atoms) * 3, 0)), 0
    basis, _ = np.linalg.qr(np.array(columns).T)
    return basis, len(columns)


def spectral_scores(atoms, config, factory):
    """Score every candidate direction by overlap with the softest non-rotor modes."""
    eigenvalues, vectors, floor, provenance = curvature_spectrum(
        atoms, factory, config, source=config.hessian_source
    )
    graph = encode(atoms, config.bond_scale, active=config.active_atoms)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    rotors, rotor_count = rotor_subspace(atoms, graph, labels)

    # Two subspaces come out first, and forgetting the first one is how the first draft of
    # this probe scored every direction at exactly zero. `curvature_spectrum` PROJECTS the
    # rigid-body modes out rather than deleting them, so they stay in the spectrum as a
    # six-dimensional null block at the soft end -- "the four softest modes" are the
    # translations and rotations, and they are orthogonal to every internal coordinate by
    # construction. They are identified the way runner.py already does it, by how much of
    # the eigenvector the trivial basis accounts for.
    trivial = _trivial_modes(atoms.positions, atoms.get_masses())
    retained = []
    for column in range(vectors.shape[1]):
        mode = vectors[:, column]
        if float(np.linalg.norm(trivial @ mode)) > 0.5:
            continue  # a rigid-body mode wearing a zero eigenvalue
        if rotors.shape[1]:
            mode = mode - rotors @ (rotors.T @ mode)
        survival = float(np.linalg.norm(mode))
        if survival < 1e-8:
            continue  # entirely rotor: the direction this projection exists to drop
        retained.append((float(eigenvalues[column]), mode / survival, survival))
    retained = retained[:SOFT_MODES]

    families = tuple(config.families)
    directions = (
        list(_bond_directions(atoms, graph, config, families))
        + list(_angle_directions(graph, labels, families))
        + list(_dihedral_directions(graph, labels, families))
    )
    scored = {}
    for entry in directions:
        family, indices, sign = entry[0], entry[1], entry[2]
        kind = FAMILIES[family]["kind"]
        tangent = mass_weighted_tangent(atoms, kind, indices)
        if tangent is None:
            continue
        overlap = max(abs(float(tangent @ mode)) for _, mode, _ in retained) if retained else 0.0
        scored[candidate_id(family, indices, sign)] = overlap
    return scored, retained, eigenvalues, rotor_count, provenance, floor


def main(argv):
    run = Path(argv[1] if len(argv) > 1 else "runs/b04_probe4_ledger")
    network = json.loads((run / "network.json").read_text())
    settings = json.loads((run / "config.json").read_text())
    config = SearchConfig.from_dict(settings["config"])

    from ase.io import read

    source_id = network["amplitude_scans"][0]["source_microstate"]
    atoms = read(str(run / "chemical" / f"{source_id.replace('/', '_')}.extxyz"))
    if config.charge_sensitive:
        atoms.info["charge"] = config.total_charge
        atoms.info["spin"] = config.multiplicity

    module, _, name = settings["calculator"]["factory"].partition(":")
    sys.path.insert(0, str(ROOT))
    factory = getattr(__import__(module, fromlist=[name]), name)

    scored, retained, eigenvalues, rotor_count, provenance, floor = spectral_scores(
        atoms, config, factory
    )

    print(f"run            {run}")
    print(f"source         {source_id}  ({len(atoms)} atoms, charge {atoms.info.get('charge')})")
    print(f"hessian        {provenance}   trivial-mode floor worst {floor_residual(floor):.2e}")
    print(f"rotor basis    {rotor_count} torsional direction(s) projected out")
    print(f"eigenvalues    {[round(float(v), 4) for v in eigenvalues[:8]]}  eV/(A^2 amu)")
    print(f"retained modes {[(round(v, 4), round(s, 3)) for v, _, s in retained]}  (eigenvalue, rotor-survival)")
    print()

    executed = [
        candidate_id(scan["family"], scan["indices"], scan["sign"])
        for scan in network["amplitude_scans"]
    ]
    reacted = {
        candidate_id(scan["family"], scan["indices"], scan["sign"])
        for scan in network["amplitude_scans"]
        if scan.get("local_brackets")
    }
    ranking = sorted(scored, key=lambda i: -scored[i])

    print("spectral ranking of every candidate (overlap with softest non-rotor modes):")
    for place, identifier in enumerate(ranking, start=1):
        run_place = executed.index(identifier) + 1 if identifier in executed else None
        mark = "  <== REACTED" if identifier in reacted else ""
        print(
            f"  {place:3d}  {identifier:22} overlap {scored[identifier]:.4f}   "
            f"run order {run_place if run_place else '-- not executed'}{mark}"
        )

    print()
    for identifier in sorted(reacted):
        place = ranking.index(identifier) + 1 if identifier in ranking else None
        print(
            f"VERDICT  the direction that reacted ({identifier}) ranks {place} of "
            f"{len(ranking)} under the spectral score, and ran "
            f"{executed.index(identifier) + 1} of {len(executed)} in the run."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
