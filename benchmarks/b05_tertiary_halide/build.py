"""Build the b05 input: tert-butyl chloride and a hydroxide ion, one fragment apart.

Deterministic and offline. The geometry is a STARTING POINT, not an answer: the search
relaxes it and decides for itself what is a basin.

The hard part of this grid is not the substrate, it is where to put the base. Two
channels compete for it -- attack at the carbon, or abstraction of a beta hydrogen --
and any placement chosen to favour one of them puts the mechanism into the premise,
which is what benchmarks/README.md rule 5 forbids. So the placement rule is stated
mechanically, in terms of the proposal layer's reach and of this project's own
hydrogen-bond threshold, and never in terms of which channel it helps:

    1. on the C3 axis behind the C-Cl bond, so the three methyls stay equivalent and
       no face is favoured
    2. the oxygen outside bonding range of every substrate atom, with the bond graph
       reporting two fragments -- it must be a separated pair, not a bound adduct
    3. inside pair_cutoff_A of the central carbon AND of at least one beta hydrogen,
       so both channels can be proposed

That lands on b04's own 3.6 A and b04's own cutoff of 4.0, and NOTHING here is tuned.
Measured: O...C 3.600, the three anti-periplanar hydrogens at 2.463 each, two fragments,
nearest contact 2.463 A against an O-H bonding threshold of 1.16 A.

Two things about that 2.463 A, because both were got wrong in drafts of this file:

    It is not a pre-formed hydrogen bond. config's hbond_acceptor_A = 2.5 covers X-H
    donors with X in (N, O, F) -- perturbations.ACCEPTORS -- and a methyl C-H is not a
    donor there, so that threshold does not describe this contact and a "clearance"
    rule built on it was measuring nothing.

    The source relaxing into a contact complex is not a flaw to design out. b04 starts
    at 3.600 A and its source minimum relaxes to 3.1835 A -- its ion-dipole complex.
    The starting point is a starting point; what the quench makes of it is recorded.

The beta hydrogens do fill the hemisphere around a tertiary centre, which is the same
fact textbooks state as "SN2 is hindered here": sweeping every direction at every
radius, the best any point achieves is min(O...H) = O...C - 0.99 A. One draft answered
that by raising pair_cutoff_A to 4.6 and moving out to 4.30 A. THAT MADE THINGS WORSE,
measurably: the larger cutoff took the non-bonded candidate count from 25 to 37 while
the direction budget still selected 7, and the two pairs this grid depends on were
crowded out -- preflight gate 6 went from PASS to FAIL. Reach cannot be bought with the
cutoff; the direction budget is what binds. The cutoff stays at b04's 4.0.

Both channels are then proposable and the search picks, which is the only arrangement
under which "PRRS did not find the elimination" means anything at all. main() prints
the measured distances so the rule is checkable rather than asserted.

The staggered rotamer is used because it is tert-butyl chloride's rotamer minimum, not
because of what it does or does not make anti-periplanar. What that choice implies for
any particular pathway is a measurement (printed below), not a design goal.

Run:  python benchmarks/b05_tertiary_halide/build.py [out.extxyz]
"""

import sys
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

# main() checks rule 2 with the project's own bond criterion, so prrs has to be
# importable when this file is run as a script -- the same line preflight.py carries. It
# was missing, so `build.py out.extxyz` raised ModuleNotFoundError at the fragment check
# and exited 1 BEFORE writing, while still printing every measurement above it. The stored
# input stayed at an earlier draft's 4.30 A for two probe runs because of it, and the
# printout was read as proof the file had been written. Gate 10 now checks the artifact
# instead of the printout.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

TOTAL_CHARGE = -1
MULTIPLICITY = 1

# Standard values, not optimised here: the search relaxes the structure before it
# perturbs it, and a starting geometry that pretends to be an optimum would be a claim
# about the model's minimum rather than an input.
R_C_CL = 1.80  # Angstrom, tertiary C-Cl
R_C_C = 1.52
R_C_H = 1.09
R_O_H = 0.964  # hydroxide
TETRAHEDRAL = np.deg2rad(109.471)

APPROACH_A = 3.6  # oxygen to central carbon; b04's value, see the docstring
PAIR_CUTOFF_A = 4.0  # b04's value; config.json must agree
# Rule 2 is a bonding question, not a hydrogen-bond one: the O-H covalent threshold
# is bond_scale * (r_O + r_H) = 1.2 * (0.66 + 0.31) = 1.16 A.
BONDING_RANGE_A = 1.16


def _tert_butyl_chloride():
    """C(CH3)3Cl with Cl on +z, the three methyls staggered below.

    Built from bond lengths and angles rather than from a stored coordinate list so the
    construction is auditable: every number above appears once.
    """
    positions = [np.zeros(3)]  # 0: central carbon
    symbols = ["C"]
    positions.append(np.array([0.0, 0.0, R_C_CL]))  # 1: chlorine
    symbols.append("Cl")

    # Three methyl carbons on a cone below the carbon, tetrahedral to the C-Cl bond.
    # The polar angle IS the tetrahedral angle, measured from the C-Cl direction: cos is
    # -1/3, so the methyls sit on the far side from the chlorine. Writing pi - TETRAHEDRAL
    # here instead puts them at 70.5 degrees, on the SAME side, and Cl...C(methyl) comes
    # out at 1.93 A -- inside the bond cutoff, so the graph grows three bonds that are not
    # there and the chlorine ends up with degree 4. The first draft did exactly that and
    # preflight gate 7 caught it: the identity layer reported degree_exceeds_valence on
    # the reactants. Bond lengths and angles alone would not have; the graph is what
    # noticed.
    polar = TETRAHEDRAL
    for k in range(3):
        azimuth = 2 * np.pi * k / 3
        direction = np.array(
            [
                np.sin(polar) * np.cos(azimuth),
                np.sin(polar) * np.sin(azimuth),
                np.cos(polar),
            ]
        )
        methyl_c = direction * R_C_C
        positions.append(methyl_c)
        symbols.append("C")
        # Local frame on this methyl: its own C-C bond is the axis.
        axis = direction
        reference = np.array([0.0, 0.0, 1.0])
        perpendicular = reference - np.dot(reference, axis) * axis
        perpendicular /= np.linalg.norm(perpendicular)
        second = np.cross(axis, perpendicular)
        # Staggered. The 60 degree offset is the whole point: without it the hydrogens
        # land at Cl-C-C-H = 0, 120, 240 degrees, which is the ECLIPSED rotamer -- the
        # rotamer maximum, and the first draft of this file built it while the docstring
        # claimed otherwise. main() prints the dihedrals so that cannot happen silently.
        for j in range(3):
            phase = 2 * np.pi * j / 3 + np.pi / 3
            hydrogen = methyl_c + R_C_H * (
                np.cos(np.pi - TETRAHEDRAL) * axis
                + np.sin(np.pi - TETRAHEDRAL)
                * (np.cos(phase) * perpendicular + np.sin(phase) * second)
            )
            positions.append(hydrogen)
            symbols.append("H")
    return Atoms(symbols, positions=np.array(positions))


def build(separation_A=APPROACH_A):
    atoms = _tert_butyl_chloride()
    # Backside of C-Cl, on the C3 axis, `separation_A` from the central carbon.
    atoms += Atoms("O", positions=[[0.0, 0.0, -float(separation_A)]])
    # The hydroxide's own hydrogen points away from the substrate, so the oxygen -- the
    # nucleophilic and basic atom -- is the one facing it. Pointing it at the substrate
    # would be pre-arranging a hydrogen bond.
    atoms += Atoms("H", positions=[[0.0, 0.0, -float(separation_A) - R_O_H]])
    atoms.info["charge"] = TOTAL_CHARGE
    atoms.info["spin"] = MULTIPLICITY
    atoms.pbc = False
    return atoms


CENTRAL_C, CHLORINE = 0, 1
METHYL_C = (2, 6, 10)
BETA_H = tuple(i for i in range(3, 14) if i not in METHYL_C)
OXYGEN, HYDROXIDE_H = 14, 15


def substitution_products():
    """tert-butanol plus a chloride, far apart. Two fragments.

    Crude on purpose: it exists so preflight can ask whether the identity layer could
    tell the candidate outcomes apart, and whether their automorphism groups enumerate.
    It is NOT a claim about what is on the surface. Lives here rather than in preflight
    because this grid's convention is that structures are built by build.py -- and
    because importing a gate script to get a structure runs every gate in it.
    """
    out = build()
    out.positions[CHLORINE] = out.positions[CENTRAL_C] + np.array([0.0, 0.0, 6.0])
    direction = out.positions[OXYGEN] - out.positions[CENTRAL_C]
    direction /= np.linalg.norm(direction)
    out.positions[OXYGEN] = out.positions[CENTRAL_C] + 1.43 * direction
    out.positions[HYDROXIDE_H] = out.positions[OXYGEN] + np.array([0.9, 0.0, -0.35])
    return out


def elimination_products():
    """Isobutylene, water and a chloride. Three fragments."""
    out = build()
    beta = min(BETA_H, key=lambda h: out.get_distance(OXYGEN, h))
    methyl = min(METHYL_C, key=lambda c: out.get_distance(c, beta))
    # the beta hydrogen joins the oxygen
    direction = out.positions[OXYGEN] - out.positions[beta]
    out.positions[beta] = out.positions[OXYGEN] - 0.96 * direction / np.linalg.norm(direction)
    # the alkene forms, the chloride leaves
    axis = out.positions[methyl] - out.positions[CENTRAL_C]
    out.positions[methyl] = out.positions[CENTRAL_C] + 1.33 * axis / np.linalg.norm(axis)
    out.positions[CHLORINE] = out.positions[CENTRAL_C] + np.array([0.0, 0.0, 7.0])
    return out


def main(argv):
    atoms = build()
    cutoff = PAIR_CUTOFF_A
    print(f"formula          {atoms.get_chemical_formula()}   atoms {len(atoms)}")
    print(f"total charge     {atoms.info['charge']}   multiplicity {atoms.info['spin']}")
    print(f"r(C-Cl)          {atoms.get_distance(CENTRAL_C, CHLORINE):.4f} A")
    print(f"r(O...C)         {atoms.get_distance(OXYGEN, CENTRAL_C):.4f} A")
    print(f"r(O-H) hydroxide {atoms.get_distance(OXYGEN, HYDROXIDE_H):.4f} A")
    print()
    print(f"the placement rule, measured (pair_cutoff_A = {cutoff}):")
    within_c = atoms.get_distance(OXYGEN, CENTRAL_C) < cutoff
    beta = sorted((atoms.get_distance(OXYGEN, h), h) for h in BETA_H)
    within_h = [(h, d) for d, h in beta if d < cutoff]
    print(
        f"  O within cutoff of the central carbon      {within_c}"
        f"  ({atoms.get_distance(OXYGEN, CENTRAL_C):.3f} A)"
    )
    print(f"  beta hydrogens within cutoff of O          {len(within_h)} of {len(BETA_H)}")
    for h, d in within_h[:6]:
        print(
            f"     H{h}  {d:.3f} A   on methyl carbon "
            f"{min(METHYL_C, key=lambda c: atoms.get_distance(c, h))}"
        )
    print(f"  nearest beta hydrogen                      {beta[0][1]} at {beta[0][0]:.3f} A")
    nearest_any = min(
        (atoms.get_distance(OXYGEN, i), i)
        for i in range(len(atoms))
        if i not in (OXYGEN, HYDROXIDE_H)
    )
    print(
        f"  nearest substrate atom of any kind         "
        f"{atoms.get_chemical_symbols()[nearest_any[1]]}{nearest_any[1]} "
        f"at {nearest_any[0]:.3f} A   (rule 2 wants > {BONDING_RANGE_A}, and 2 fragments)"
    )
    print()
    # Anti-periplanar geometry is REPORTED, never targeted: see the module docstring.
    print("dihedral Cl-C-C-H, every beta hydrogen (0 = eclipsing Cl, 180 = anti):")
    for methyl in METHYL_C:
        hydrogens = sorted(h for h in BETA_H if atoms.get_distance(methyl, h) < 1.3)
        angles = {h: atoms.get_dihedral(CHLORINE, CENTRAL_C, methyl, h) for h in hydrogens}
        # ASE returns [0, 360), so "closest to anti" is a distance to 180, not a max of
        # the absolute value -- which picks 300 over 180 and did, in the first draft.
        anti = min(angles, key=lambda h: abs(angles[h] - 180.0))
        print(
            f"  methyl C{methyl}: "
            + "  ".join(f"H{h} {angles[h]:6.1f}" for h in hydrogens)
            + f"   -> closest to anti: H{anti} "
            f"({atoms.get_distance(OXYGEN, anti):.3f} A from O)"
        )
    print()
    print("fragments        2 (the search must be able to reach 1, 2 and 3 -- see README)")
    # The rule is enforced, not merely printed. Drafts of this file violated it twice:
    # the eclipsed rotamer put every beta hydrogen out of proposal reach, and the
    # attempted fix (4.30 A with a wider cutoff) put both channel pairs out of the
    # direction budget's selection instead.
    if not within_c or not within_h:
        raise SystemExit(
            "placement rule 1/3 violated: both channels must be inside the proposal reach"
        )
    from prrs.state import encode

    fragments = encode(atoms, 1.2).components
    print(f"  bond graph fragments                       {fragments}")
    if nearest_any[0] <= BONDING_RANGE_A or fragments != 2:
        raise SystemExit(
            f"placement rule 2 violated: nearest contact {nearest_any[0]:.3f} A, "
            f"{fragments} fragments -- this is a bound adduct, not a separated pair"
        )
    if len(argv) > 1:
        write(str(Path(argv[1])), atoms)
        print(f"written          {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
