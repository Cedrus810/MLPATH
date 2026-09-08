"""Chemical-state identity is not reaction-event identity.

The chemical key quotients away every redundancy of the state, and it must: transferring
one of three equivalent hydrogens has to give one product, not three. That quotient is
also blind to the path, so it cannot decide whether a reaction occurred -- a proton moving
between two equivalent oxygens comes back with the identical key over a bond set that
changed. Both questions are real and they need separate answers, which is what the bond
delta supplies alongside the key.

The event identity then needs its own quotient, for the same reason the state identity did:
raw atom indices would let equivalent atoms manufacture duplicate channels.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.io import read
from pathlib import Path
from prrs.chemistry import automorphisms, chemical_key, classify_transition, reaction_event_key
from prrs.state import encode

ETHANOL = Path(__file__).resolve().parents[1] / "runs" / "e2e_ethanol" / "input.extxyz"


def _malonaldehyde(proton_on):
    """Fixed atom identities; only the bridging proton moves between the two oxygens."""
    o1 = np.array([-1.10, 1.15, 0.0])
    o2 = np.array([1.24, 1.28, 0.0])
    skeleton = [
        o1,
        np.array([-1.19, -0.05, 0.0]),
        np.array([0.00, -0.80, 0.0]),
        np.array([1.19, -0.05, 0.0]),
        o2,
    ]
    host, other = (o1, o2) if proton_on == "O1" else (o2, o1)
    bridging = host + 0.98 * (other - host) / np.linalg.norm(other - host)
    hydrogens = [
        np.array([-2.15, -0.55, 0.0]),
        np.array([0.00, -1.88, 0.0]),
        np.array([2.15, -0.55, 0.0]),
    ]
    return Atoms("OCCCOHHHH", positions=skeleton + [bridging] + hydrogens)


@pytest.mark.parametrize(
    "same_key,broken,formed,expected",
    [
        (True, [], [], "conformational_transition"),
        (True, [[4, 5]], [[0, 5]], "degenerate_reaction"),
        (False, [[4, 5]], [[0, 5]], "reaction"),
        (False, [], [], "key_change_without_bond_change"),
    ],
)
def test_every_cell_of_the_classification_table_is_named(same_key, broken, formed, expected):
    """The fourth cell cannot happen if the key is a function of the bond graph, so it is
    named to be diagnosed rather than silently folded into one of the other three."""
    assert classify_transition("a", "a" if same_key else "b", broken, formed) == expected


def test_a_proton_transfer_changes_bonds_even_when_the_key_may_not():
    """The whole point: a vanishing bond delta is a conformational change, and this is not
    one, whatever the key does."""
    before, after = _malonaldehyde("O2"), _malonaldehyde("O1")
    event = reaction_event_key(before, after)
    assert event["broken"] == [[4, 5]]
    assert event["formed"] == [[0, 5]]
    # Same graph up to relabelling -- the flip O1<->O2, C1<->C3 is an isomorphism -- so the
    # graph hash agrees while the bonds this particular proton holds do not.
    assert chemical_key(before)["graph_hash"] == chemical_key(after)["graph_hash"]
    assert (
        classify_transition("k", "k", event["broken"], event["formed"]) == "degenerate_reaction"
    )


def test_the_event_key_is_the_same_in_both_directions():
    """Forward and reverse are two observations of one channel, not two channels, which is
    the same undirected treatment a saddle's two endpoints already get."""
    before, after = _malonaldehyde("O2"), _malonaldehyde("O1")
    forward = reaction_event_key(before, after)
    reverse = reaction_event_key(after, before)
    assert forward["key"] == reverse["key"]
    assert forward["broken"] == reverse["formed"]
    assert forward["formed"] == reverse["broken"]


def test_equivalent_atoms_do_not_manufacture_duplicate_channels():
    """Breaking any one of three equivalent methyl C-H bonds is one channel, not three.

    Without the quotient the raw indices differ and each transfer opens its own channel,
    which is the reaction-layer version of the duplicate-conformer problem the chemical key
    and the symmetry-aware RMSD already solve one level up.
    """
    ethanol = read(ETHANOL)
    graph = encode(ethanol)
    permutations_, info = automorphisms(ethanol)
    assert info["enumerated"] and len(permutations_) > 1, info

    carbons = [i for i, z in enumerate(ethanol.numbers) if z == 6]
    methyl = max(
        carbons,
        key=lambda c: sum(
            1 for e in graph.edges if c in e and ethanol.numbers[sum(e) - c] == 1
        ),
    )
    equivalent = [
        sum(e) - methyl
        for e in graph.edges
        if methyl in e and ethanol.numbers[sum(e) - methyl] == 1
    ]
    assert len(equivalent) == 3, equivalent

    keys = set()
    for hydrogen in equivalent:
        pulled = ethanol.copy()
        away = pulled.positions[hydrogen] - pulled.positions[methyl]
        pulled.positions[hydrogen] = (
            pulled.positions[methyl] + away / np.linalg.norm(away) * 4.0
        )
        event = reaction_event_key(ethanol, pulled)
        assert event["broken"] == [sorted([methyl, hydrogen])]
        keys.add(event["key"])
    assert len(keys) == 1, "three equivalent hydrogens produced three channels"


def test_a_methyl_rotation_is_not_a_reaction():
    """Ethanol's measured self-connected barrier must stay a conformational transition."""
    from prrs.perturbations import rotatable_torsions, step_torsion
    from prrs.chemistry import canonical_labels

    ethanol = read(ETHANOL)
    graph = encode(ethanol)
    labels = canonical_labels(ethanol.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    assert rotors, "expected a symmetric rotor"
    rotated, _ = step_torsion(ethanol, rotors[0], 2 * np.pi / 3, in_place=False)
    event = reaction_event_key(ethanol, rotated)
    assert event["broken"] == [] and event["formed"] == []
    before, after = chemical_key(ethanol)["key"], chemical_key(rotated)["key"]
    assert (
        classify_transition(before, after, event["broken"], event["formed"])
        == "conformational_transition"
    )


def test_the_same_bond_pattern_in_a_different_substance_is_a_different_channel():
    """Measured collision: malonaldehyde and 3-oxobutanal both break O-H at the same
    indices and form O-H at the same indices, so the bond pattern alone gave them one key.

    A channel identity has to name the substance as well as the rearrangement, which the
    unordered pair of endpoint graph hashes supplies without giving the key a direction.
    """
    before, after = _malonaldehyde("O2"), _malonaldehyde("O1")
    substituted_before, substituted_after = before.copy(), after.copy()
    for structure in (substituted_before, substituted_after):
        numbers = structure.numbers.copy()
        numbers[8] = 9  # one aldehyde H becomes fluorine
        structure.set_atomic_numbers(numbers)

    plain = reaction_event_key(before, after)
    substituted = reaction_event_key(substituted_before, substituted_after)
    assert plain["broken"] == substituted["broken"]
    assert plain["formed"] == substituted["formed"]
    assert plain["key"] != substituted["key"]
    assert plain["endpoint_graphs"] != substituted["endpoint_graphs"]
    # Still direction-free after the tagging.
    assert (
        substituted["key"] == reaction_event_key(substituted_after, substituted_before)["key"]
    )


# ---------------------------------------------------------------------------
# Which bonds can bear E/Z is an energetic question, not a bond length.
# ---------------------------------------------------------------------------

P1 = Path(__file__).resolve().parents[1] / "runs" / "p1_malonaldehyde"


def test_bond_order_semantics_give_the_two_structures_one_key():
    """The defect P1 found, and the fix, pinned against the structures that found it.

    These two relaxed geometries differ by less than 0.0004 A in every bond and by 0.07 meV
    in energy, and the enol C-OH straddles the old 0.93*(r_C + r_O) = 1.3206 A line by
    +0.00021 and -0.00011 A. Deciding E/Z eligibility from the bond graph instead makes the
    question conformation-independent, so the same answer comes out of both.
    """
    from ase.data import covalent_radii

    a = read(P1 / "c0000_m0000.extxyz")
    b = read(P1 / "c0001_m0000.extxyz")
    limit = 0.93 * (covalent_radii[6] + covalent_radii[8])
    lengths = [np.linalg.norm(s.positions[3] - s.positions[4]) for s in (a, b)]
    assert lengths[0] > limit > lengths[1], lengths  # the geometry that split them
    assert abs(lengths[0] - lengths[1]) < 4e-4
    assert chemical_key(a)["key"] == chemical_key(b)["key"]
    assert chemical_key(a)["stereo_unresolved"] is None


def test_the_admissible_set_is_conformation_independent():
    """The measurement that killed the stiffness proposal: the same C=C spans 0.66 to 4.59
    eV/rad^2 across conformers of one substance, because opening the intramolecular
    hydrogen bond destroys the conjugation. A graph-level assignment does not move at all.
    """
    from prrs.chemistry import admissible_bond_orders, locked_edges
    from prrs.state import encode

    sets = []
    for name in ("c0000_m0000.extxyz", "c0001_m0000.extxyz"):
        atoms = read(P1 / name)
        graph = encode(atoms)
        solutions, info = admissible_bond_orders(atoms.numbers, graph.edges, graph.index)
        assert info["count"] == 1, info  # the enol has one Lewis structure
        locked, _ = locked_edges(atoms.numbers, graph.edges, graph.index)
        sets.append(locked)
    assert sets[0] == sets[1] == frozenset({(0, 1), (2, 3)})


def test_a_bond_that_is_double_only_in_some_kekule_structures_carries_no_parity():
    """Benzene has two admissible assignments and no edge is double in both, so no ring
    bond gets an E/Z. Reading one assignment would have made that an artefact of the pick.
    """
    from ase.build import molecule
    from prrs.chemistry import admissible_bond_orders, locked_edges
    from prrs.state import encode

    for name, expected_solutions, expected_locked in (
        ("C6H6", 2, 0),
        ("C2H4", 1, 1),
        ("C2H6", 1, 0),
    ):
        atoms = molecule(name)
        graph = encode(atoms)
        solutions, info = admissible_bond_orders(atoms.numbers, graph.edges, graph.index)
        assert info["count"] == expected_solutions, (name, info)
        locked, _ = locked_edges(atoms.numbers, graph.edges, graph.index)
        assert len(locked) == expected_locked, (name, locked)


def test_an_unsupported_valence_is_stated_rather_than_guessed():
    """Unknown is recorded as unknown, not as "no stereochemical difference".

    This assertion was `locked_bond_parity == []` until 2026-09-04, which is the guess the
    test's own name forbids: [] is exactly what a molecule with genuinely no locked bonds
    produces, so the field could not tell "there are none" from "I could not tell". The
    key already separated them through `stereo_unresolved`, so no identity was ever wrong
    -- but the field was read by the v2 digest, and a regression net cannot review a
    migration through a field that conflates the two states.

    Nitrogen is deliberately outside `VALENCE` (its valence is not single-valued in that
    setting), so an amine is the cheapest way to reach the unresolved branch.
    """
    from ase import Atoms
    from prrs.chemistry import locked_edges
    from prrs.state import encode

    amine = Atoms("NH3", positions=[[0, 0, 0], [1.0, 0, 0], [-0.3, 0.95, 0], [-0.3, -0.5, 0.8]])
    graph = encode(amine)
    locked, info = locked_edges(amine.numbers, graph.edges, graph.index)
    assert locked is None
    assert info["reason"] == "unsupported_elements" and info["elements"] == ["N"]
    key = chemical_key(amine)
    assert key["stereo_unresolved"] == "unsupported_elements"
    assert key["locked_bond_parity"] is None, "unknown must not be reported as empty"

    # And the contrast that gives that assertion its meaning: a molecule the annotator
    # CAN read, with no locked bonds, still reports [].
    methane = Atoms(
        "CH4",
        positions=[
            [0, 0, 0],
            [0.63, 0.63, 0.63],
            [0.63, -0.63, -0.63],
            [-0.63, 0.63, -0.63],
            [-0.63, -0.63, 0.63],
        ],
    )
    resolved = chemical_key(methane)
    assert resolved["stereo_unresolved"] is None
    assert resolved["locked_bond_parity"] == []


def test_an_open_shell_endpoint_is_refused_rather_than_admitted():
    """A radical pair is outside a potential trained on neutral closed-shell molecules.

    Measured: seeded from a malonaldehyde probe frame, an unbounded climb reached an O-H
    homolysis 2.34 eV up, and descending it produced C3H3O2 plus a hydrogen atom. The
    potential ignores total charge outright and was never trained on open shells, so
    admitting that as a product would put a structure it cannot describe into the network.
    The bond-order solver already says so -- no assignment satisfies the valences -- and the
    registry now refuses on that verdict rather than recording it and moving on.
    """
    from ase.calculators.singlepoint import SinglePointCalculator
    from prrs.config import SearchConfig
    from prrs.network import Registry

    def scored(atoms):
        result = atoms.copy()
        result.calc = SinglePointCalculator(
            result, energy=0.0, forces=np.zeros((len(atoms), 3))
        )
        return result

    radical = _malonaldehyde("O2")
    positions = radical.positions.copy()
    positions[5] += np.array([0.0, 6.0, 0.0])  # pull the proton clean off
    radical.set_positions(positions)
    key = chemical_key(radical)
    assert key["stereo_unresolved"] == "no_assignment_satisfies_valences"
    assert len(key["fragments"]) == 2

    registry = Registry(SearchConfig(closed_shell_only=True))
    admission = registry.admit(scored(radical), None, "t0")
    assert admission.outcome == "rejected"
    assert admission.reason == "out_of_domain"
    assert registry.nodes == []

    # The switch exists because a charge- and spin-aware potential could handle it.
    permissive = Registry(SearchConfig(closed_shell_only=False))
    assert permissive.admit(scored(radical), None, "t0").outcome == "new_chemical_node"

    # An element the solver does not cover is a different verdict and must not be refused.
    from ase import Atoms

    amine = Atoms("NH3", positions=[[0, 0, 0], [1.0, 0, 0], [-0.3, 0.95, 0], [-0.3, -0.5, 0.8]])
    assert chemical_key(amine)["stereo_unresolved"] == "unsupported_elements"
    assert (
        Registry(SearchConfig(closed_shell_only=True)).admit(scored(amine), None, "t0").outcome
        == "new_chemical_node"
    )


def test_the_fourth_cell_is_split_by_the_key_components():
    """Key changed, no bond changed: real E/Z, or a defect in the key.

    Named `key_change_without_bond_change` on the assumption that the cell is unreachable
    if the key is a function of the bond graph. The key is deliberately not that -- it
    carries configurational stereochemistry too -- so E/Z isomerisation lands there
    legitimately. Measured twice independently on 3-oxobutanal: same locked bond, parity
    -1 -> +1, graph hash unchanged, 365-404 meV above the reactant.

    The split needs the components. Given only the digests the conservative name is kept,
    because "I cannot tell" must not be reported as "it is fine".
    """
    base = {
        "graph_hash": "g",
        "fragments": ["C4H6O2"],
        "tetrahedral_parity": [],
        "locked_bond_parity": [[["x", "y"], -1]],
        "stereo_unresolved": None,
    }
    flipped = dict(base, locked_bond_parity=[[["x", "y"], 1]])
    other_graph = dict(base, graph_hash="h")
    relabelled = dict(base, locked_bond_parity=[[["x", "z"], -1]])

    assert (
        classify_transition("a", "b", [], [], base, flipped) == "configurational_isomerisation"
    )
    # No components: cannot tell, so keep the name that says "diagnose me".
    assert classify_transition("a", "b", [], []) == "key_change_without_bond_change"
    # A key change with something else also different is not stereochemistry.
    assert (
        classify_transition("a", "b", [], [], base, other_graph)
        == "key_change_without_bond_change"
    )
    # Parity defined on a different bond is not a flip of the same one.
    assert (
        classify_transition("a", "b", [], [], base, relabelled)
        == "key_change_without_bond_change"
    )
    # The other three cells are untouched.
    assert classify_transition("a", "b", [[0, 5]], [[4, 5]], base, flipped) == "reaction"
    assert classify_transition("a", "a", [[0, 5]], [[4, 5]]) == "degenerate_reaction"
    assert classify_transition("a", "a", [], []) == "conformational_transition"


def test_the_model_and_the_acceptance_checker_share_one_stereo_test():
    """Two copies of this test is how the checker and the model come to disagree about
    what a class means. The checker imports the model's implementation."""
    import importlib.util
    from prrs.chemistry import stereo_only_difference

    spec = importlib.util.spec_from_file_location(
        "p2_check_probe", Path(__file__).resolve().parents[1] / "docs/experiments/p2_check.py"
    )
    source = spec.origin
    text = Path(source).read_text(encoding="utf-8")
    assert "from prrs.chemistry import stereo_only_difference" in text
    assert "def stereo_only_difference" not in text, "the checker has its own copy again"
    ok, why = stereo_only_difference(
        {
            "graph_hash": "g",
            "fragments": ["x"],
            "locked_bond_parity": [[["a", "b"], -1]],
            "tetrahedral_parity": [],
        },
        {
            "graph_hash": "g",
            "fragments": ["x"],
            "locked_bond_parity": [[["a", "b"], 1]],
            "tetrahedral_parity": [],
        },
    )
    assert ok and "parity flipped" in why


def test_both_callers_of_the_stereo_test_get_the_same_answer():
    """The two callers do not pass the same shape, and that broke the split once.

    `publish()` strips "key" out of `key_components` before writing network.json, so the
    acceptance checker reads a dict without it while the model passes `node.components`
    with it still present. Measured on P2 T=300 seed 71: the checker reported "parity
    flipped, everything else identical" and the model reported "key also differs" about the
    same two nodes, so a real E/Z isomerisation was recorded under the defect name while the
    checker insisted it was fine.

    Sharing one implementation is not the same as sharing one input convention.
    """
    from prrs.chemistry import stereo_only_difference

    stored = {
        "graph_hash": "g",
        "fragments": ["C4H6O2"],
        "tetrahedral_parity": [],
        "locked_bond_parity": [[["x", "y"], -1]],
        "stereo_unresolved": None,
        "charge_sensitive": False,
        "note": "...",
        "automorphisms": {"count": 6},
    }
    flipped = dict(stored, locked_bond_parity=[[["x", "y"], 1]])
    from_checker = stereo_only_difference(stored, flipped)
    # The runtime shape: the same components plus the digest, which necessarily differs.
    from_model = stereo_only_difference(dict(stored, key="aaa"), dict(flipped, key="bbb"))
    assert from_checker[0] is True, from_checker
    assert from_model[0] is True, from_model
    assert from_checker[0] == from_model[0]
    # And a key change that is NOT a parity flip is still refused in both shapes.
    other = dict(stored, graph_hash="h")
    assert stereo_only_difference(stored, other)[0] is False
    assert stereo_only_difference(dict(stored, key="aaa"), dict(other, key="bbb"))[0] is False


def test_a_degenerate_reaction_is_published_as_a_self_loop_channel(tmp_path):
    """A degenerate transfer must reach the network as a channel, not vanish into it.

    The gap this closes was found by mutation on 2026-09-04. Suppressing the self-loop
    branch of `record_channel` -- so a degenerate reaction opens no channel at all --
    left the whole suite green (225 passed) and left the digest green too. Every existing
    guard was one level too high:

      * `classify_transition` returning "degenerate_reaction" tests the CLASSIFIER;
      * "P1's two structures share one key" tests STATE identity;
      * neither says the event was ever PUBLISHED.

    A degenerate reaction is precisely the case where state identity cannot carry the
    evidence: both ends are the same substance, so the only record that the reaction
    happened is the channel. If it is not published there is nothing left.

    The system is the SharedProton double well from test_boundary_continuation: two
    identical heavy atoms with a proton between them, so the two minima are mirror
    images with the same chemical key and different bond sets.
    """
    import sys
    from pathlib import Path as _Path

    sys.path.insert(0, str(_Path(__file__).parent))
    from test_boundary_continuation import SharedProton, shared_proton_atoms

    from prrs.config import SearchConfig
    from prrs.search import ReactionSearch

    config = SearchConfig(
        seed=11,
        max_trials=24,
        max_depth=1,
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.3, 0.6, 0.9),
        quench_fmax_eV_A=0.01,
        closed_shell_only=False,
        minimum_check="probe",
        temperature_K=0.0,
        # Off by construction, not for convenience: this is C-H-C with the sum r01 + r12
        # held by the potential, so driving either bond necessarily moves the other. The
        # collateral gate exists to catch a probe that tore something it was not aiming
        # at; here the coupling IS the coordinate.
        collateral_gate=False,
    )
    network = ReactionSearch(SharedProton, config).run(
        shared_proton_atoms(-0.4), tmp_path / "degenerate"
    )

    channels = network["reaction_channels"]
    assert channels, (
        "a degenerate transfer produced no reaction_channel at all; both ends are the "
        "same substance, so the channel is the only place the event can be recorded"
    )
    loops = [c for c in channels if c["self_loop"]]
    assert loops, f"no self-loop channel among {[c['class'] for c in channels]}"
    loop = loops[0]
    assert loop["class"] == "degenerate_reaction", loop["class"]
    assert loop["ends"][0] == loop["ends"][1], loop["ends"]
    assert loop["broken"] and loop["formed"], (
        "a self-loop with no bond delta is a conformational change, not a reaction"
    )
    # The directed-edge layer must stay empty: same chemical state at both ends.
    assert not network["reactions"], (
        "a degenerate reaction opens no directed edge -- see HANDOFF.md section 3"
    )


def test_a_budget_overflow_is_refused_rather_than_named():
    """An identity that depends on the enumeration budget is not an identity.

    Two failure modes were being reported through one field. "This molecule's
    stereochemistry cannot be resolved" -- an element outside the valence table, a graph
    no neutral closed-shell assignment satisfies -- is a property of the molecule, and it
    belongs in the key so a resolved substance is never compared against an unresolved one
    on the quiet. "The enumerator ran out of budget" is a property of the RUN, and putting
    it in the key makes chemical identity a function of how much compute was allowed.

    Measured on the P3 A structure before the fix:

        bond_order_limit = 1000  ->  key 8710064d12ce6a5a
        bond_order_limit =  100  ->  key 38584b6ee5c4eca2

    Same molecule, same geometry, same code. The event key had the same defect more
    quietly: when `automorphisms` gives up it returns the identity permutation alone, so
    the bond delta is canonicalised against an incomplete group and the key moves with no
    field in the hash recording that it did.

    No recorded result is affected -- every frame in the digest corpus resolves inside the
    default budget and all four recorded event pairs enumerate completely -- so this pins
    behaviour in a case that has not occurred rather than changing one that has.
    """
    from prrs.chemistry import IdentityUnavailable

    a, b = Path("runs/p3_A_source.extxyz"), Path("runs/p3_B_source.extxyz")
    if not (a.exists() and b.exists()):
        pytest.skip("P3 corpus absent")
    source, endpoint = read(str(a)), read(str(b))

    # The generous budget is unchanged: this must keep working exactly as before.
    assert chemical_key(source, 1.2)["key"].startswith("8710064d12ce6a5a")
    assert reaction_event_key(source, endpoint, 1.2)["key"] == "6e1ee45ce5a39fd8"

    # A budget that cannot finish must refuse, not name.
    with pytest.raises(IdentityUnavailable, match="enumeration_limit"):
        chemical_key(source, 1.2, bond_order_limit=100)
    with pytest.raises(IdentityUnavailable, match="automorphism"):
        reaction_event_key(source, endpoint, 1.2, None, 1)

    # A CHEMICAL resolution failure is still recorded, not raised: that one is a fact
    # about the molecule and belongs in the key.
    chloride = Atoms("Cl", positions=[[0, 0, 0]])
    assert (
        chemical_key(chloride, 1.2)["stereo_unresolved"] == "no_assignment_satisfies_valences"
    )


def test_a_confirmed_minimum_refused_by_a_label_is_still_retrievable():
    """A verified stationary point may be kept out of the network, never lost.

    `Registry.admit` can refuse an endpoint on a chemical label -- no neutral closed-shell
    valence assignment, under a model that declares it needs one. The refusal happens
    AFTER `confirm_minimum` has passed, so what is excluded is a confirmed stationary
    point on the surface the search ran on.

    This pins the recorded cases rather than a synthetic one. A synthetic C-H-C or O-H-O
    double well is itself a radical: its SOURCE fails the same check and the search never
    starts, so the real occurrences are the honest fixture.

    Measured across every recorded run: three refusals, in p2_T300 seeds 31, 67 and 73.
    In all three the endpoint structure is on disk -- but before `withheld_candidates` it
    was reachable only by reading attempt records one at a time, so a QM verification pass
    had no list to work from. Structures and stationary points are the primary results,
    the chemical label is auxiliary, and "the annotator could not name it" is not grounds
    for the minimum to stop existing.

    NOT covered here: that `Search.run` populates the new list. No recorded run predates
    it and the synthetic route is blocked above; it is exercised the first time a run
    refuses an endpoint.
    """
    import json

    found = []
    for seed in (31, 67, 73):
        run = Path(f"runs/p2_T300_seed{seed}")
        if not (run / "attempts.jsonl").exists():
            continue
        for line in (run / "attempts.jsonl").read_text().splitlines():
            record = json.loads(line)
            if record.get("network_exclusion") != "out_of_domain":
                continue
            found.append(seed)
            # It reached admission, so confirm_minimum had already accepted it.
            assert record["status"] == "completed", record["status"]
            # The structure survives the refusal.
            assert (run / record["endpoint"]).exists(), record["endpoint"]
            assert (run / record["trajectory"]).exists(), record["trajectory"]
    if not found:
        pytest.skip("p2_T300 corpus absent")
    assert sorted(found) == [31, 67, 73], f"expected the three recorded refusals, got {found}"


def test_every_rejection_path_keeps_the_geometry_and_says_where_it_came_from():
    """Whatever keeps a geometry out of the network, the geometry survives with its source.

    Three paths reach `Registry.admit`, and they used to end differently. `execute` listed
    the candidate; the two-sided saddle descent set a reason string on its own entry and
    dropped the endpoint; the boundary continuation did the same. So a minimum that
    `confirm_minimum` had already accepted could disappear because a chemical label refused
    it, and nothing downstream could learn it had existed -- which is exactly what a QM
    verification pass needs to look at.

    All four now go through one `withhold()`: geometry written under `withheld/`,
    provenance in `source`, reason and detail kept. A descent endpoint and a trial endpoint
    are not the same kind of evidence, so a reader must not have to guess which is which.

    This pins the contract at the source level, since reaching every branch needs a system
    whose source is admitted while some endpoint is not -- see the sibling test for why a
    synthetic one is not available.
    """
    import inspect
    import re

    from prrs import search

    body = inspect.getsource(search.ReactionSearch.run)
    assert body.count('network["withheld_candidates"].append') == 1, (
        "more than one place appends candidates; they will drift apart"
    )
    sources = set(re.findall(r'withhold\(\s*[\w.]+,\s*"(\w+)"', body))
    assert sources == {"trial_endpoint", "saddle_descent", "boundary_continuation"}, sources

    # The two paths that used to drop the endpoint now hand it to withhold().
    for marker in (
        'entry["reason"] = admission.reason',
        'entry["admission_reason"] = admission.reason',
    ):
        after = body[body.index(marker) : body.index(marker) + 400]
        assert "withhold(" in after, f"rejection at {marker!r} still drops the geometry"
