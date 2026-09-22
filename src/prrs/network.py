"""Two-layer reaction network: chemical macrostates holding conformer microstates.

Flattening every distinct minimum into one node does not survive contact with real
molecules. Ethanol has three conformer minima within 47 meV; anything larger has
combinatorially many, and a flat network spends its whole budget enumerating them and
draws reaction edges between conformers of one substance. So:

    chemical node    one substance, identified by chemical_key
      microstates    distinct conformers of that substance, a bounded reservoir

Admission is decided by chemical_key, never by the response label. An equivalent
hydrogen transferring is classified "reactive" because the graph changed, yet it is the
same chemistry; only a changed key opens a chemical node.

Budgets are separate. Discovering a ninth conformer of a known substance must not
consume a chemical node slot, because conformational richness is not chemical
discovery. But conformers stay searchable: reaction accessibility is conditional on
conformer -- anti may be unreactive where gauche is not -- so a bounded set of them
remains available as perturbation origins.
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import ast
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from .chemistry import BUDGET_REASONS, automorphisms, chemical_key, free_aligned_symmetric_rmsd
from .perturbations import rotatable_torsions
from . import internal
from .state import encode
from .canonical import canonical_form

# Prose that ships inside the records this module writes. Held as named constants
# rather than inline in the dict literals so the shape of each record is visible at
# a glance; the text is part of the output contract, so changing one is a schema
# change and not a comment edit.
_ADMIT_MEANING = (
    "the Lewis annotator found no neutral closed-shell valence assignment, and this model "
    "declares it requires one (closed_shell_only). Two statements, not one: the annotator "
    "supplied the evidence, the model's declaration made it a refusal. A charge-aware "
    "potential declares closed_shell_only = False and the same structure is admitted with "
    "its stereochemistry recorded as unresolved"
)


def torsion_profile(atoms, bond_scale=1.2, active=None):
    """Torsions that carry real conformational freedom.

    Symmetry rotors are excluded on purpose. A methyl turn changes its own torsion by
    120 degrees while producing a structure that is a symmetry copy of the original, so
    counting it would make copies look diverse and let them win reservoir slots.
    """
    graph = encode(atoms, bond_scale, active=active)
    form = canonical_form(atoms.numbers, graph.edges, graph.index)
    labels = form.colours
    orbits = form.orbits if form.info.get("enumerated") else labels
    genuine, _rotors = rotatable_torsions(graph, labels, orbits)
    return {
        tuple(t.indices): internal.coordinate(atoms.positions, "dihedral", t.indices)
        for t in genuine
    }


def torsion_distance(a, b):
    """RMS circular difference. Angles are periodic; a plain difference is not a metric."""
    shared = sorted(set(a) & set(b))
    if not shared:
        return 0.0
    deltas = [abs((a[k] - b[k] + np.pi) % (2 * np.pi) - np.pi) for k in shared]
    return float(np.sqrt(np.mean(np.square(deltas))))


def free_bonds_of(structure):
    """Bonds whose torsion the quench found flat, read off the structure that was quenched.

    The mode-resolved polish already measures every rotatable torsion's stiffness and
    labels the ones below the floor "free"; that verdict belongs to the geometry it was
    measured on, so it travels with it rather than being recomputed. A structure that
    arrived without a quench report contributes nothing, and comparison falls back to the
    plain symmetry-aware RMSD -- which is the behaviour everywhere no coordinate is free.
    """
    rounds = (structure.info.get("quench") or {}).get("rounds") or []
    if not rounds:
        return ()
    bonds = set()
    for mode in rounds[-1].get("modes", []):
        # The polish reports two kinds of coordinate now, and only one of them has a bond.
        # `indices[1:3]` off a fragment rotation is the fragment's second and third ATOM,
        # which need not be bonded to each other and certainly is not a free rotor. Entries
        # written before the field existed are torsions, which is what the default says.
        if mode.get("coordinate", "torsion") != "torsion":
            continue
        if mode.get("tier") == "free":
            indices = mode["indices"]
            bonds.add(tuple(sorted(indices[1:3])))
    return tuple(sorted(bonds))


def measured_bonds_of(structure):
    """Every torsion the quench measured, with the verdict it reached.

    free_bonds_of answers "which are flat"; this answers "which were looked at", and the
    difference matters when two structures are compared. A bond absent here was never
    measured on this geometry, which is not the same as having been measured and found
    stiff -- the same measured/not-measured split the identity layer had to make
    (docs/DECISIONS_PENDING.md, "unknown is not empty").
    """
    rounds = (structure.info.get("quench") or {}).get("rounds") or []
    if not rounds:
        return {}
    return {
        tuple(sorted(mode["indices"][1:3])): mode.get("tier")
        for mode in rounds[-1].get("modes", [])
        if mode.get("coordinate", "torsion") == "torsion"  # same reason as free_bonds_of
    }


@dataclass
class Microstate:
    id: str
    structure: object
    energy_eV: float
    discovered_by: str
    torsions: dict
    free_bonds: tuple = ()
    trials_spent: int = 0
    outcomes: Counter = field(default_factory=Counter)
    score: float = None

    def summary(self, relative_path):
        return {
            "id": self.id,
            "structure": relative_path,
            "energy_eV": self.energy_eV,
            "discovered_by": self.discovered_by,
            "trials_spent": self.trials_spent,
            "torsions_rad": {str(list(k)): v for k, v in self.torsions.items()},
            "free_bonds": [list(bond) for bond in self.free_bonds],
            "outcomes": dict(self.outcomes),
            "reservoir_score": self.score,
        }


@dataclass
class ChemicalNode:
    id: str
    key: str
    components: dict
    depth: int
    discovered_by: str
    permutations: list
    microstates: list = field(default_factory=list)
    evicted: list = field(default_factory=list)
    trials_spent: int = 0
    # Saddle-connected transitions between conformers of this substance. They are not
    # reactions and do not belong in the chemical network, but they are real information
    # about how the reservoir's members interconvert.
    conformer_transitions: list = field(default_factory=list)


@dataclass
class Admission:
    outcome: str
    node: ChemicalNode = None
    microstate: Microstate = None
    reason: str = None
    detail: dict = field(default_factory=dict)


class Registry:
    """Owns the two layers, the reservoir policy and both budgets."""

    def __init__(self, config):
        self.config = config
        self.nodes = []

    # ---- identity -------------------------------------------------------------
    def key_of(self, structure):
        return chemical_key(
            structure,
            self.config.bond_scale,
            self.config.active_atoms,
            charge_sensitive=self.config.charge_sensitive,
            charge=self.config.total_charge,
            multiplicity=self.config.multiplicity,
            parity_tolerance=self.config.parity_tolerance,
        )

    def find_node(self, key):
        return next((node for node in self.nodes if node.key == key), None)

    def identify(self, structure):
        """Which (node, microstate) this geometry already IS. Reads; never inserts.

        `admit` is the wrong tool for asking "where did this land": it mints nodes and
        microstates, so using it to judge a probe would write the probe into the very
        network the probe is supposed to observe. The committor needs exactly this
        question and nothing else (PLAN item 2.2, decision 3, 2026-09-20: the shots are
        an independent probe, model report only, never a criterion) -- and asking it
        through `admit` was how the first production run produced A 0 / B 0 / other 16
        on every bracket: nothing set `target_microstate`, so both named outcomes were
        unreachable by construction.

        Returns (None, None) for a geometry no existing node matches, which is a third
        basin or an unrecognised one -- never silently the nearest of the two.
        """
        # key_of returns the decomposed identity; the node carries only its hash, so the
        # ["key"] is the comparable half -- same extraction `admit` makes. Handing the
        # whole dict to find_node compares a dict against a string and silently matches
        # nothing, which is the second way this judging could report "no basin" for a
        # geometry sitting squarely in one.
        node = self.find_node(self.key_of(structure)["key"])
        if node is None:
            return None, None
        return node, self._match_microstate(node, structure)

    # ---- reservoir policy -----------------------------------------------------
    def _energy_term(self, node, energy):
        floor = (
            min(state.energy_eV for state in node.microstates) if node.microstates else energy
        )
        from ase import units

        scale = units.kB * max(self.config.conformer_temperature_K, 1e-6)
        return float(np.exp(-max(energy - floor, 0.0) / scale))

    def _diversity_term(self, node, torsions, exclude=None):
        others = [s for s in node.microstates if s is not exclude]
        if not others:
            return 1.0
        nearest = min(torsion_distance(torsions, s.torsions) for s in others)
        return float(min(nearest / np.pi, 1.0))

    def _response_term(self, node, state):
        """Optimistic until searched, then rewards a conformer that responds differently.

        The response term cannot be known before budget is spent on a microstate, so an
        unsearched one is scored as maximally informative. That is deliberate: it makes
        the reservoir an online selection under uncertainty rather than a clustering of
        what is already known.
        """
        if state is None or state.trials_spent == 0:
            return 1.0
        siblings = Counter()
        for other in node.microstates:
            if other is not state:
                siblings.update(other.outcomes)
        total = sum(state.outcomes.values())
        if not total:
            return 1.0
        unique = sum(
            count for outcome, count in state.outcomes.items() if outcome not in siblings
        )
        return float(unique / total)

    def score(self, node, energy, torsions, state=None):
        cfg = self.config
        return (
            cfg.reservoir_weight_energy * self._energy_term(node, energy)
            + cfg.reservoir_weight_diversity * self._diversity_term(node, torsions, state)
            + cfg.reservoir_weight_response * self._response_term(node, state)
        )

    def rescore(self, node):
        for state in node.microstates:
            state.score = self.score(node, state.energy_eV, state.torsions, state)

    # ---- admission ------------------------------------------------------------
    def admit(self, structure, parent, trial_id):
        """Route a confirmed minimum to a microstate or a chemical node.

        Notes
        -----
        [1] Three separate things used to be answered by one field, and separating them is
            what this branch now does.

            **Identity content** -- the graph, the fragments, the parities. What the
            substance is.
            **Identity resolution status** -- whether the stereochemistry could be worked
            out at all. A statement about the annotator, not about the substance.
            **Model domain evidence** -- whether the potential can describe this species.
            A statement about the model.

            The branch below is the third one, and it is gated on `closed_shell_only`
            because that is a MODEL declaration. Its original justification was explicit
            about this: MACE-OFF is trained on neutral closed-shell molecules and does not
            respond to total charge at all, measured bit for bit, so admitting a charged or
            open-shell species would put a structure that potential cannot describe into
            the network as a product.

            That justification is a property of THAT model. Charge-aware potentials do not
            carry it -- MACE-POLAR-1-M's energy spans 16.67 eV across total charge 0/-1/+1
            -- and for them `closed_shell_only = False` is the correct declaration rather
            than a workaround. The Lewis solver is then used as what it is: an annotator
            whose failure is recorded, never a veto over a minimum the potential can
            evaluate.

            What must NOT happen is the inverse reading: an unresolved annotation is not
            evidence that a species is out of domain, and a resolved one is not evidence
            that it is in. Domain belongs to `benchmarks/*/model.json` and the preflight
            gates -- element coverage read from the checkpoint, measured charge response,
            uncertainty when a committee is available.

            Unsupported elements were never refused here and still are not: that only means
            the stereochemistry is unresolved.
        """
        cfg = self.config
        components = self.key_of(structure)
        # [1] the model's declared domain, not the annotator's verdict
        if cfg.closed_shell_only and components.get("stereo_unresolved") in (
            "no_assignment_satisfies_valences",
            "degree_exceeds_valence",
        ):
            return Admission(
                "rejected",
                reason="out_of_domain",
                detail={
                    "chemical_key": components["key"],
                    "fragments": components["fragments"],
                    "verdict": components["stereo_unresolved"],
                    # Which layer refused, so a reader never has to infer it. The Lewis
                    # solver supplied the evidence; the model's declaration is what made
                    # it a refusal.
                    "refused_by": "model_domain_declaration",
                    "declaration": "closed_shell_only",
                    "meaning": _ADMIT_MEANING,
                },
            )
        node = self.find_node(components["key"])
        energy = float(structure.get_potential_energy())
        torsions = torsion_profile(structure, cfg.bond_scale, cfg.active_atoms)

        if node is None:
            if len(self.nodes) >= cfg.chemical_max_nodes:
                return Admission(
                    "rejected",
                    reason="chemical_budget",
                    detail={"chemical_key": components["key"]},
                )
            permutations_, info = automorphisms(
                structure, cfg.bond_scale, cfg.active_atoms, cfg.automorphism_limit
            )
            if info.get("reason") in BUDGET_REASONS:
                # P1-4 (decision 7, 2026-09-20): on overflow automorphisms() hands back
                # [identity] with the reason in info, and the event path already refuses
                # to canonicalise against that ("An incomplete group is never used").
                # This path used to build the node anyway, and conformer matching then
                # silently degraded to naive RMSD -- a methyl rotation would look like a
                # new conformer forever. Refuse the node instead: the run keeps going,
                # the geometry is withheld with this reason, and nothing downstream
                # reads a group that was never certified.
                return Admission(
                    "rejected",
                    reason="automorphism_budget",
                    detail=dict(info),
                )
            node = ChemicalNode(
                id=f"c{len(self.nodes):04d}",
                key=components["key"],
                components=components,
                depth=0 if parent is None else parent.depth + 1,
                discovered_by=trial_id,
                permutations=permutations_,
            )
            node.components["automorphisms"] = info
            self.nodes.append(node)
            state = self._insert(node, structure, energy, torsions, trial_id)
            return Admission("new_chemical_node", node, state)

        existing = self._match_microstate(node, structure)
        if existing is not None:
            return Admission("existing_microstate", node, existing)

        if len(node.microstates) < cfg.conformer_max_per_node:
            state = self._insert(node, structure, energy, torsions, trial_id)
            outcome = (
                "new_chemical_node" if node.discovered_by == trial_id else "new_microstate"
            )
            return Admission(outcome, node, state)

        # Reservoir full: keep the better set rather than the earlier one, but never
        # evict a microstate that still has unspent search budget.
        candidate = self.score(node, energy, torsions, None)
        self.rescore(node)
        evictable = [
            s for s in node.microstates if s.trials_spent >= cfg.conformer_trials_per_microstate
        ]
        if not evictable:
            return Admission(
                "rejected",
                node,
                reason="conformer_reservoir_full",
                detail={"capacity": cfg.conformer_max_per_node},
            )
        worst = min(evictable, key=lambda s: s.score)
        if candidate <= worst.score:
            return Admission(
                "rejected",
                node,
                reason="conformer_reservoir_outscored",
                detail={"candidate_score": candidate, "weakest_resident_score": worst.score},
            )
        node.microstates.remove(worst)
        node.evicted.append({"id": worst.id, "score": worst.score, "replaced_by": trial_id})
        state = self._insert(node, structure, energy, torsions, trial_id)
        return Admission(
            "new_microstate",
            node,
            state,
            detail={"evicted": worst.id, "candidate_score": candidate},
        )

    def _match_microstate(self, node, structure):
        cfg = self.config
        candidate_free = free_bonds_of(structure)
        candidate_tiers = measured_bonds_of(structure)
        for state in node.microstates:
            if abs(state.energy_eV - structure.get_potential_energy()) > cfg.basin_energy_eV:
                continue
            # The union, because either structure may be the one that revealed the
            # coordinate is flat -- but never against a measurement. Flatness belongs to
            # the geometry it was measured on, so when the other structure measured the
            # same torsion and did NOT find it free, the two geometries disagree and the
            # coordinate is not quotiented. Quotienting it anyway merges two real
            # conformers, which free_aligned_symmetric_rmsd's own docstring names as the
            # thing nothing may infer: measured on ethanol, anti and gauche sit 0.4383 A
            # apart and collapse to 0.0614 A -- one side's report deciding for both.
            # A bond nobody measured raises no objection, which keeps the original
            # purpose: a genuinely free rotor must not mint conformers just because the
            # structure it was measured on is not the one being compared.
            state_tiers = measured_bonds_of(state.structure)
            free = tuple(
                sorted(
                    bond
                    for bond in set(state.free_bonds) | set(candidate_free)
                    if candidate_tiers.get(bond, "free") == "free"
                    and state_tiers.get(bond, "free") == "free"
                )
            )
            distance = free_aligned_symmetric_rmsd(
                state.structure,
                structure,
                node.permutations,
                cfg.active_atoms,
                free_bonds=free,
                bond_scale=cfg.bond_scale,
                samples=cfg.free_alignment_samples,
            )
            if distance <= cfg.basin_rmsd_A:
                return state
        return None

    def _insert(self, node, structure, energy, torsions, trial_id):
        state = Microstate(
            id=f"{node.id}/m{len(node.microstates) + len(node.evicted):04d}",
            structure=structure,
            energy_eV=energy,
            discovered_by=trial_id,
            torsions=torsions,
            free_bonds=free_bonds_of(structure),
        )
        node.microstates.append(state)
        self.rescore(node)
        return state

    # ---- budgets --------------------------------------------------------------
    def may_search(self, node, state):
        cfg = self.config
        if node.trials_spent >= cfg.chemical_trials_per_node:
            return False, "chemical_trials_per_node"
        if state.trials_spent >= cfg.conformer_trials_per_microstate:
            return False, "conformer_trials_per_microstate"
        return True, None

    def note_trial(self, node, state, outcome):
        node.trials_spent += 1
        state.trials_spent += 1
        state.outcomes[outcome] += 1
        self.rescore(node)

    @classmethod
    def from_network(cls, network, output, config):
        """Rebuild a live Registry from a published network.json (PLAN item 6.2).

        Everything a microstate needs is in `Microstate.summary()` plus the structure
        file it points at; everything a node needs is in the published entry except
        one thing -- `permutations`, the automorphism arrays, which are not serialised.
        They are re-enumerated from the node's first structure instead, which is
        legitimate exactly because the enumeration is deterministic and index-ordered:
        same structure, same config, same group, same order. If the enumeration cannot
        be certified inside the budget, reconstruction REFUSES rather than handing the
        reservoir an identity-only group -- the same fail-closed rule `admit` follows
        (P1-4), because a resumed run would otherwise compare conformers with a group
        that was never certified.

        `output` is the run directory the network.json lives in; structure paths in
        the summary are relative to it.
        """
        output = Path(output)
        registry = cls(config)
        for published in network["chemical_nodes"]:
            first_path = output / published["microstates"][0]["structure"]
            if not first_path.exists():
                raise ValueError(
                    f"cannot rebuild node {published['id']}: {first_path} is missing"
                )
            first = read(str(first_path))
            permutations_, info = automorphisms(
                first, config.bond_scale, config.active_atoms, config.automorphism_limit
            )
            if info.get("reason") in BUDGET_REASONS:
                raise ValueError(
                    f"cannot rebuild node {published['id']}: the automorphism group "
                    f"does not certify inside automorphism_limit ({info['reason']})"
                )
            components = dict(published["key_components"])
            components["key"] = published["chemical_key"]
            node = ChemicalNode(
                id=published["id"],
                key=published["chemical_key"],
                components=components,
                depth=published["depth"],
                discovered_by=published["discovered_by"],
                permutations=permutations_,
                evicted=list(published["reservoir"]["evicted"]),
                trials_spent=published["trials_spent"],
                conformer_transitions=list(published["conformer_transitions"]),
            )
            for state in published["microstates"]:
                path = output / state["structure"]
                if not path.exists():
                    raise ValueError(
                        f"cannot rebuild microstate {state['id']}: {path} is missing"
                    )
                structure = read(str(path))
                sidecar = path.with_suffix(".npz")
                if sidecar.exists():
                    # bit-exact positions and forces for resume: extxyz text would
                    # cost the last ulps, and every recomputation downstream of the
                    # rebuild (spectra, references, Hessians) would drift with them
                    payload = np.load(sidecar)
                    structure.set_positions(payload["positions"])
                    forces = payload["forces"]
                    structure.calc = SinglePointCalculator(
                        structure,
                        energy=float(payload["energy"]),
                        forces=forces if len(forces) else None,
                    )
                else:
                    # runs older than resume have no sidecar: the extxyz carries what
                    # it carries, and the summary's energy stands in when it is bare
                    try:
                        structure.get_potential_energy()
                    except Exception:
                        structure.calc = SinglePointCalculator(
                            structure, energy=state["energy_eV"]
                        )
                node.microstates.append(
                    Microstate(
                        id=state["id"],
                        structure=structure,
                        energy_eV=state["energy_eV"],
                        discovered_by=state["discovered_by"],
                        trials_spent=state["trials_spent"],
                        torsions={
                            tuple(ast.literal_eval(key)): value
                            for key, value in state["torsions_rad"].items()
                        },
                        free_bonds=tuple(tuple(bond) for bond in state["free_bonds"]),
                        outcomes=Counter(state["outcomes"]),
                        score=state["reservoir_score"],
                    )
                )
            registry.nodes.append(node)
        return registry
