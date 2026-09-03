"""Budgeted search over a two-layer network, with conformers held inside chemistry.

The controller spends two budgets. A chemical node slot is spent only when a confirmed
minimum has a chemical key no existing node has; finding another conformer of a known
substance spends that node's conformer budget instead. Every confirmed minimum stays
searchable within its budget, because whether A reaches B depends on which conformer of
A was perturbed.

A rejected endpoint is not discarded. Exactly one negative Hessian eigenvalue is a
first-order saddle -- a transition state candidate whose unstable mode is the reaction
coordinate -- so it goes to its own pool rather than into the reservoir or the bin.
"""
from collections import Counter, deque
from dataclasses import replace
from importlib.metadata import PackageNotFoundError, version
from types import SimpleNamespace
from pathlib import Path
import hashlib
import importlib.util
import json
import platform
import sys
import numpy as np
from ase.io import read, write
from .config import SearchConfig
from .io import append_json, atomic_json, structure_hash
from .network import Registry
from . import internal
from .perturbations import FAMILIES, propose
from .chemistry import classify_transition, reaction_event_key
from .runner import (PathRecorder, confirm_minimum, continue_across_boundary,
                     descend_saddle, follow_min_mode, relax_source, run_trial)
from .state import classify, validate_atoms

# Distribution name paired with the module it installs.
_VERSION_PROBES = (("ase", "ase"), ("numpy", "numpy"), ("scipy", "scipy"),
                   ("torch", "torch"), ("mace-torch", "mace"), ("openmm", "openmm"),
                   ("openmm-ml", "openmmml"), ("openmm-torch", "openmmtorch"),
                   ("nnpops", "NNPOps"))


def versions():
    """Record the software stack as it actually is.

    Some conda packages ship no Python distribution metadata, so importlib.metadata
    alone reports None for openmm-ml, openmm-torch and nnpops while all three are
    installed and importable. It is not a property of conda packaging in general --
    openmm, numpy and pymbar-core do carry metadata -- which is exactly why metadata
    cannot be used as a presence test. A manifest that understates the stack cannot
    support reproduction, so fall back to module presence. find_spec is used rather
    than an import so that writing a manifest never pulls torch into the process.
    """
    result = {"python": sys.version, "platform": platform.platform()}
    for distribution, module in _VERSION_PROBES:
        try:
            result[distribution] = version(distribution)
            continue
        except PackageNotFoundError:
            pass
        try:
            found = importlib.util.find_spec(module) is not None
        except (ImportError, ValueError):
            found = False
        result[distribution] = "installed; no version metadata" if found else None
    return result


# sin of the angle below which an angle or dihedral coordinate has no usable
# gradient. Not a tuned number: it is the point where the formula starts dividing
# by a cross product that has gone to zero.
config_singular_sin = 1e-3


def _probe_tangent(atoms, probe):
    """The tangent of the coordinate the probe drives, from the proposal itself.

    The degraded fallback for when no adjacent probe frames were kept. Deliberately not the
    source-to-endpoint displacement: that sum contains the relaxation and whatever rotors
    turned along the way, so following it would aim the climb at a mixture rather than at
    the coordinate that was driven.

    `internal.gradient` returns dq/dx, a covector. The Cartesian displacement that changes q
    fastest under the mass metric is not that gradient but

        d_x  ~  M^-1 g_x,     so that   d_q = M^(1/2) d_x  ~  M^(-1/2) g_x,

    and the difference is not cosmetic. For an O-H bond the gradient has equal magnitude on
    both atoms; handing it over directly would have the walk weight the oxygen by sqrt(16)
    against the hydrogen's sqrt(1), aiming four times more strongly at the atom that barely
    moves. Frame differences are already displacements and are left alone -- only this
    gradient-derived fallback needs the inverse mass.
    """
    kind = FAMILIES[probe["family"]]["kind"] if probe.get("family") in FAMILIES else None
    if kind is None:
        return None, "no_family"
    indices = tuple(int(i) for i in probe["indices"])
    try:
        gradient = np.asarray(internal.gradient(atoms.positions, kind, indices),
                              dtype=float)
    except Exception:                                       # noqa: BLE001 - degraded path
        return None, "gradient_unavailable"
    if not np.isfinite(gradient).all():
        return None, "gradient_nonfinite"
    # An angle at 0 or pi, and a dihedral whose either end is collinear, have no defined
    # gradient. The formula does not raise there -- measured on a collinear dihedral it
    # returned a finite vector of norm 3.3e6 spread over every atom in the molecule, because
    # it divided by a cross product of 3.5e-18. Both checks below catch that, and either way
    # the answer is to fail closed: substituting the global lowest mode would aim the climb
    # at a different reaction.
    positions = atoms.positions
    legs = [positions[b] - positions[a]
            for a, b in zip(indices, indices[1:])]
    for first, second in zip(legs, legs[1:]):
        scale = np.linalg.norm(first) * np.linalg.norm(second)
        if scale < 1e-12:
            return None, "coordinate_degenerate"
        if np.linalg.norm(np.cross(first, second)) / scale < config_singular_sin:
            return None, "coordinate_singular"
    # An internal coordinate's gradient is supported on its own atoms. Anything else is
    # numerical debris from the same singularity.
    outside = [i for i in range(len(atoms)) if i not in indices]
    if outside:
        leak = float(np.linalg.norm(gradient[outside]))
        if leak > 1e-6 * max(float(np.linalg.norm(gradient[list(indices)])), 1e-30):
            return None, "gradient_leaks_outside_coordinate"
    tangent = gradient / atoms.get_masses()[:, None]
    if np.linalg.norm(tangent) < 1e-10:
        return None, "gradient_degenerate"
    return tangent, f"proposal_{kind}"


def rung_verdict(source_structure, endpoint, source_components,
                 endpoint_components, config):
    """Which of the three identity layers a continuation rung's endpoint differs in.

    Module level rather than inline in the search, so the classification that decides
    whether a rung counts as an observation can be tested against a known-degenerate
    system instead of only through a whole run.

    The chemical key alone is the wrong test, and the malonaldehyde run showed it
    immediately: its proton transfer is degenerate, so the continuation quenched onto the
    product -- q_PT went from +0.06 to -0.77, the proton demonstrably on the other oxygen --
    while the key was unchanged, and every rung was called "returned to the source". State
    identity and event identity are different questions, which is the whole reason the third
    layer exists. So the bond delta decides whether anything happened and the key decides
    whether it opened a node.

    Returns a dict with "verdict" in returned_to_source / reached_degenerate_product /
    reached_other_state, plus the bond delta and the transition class.
    """
    source_key = source_components["key"]
    endpoint_key = endpoint_components["key"]
    event = reaction_event_key(source_structure, endpoint, config.bond_scale,
                               config.active_atoms, config.automorphism_limit)
    # Components, not just digests: without them the fourth cell cannot be split and a real
    # E/Z isomerisation is reported as an unexplained key change.
    kind = classify_transition(source_key, endpoint_key, event["broken"], event["formed"],
                               source_components, endpoint_components)
    answer = {"chemical_key_matches_source": endpoint_key == source_key,
              "broken": event["broken"], "formed": event["formed"],
              "transition_class": kind, "new_state": endpoint_key != source_key}
    if kind == "conformational_transition":
        answer.update(verdict="returned_to_source",
                      note="same chemical state and no bond change, so this rung had not "
                           "entered another basin")
        return answer
    answer["verdict"] = ("reached_other_state" if answer["new_state"]
                         else "reached_degenerate_product")
    return answer


def _dividing_frame(directory, probe=None):
    """A seed frame on the ridge, and the tangent of the coordinate that drove it there.

    A Newton walk reaches the nearest stationary point, so the seed decides which one; and
    a mode-following walk needs to be told which coordinate to follow, or it climbs the
    softest thing in the molecule -- on 3-oxobutanal that was eight rotor barriers at -65 to
    -80 cm^-1 while the proton transfer sits at -3186.

    Tangent sources, in order of how local they are:

      topology switch   the difference between the switching frame and the one before it
      last probe step   the difference across the final frame of the driven phase
      proposal          the coordinate's own gradient, when no adjacent frames were kept

    Returns (atoms, description, tangent, tangent_source).
    """
    observations = directory / "observations.jsonl"
    trajectory = directory / "response.traj"
    if not observations.exists() or not trajectory.exists():
        return None, "no_trajectory", None, None
    rows = [json.loads(line) for line in observations.read_text().splitlines()]
    frames = read(str(trajectory), index=":")

    def tangent_between(later, earlier):
        if later is None or earlier is None:
            return None
        if later >= len(frames) or earlier >= len(frames) or later == earlier:
            return None
        step = frames[later].positions - frames[earlier].positions
        return step if np.linalg.norm(step) > 1e-8 else None

    switches = [row for row in rows
                if row["phase"] != "quench"
                and ((row.get("events") or {}).get("formed")
                     or (row.get("events") or {}).get("broken"))]
    if switches:
        row = min(switches, key=lambda r: r["frame"])
        index = row["frame"]
        previous = max((r["frame"] for r in rows if r["frame"] < index), default=None)
        tangent = tangent_between(index, previous)
        source = "topology_switch_difference"
        if tangent is None and probe is not None:
            tangent, source = _probe_tangent(frames[index], probe)
        if index < len(frames):
            return frames[index], f"topology_switch:{row['phase']}:{index}", tangent, source

    driven = [row["frame"] for row in rows if row["phase"] in ("perturbed", "pulse")]
    fallback = directory / "raw_endpoint.extxyz"
    if fallback.exists():
        atoms = read(str(fallback))
        tangent, source = None, None
        if len(driven) >= 2:
            tangent = tangent_between(driven[-1], driven[-2])
            source = "last_probe_step"
        if tangent is None and probe is not None:
            tangent, source = _probe_tangent(atoms, probe)
        return atoms, "raw_endpoint_no_switch_frame", tangent, source
    return None, "no_switch_frame", None, None



class ReactionSearch:
    """Serial reference implementation. Every reaction edge needs observed evidence."""

    def __init__(self, calculator_factory, config=None, calculator_metadata=None):
        self.factory = calculator_factory
        self.config = config or SearchConfig()
        self.calculator_metadata = calculator_metadata or {"factory": "python_callable"}

    def run(self, atoms, output):
        cfg = self.config
        validate_atoms(atoms, cfg.active_atoms)
        output = Path(output)
        if output.exists() and any(output.iterdir()):
            raise FileExistsError(f"Output must be a new or empty directory: {output}")
        for name in ("chemical", "ts_candidates"):
            (output / name).mkdir(parents=True, exist_ok=True)

        source_hash = hashlib.sha256()
        for path in sorted(Path(__file__).parent.glob("*.py")):
            source_hash.update(path.name.encode())
            source_hash.update(path.read_bytes())
        manifest = {
            "schema_version": 2, "status": "running", "config": cfg.to_dict(),
            "calculator": self.calculator_metadata, "versions": versions(),
            "input_structure_sha256": structure_hash(atoms),
            "implementation_sha256": source_hash.hexdigest(),
            "units": {"energy": "eV", "length": "Angstrom", "force": "eV/Angstrom",
                      "time": "fs", "mass": "amu", "angle": "radian"},
            "scope": "finite nonperiodic, fixed atom identities, single PES",
            "uncertainty_policy": "required" if cfg.max_uncertainty_eV_A else "not_required",
            "identity_policy": ("chemical node keyed on automorphism-invariant graph, "
                                "fragments and configurational stereochemistry; "
                                "conformation excluded"),
            "budget_policy": "chemical node slots and conformer slots are separate",
            "seed_policy": "seeded directions; common thermal seed across each amplitude scan",
        }
        manifest["experiment_sha256"] = hashlib.sha256(
            json.dumps(manifest, sort_keys=True).encode()).hexdigest()
        atomic_json(output / "config.json", manifest)
        write(str(output / "input.extxyz"), atoms)

        registry = Registry(cfg)
        records, ts_candidates = [], []
        saddle_searches = 0
        continuations = 0
        queue = deque()
        network = {"schema_version": 2, "status": "running", "chemical_nodes": [],
                   "reactions": [], "reaction_channels": [],
                   "ts_candidates": [], "ts_connected_channels": [],
                   "amplitude_scans": [], "attempts": 0,
                   "interpretation": ("Observed connectivity candidates between chemical "
                                      "macrostates; no TS, barrier or rate inferred")}

        def structure_path(state):
            return f"chemical/{state.id.replace('/', '_')}.extxyz"

        def publish():
            network["chemical_nodes"] = [{
                "id": node.id, "chemical_key": node.key, "depth": node.depth,
                "discovered_by": node.discovered_by,
                "fragments": node.components["fragments"],
                "key_components": {k: v for k, v in node.components.items()
                                   if k != "key"},
                "trials_spent": node.trials_spent,
                "reservoir": {"capacity": cfg.conformer_max_per_node,
                              "occupied": len(node.microstates),
                              "evicted": node.evicted},
                "microstates": [s.summary(structure_path(s)) for s in node.microstates],
            } for node in registry.nodes]
            network["ts_candidates"] = ts_candidates
            for published, node in zip(network["chemical_nodes"], registry.nodes):
                published["conformer_transitions"] = node.conformer_transitions
            # Derived here rather than stamped when a saddle is found: a directed edge
            # discovered later would otherwise miss saddle evidence recorded earlier.
            for edge in network["reactions"]:
                pair = {edge["source"], edge["target"]}
                support = [channel["ts"] for channel in network["ts_connected_channels"]
                           if set(channel["ends"]) == pair]
                if support:
                    edge["ts_support"] = support
            network["attempts"] = len(records)
            atomic_json(output / "network.json", network)

        try:
            calculator = self.factory()

            def physical_factory():
                return calculator

            root = relax_source(atoms, physical_factory, cfg)
            admission = registry.admit(root, None, None)
            if admission.outcome != "new_chemical_node":
                raise RuntimeError(f"source structure was not admitted: {admission.reason}")
            write(str(output / structure_path(admission.microstate)), root)
            queue.append((admission.node, admission.microstate))
            publish()

            def record_channel(source_structure, endpoint_structure, source_node,
                               target_node, evidence):
                """Record a reaction event, which is not the same thing as a new state.

                A chemical key quotients away every redundancy of the state, including the
                automorphism that makes a degenerate rearrangement end where it began. That
                is right for deciding whether a product is new, and blind to whether a
                reaction happened: a proton that moves from one oxygen to its equivalent
                partner returns the identical key over a bond set that changed. So the key
                decides the node and the bond delta decides the event, and a channel is
                allowed to be a self-loop.

                A self-loop opens no node and so spends no chemical budget, which keeps the
                two-layer network and its two budgets exactly as they were.
                """
                event = reaction_event_key(source_structure, endpoint_structure,
                                           cfg.bond_scale, cfg.active_atoms,
                                           cfg.automorphism_limit)
                kind = classify_transition(source_node.key, target_node.key,
                                           event["broken"], event["formed"],
                                           source_node.components,
                                           target_node.components)
                if kind == "conformational_transition":
                    return None, kind
                channel = next((c for c in network["reaction_channels"]
                                if c["event_key"] == event["key"]), None)
                if channel is None:
                    channel = {
                        "id": f"rc{len(network['reaction_channels']):04d}",
                        "event_key": event["key"], "class": kind,
                        "ends": sorted([source_node.id, target_node.id]),
                        "self_loop": source_node.id == target_node.id,
                        "broken": event["broken"], "formed": event["formed"],
                        "symmetry_enumerated": event["symmetry_enumerated"],
                        "undirected": True, "observed_by": [], "ts_support": [],
                        "meaning": ("bond-changing event, canonical under the endpoints' "
                                    "automorphisms and undirected; a self-loop is a "
                                    "degenerate reaction, not a conformational change")}
                    network["reaction_channels"].append(channel)
                for field_, value in evidence.items():
                    if value is not None and value not in channel[field_]:
                        channel[field_].append(value)
                return channel, kind

            def record_reaction(source_node, target_node, trial_id, record, probe,
                                mechanism="free_response"):
                """A directed observation: this source reached this target.

                `mechanism` is on the record because two different observations can produce
                the same edge and they are not equally direct. "free_response" is the trial
                itself -- probe, unbiased NVE, quench. "boundary_continuation" is the same
                probe advanced past the frame where the topology switched, quenched at a
                rung; the motion is the probe's own, the advance is not dynamics. Reading a
                continuation edge as a free-response observation would overstate it, and
                leaving the two indistinguishable would make that mistake unavoidable.
                """
                edge = next((e for e in network["reactions"]
                             if e["source"] == source_node.id
                             and e["target"] == target_node.id), None)
                if edge is None:
                    edge = {"id": f"r{len(network['reactions']):04d}",
                            "source": source_node.id, "target": target_node.id,
                            "observed_only": True, "attempts": [],
                            "classifications": [], "best_by_family": {},
                            "source_microstates": [], "mechanisms": [],
                            "continuations": []}
                    network["reactions"].append(edge)
                if mechanism not in edge.setdefault("mechanisms", []):
                    edge["mechanisms"].append(mechanism)
                if mechanism == "free_response":
                    edge["attempts"].append(trial_id)
                else:
                    # Not an attempt. `attempts` means "this trial, on its own, went from
                    # source to target", and the seed trial of a continuation did not: it
                    # stopped at the boundary, which is why the continuation was run. Listing
                    # it there would make the seed trial claim an observation it never made,
                    # and every consumer that reads a trial's own record to check the edge
                    # would find the two disagreeing.
                    edge.setdefault("continuations", []).append(
                        {"seed_attempt": trial_id, "mechanism": mechanism})
                if record["classification"] not in edge["classifications"]:
                    edge["classifications"].append(record["classification"])
                if record["source_microstate"] not in edge["source_microstates"]:
                    edge["source_microstates"].append(record["source_microstate"])
                best = edge["best_by_family"].get(probe.family)
                if best is None or probe.amplitude < best["amplitude"]:
                    edge["best_by_family"][probe.family] = {
                        "amplitude": probe.amplitude, "unit": probe.unit,
                        "attempt": trial_id,
                        "meaning": "minimum observed success, not a barrier"}
                record["reaction"] = edge["id"]

            def register_saddle(structure, diagnostics, attempt_id, reached_from, node,
                                seed, found_by, quench_rounds=()):
                """Record a first-order saddle and resolve both of its endpoints.

                Shared by the two ways one is reached. A quench that stalls on the ridge is
                a soft-coordinate coincidence and produced every saddle the ethanol run
                found; climbing to one from a probe frame is what works when the reaction
                coordinate is stiff, which is every case where the ridge is sharp. How it
                was reached is recorded because the two do not have the same reliability.
                """
                candidate_id = f"ts{len(ts_candidates):04d}"
                relative = f"ts_candidates/{candidate_id}.extxyz"
                write(str(output / relative), structure)
                mode = diagnostics.pop("unstable_mode", None)
                if mode is not None:
                    atomic_json(output / f"ts_candidates/{candidate_id}_mode.json",
                                {"unstable_mode": mode,
                                 "meaning": "mass-unweighted reaction coordinate"})
                # The polish already reported which coordinate carries the negative
                # curvature, so the reaction coordinate has a chemically readable name
                # without projecting the eigenvector. A climbed saddle has no polish
                # report of its own, so the list is simply empty there.
                negative = [tuple(report["indices"])
                            for round_ in quench_rounds
                            for report in round_["modes"]
                            if report["tier"] == "negative_curvature"]
                ts_entry = {
                    "id": candidate_id, "attempt": attempt_id,
                    "structure": relative, "reached_from": reached_from,
                    "saddle_order": 1, "found_by": found_by,
                    # Only the a -> 0 path produces genuine eigenvalues, so only it may
                    # quote a frequency. A band-regressed response curvature is a different
                    # quantity and is reported as one.
                    "imaginary_wavenumbers_icm":
                        diagnostics.get("imaginary_wavenumbers_icm"),
                    "curvature_source": (diagnostics.get("provenance") or {}).get("source"),
                    "trivial_mode_floor_worst_residual":
                        diagnostics.get("trivial_mode_floor_worst_residual"),
                    "response": diagnostics.get("response"),
                    "energy_eV": float(structure.get_potential_energy()),
                    "negatively_curved_torsions": [list(i) for i in negative],
                    "status": "unrefined first-order saddle candidate"}
                ts_candidates.append(ts_entry)
                if (cfg.irc_enabled and mode is not None
                        and len(ts_candidates) <= cfg.irc_max_per_run):
                    connect_saddle(ts_entry, structure, mode, node, attempt_id, seed)
                return ts_entry

            def connect_saddle(ts_entry, saddle, mode, node, trial_id, seed):
                """Name the two minima a first-order saddle joins, by descending both sides.

                A saddle connection is undirected structural evidence and is kept apart
                from the observed directed edges: seeing A reach B does not license B to A,
                and a shared saddle is not an observation of either direction. Endpoints in
                one chemical node are a conformer interconversion, which belongs inside that
                node rather than in the reaction network.
                """
                ends, admissions = [], []
                for sign in (1, -1):
                    entry = {"sign": sign}
                    # Local connectivity, established one ridge at a time. The step off
                    # the saddle stays small, because a large displacement can cross several
                    # dividing surfaces and land in a basin that is not adjacent -- and then
                    # "these two minima are joined through this saddle" would be a claim the
                    # walk never checked. Growing the original displacement until something
                    # is a minimum was exactly that mistake.
                    #
                    # When the descent lands on another saddle instead, the walk continues
                    # from THAT saddle along its own unstable mode. On 3-oxobutanal the soft
                    # torsions put several ridges within one step of each other and every
                    # descent from all eight candidates stopped on a further first-order
                    # saddle; relaying through them keeps every hop local. How many hops it
                    # took is recorded, because a chain of them is weaker evidence than one.
                    descent = PathRecorder(output / "paths" / ts_entry["id"], cfg, "descent")
                    descent.context = {"sign": sign}
                    endpoint, checked, relays = descend_saddle(
                        saddle, mode, physical_factory, cfg, sign,
                        seed=seed + 7919 * sign, recorder=descent)
                    ok = endpoint is not None
                    entry["relays"] = relays
                    entry["relay_meaning"] = (
                        "each hop is a local descent from one saddle; a chain of them is "
                        "weaker connectivity evidence than a single step")
                    admission = None
                    if not ok:
                        if checked.get("saddle_order") is None:
                            entry.update(status="quench_failed",
                                         reason=checked.get("reason"))
                        else:
                            entry.update(status="not_a_minimum",
                                         saddle_order=checked.get("saddle_order"))
                    else:
                        admission = registry.admit(endpoint, node, trial_id)
                        entry["status"] = admission.outcome
                        if admission.microstate is None:
                            entry["reason"] = admission.reason
                            admission = None
                        else:
                            entry.update(chemical_node=admission.node.id,
                                         microstate=admission.microstate.id,
                                         energy_eV=admission.microstate.energy_eV)
                            path = output / structure_path(admission.microstate)
                            if not path.exists():
                                write(str(path), admission.microstate.structure)
                            if (admission.outcome in ("new_chemical_node",
                                                      "new_microstate")
                                    and admission.node.depth < cfg.max_depth):
                                queue.append((admission.node, admission.microstate))
                    ends.append(entry)
                    admissions.append(admission)
                ends_note = ("steepest-descent endpoints from the saddle; no path, "
                             "barrier or direction is implied")
                ts_entry.update(connects=ends, connects_meaning=ends_note)
                if not all(admissions):
                    return
                left, right = admissions
                # Whether the two ends are one node decides which network they belong to;
                # whether any bond changed decides whether this was a reaction at all.
                # Ethanol's methyl barrier answers no to both and stays a conformer
                # transition; a degenerate proton transfer answers no to the first and yes
                # to the second, and is a reaction with a self-loop.
                channel, kind = record_channel(left.microstate.structure,
                                               right.microstate.structure,
                                               left.node, right.node,
                                               {"ts_support": ts_entry["id"]})
                ts_entry["transition_class"] = kind
                if channel is not None:
                    ts_entry["reaction_channel"] = channel["id"]
                if kind == "conformational_transition":
                    left.node.conformer_transitions.append({
                        "ts": ts_entry["id"],
                        "microstates": [left.microstate.id, right.microstate.id],
                        "self_connected": left.microstate is right.microstate,
                        "negatively_curved_torsions":
                            ts_entry["negatively_curved_torsions"],
                        "meaning": "conformer interconversion through a saddle, not a reaction"})
                if left.node is not right.node:
                    network["ts_connected_channels"].append({
                        "ts": ts_entry["id"], "undirected": True,
                        "ends": [left.node.id, right.node.id],
                        "end_microstates": [left.microstate.id, right.microstate.id],
                        "transition_class": kind,
                        "meaning": ("chemical states joined through a saddle; undirected, "
                                    "and no barrier is implied")})


            def seek_saddle(node, state, observed):
                """Climb to the saddle of a reaction this scan already observed.

                Seeds are ranked by how far they are from a basin rather than by how close
                they are to the dividing surface, because a Newton walk converges to the
                nearest stationary point and a frame that barely crossed is nearest to the
                minimum it came from. A transient crossing -- the topology switched and did
                not hold -- is the frame that sits on the ridge.
                """
                nonlocal saddle_searches
                if not cfg.saddle_search_enabled:
                    return {"attempted": False, "reason": "disabled"}
                reacted = [r for r in observed if r.get("reaction_channel")]
                boundaries = [r for r in observed if r["status"] == "ts_candidate"]
                crossings = [r for r in observed
                             if (r.get("crossing") or {}).get("observed")]
                if not (reacted or boundaries or crossings):
                    return {"attempted": False, "reason": "nothing_left_the_basin"}
                if saddle_searches >= cfg.saddle_search_per_run:
                    return {"attempted": False, "reason": "saddle_search_budget"}

                def rank(record):
                    return (0 if (record.get("crossing") or {}).get("observed") else
                            1 if record.get("reaction_channel") else 2,
                            record["probe"]["amplitude"])
                seeds = sorted(crossings + reacted + boundaries, key=rank)
                # requested_channel is the search's intent and only exists when a channel
                # was actually observed in this scan. A scan that produced only basin
                # boundaries has nothing to request, and inventing one from the direction
                # that was probed would manufacture the intent the metric is supposed to
                # measure against.
                requested = reacted[0]["reaction_channel"] if reacted else None
                report = {"attempted": True, "channel": requested,
                          "probed": list(observed[0]["probe"]["indices"]),
                          "seeds": []}
                for candidate in seeds[:cfg.saddle_search_seeds]:
                    frame, origin, tangent, tangent_source = _dividing_frame(
                        output / "trials" / candidate["attempt_id"],
                        probe=candidate.get("probe"))
                    if frame is None:
                        continue
                    if tangent is None:
                        # Refusing rather than climbing the softest mode. Without a tangent
                        # the walk answers a different question, and on this molecule the
                        # answer is a methyl rotor.
                        report["seeds"].append(
                            {"seed_attempt": candidate["attempt_id"], "seed_frame": origin,
                             "reason": "missing_target_direction",
                             "amplitude": candidate["probe"]["amplitude"],
                             "transient": bool((candidate.get("crossing") or {}).get("observed"))})
                        continue
                    saddle_searches += 1
                    climb = PathRecorder(output / "paths" / f"climb{saddle_searches:03d}",
                                         cfg, "climb")
                    climb.context = {"seed_attempt": candidate["attempt_id"],
                                     "direction_source": tangent_source}
                    structure, diagnostics = follow_min_mode(
                        frame, physical_factory, cfg, direction=tangent,
                        ceiling=state.energy_eV + cfg.saddle_search_max_rise_eV,
                        reference=state.energy_eV, recorder=climb)
                    entry = {"seed_attempt": candidate["attempt_id"],
                             "seed_frame": origin,
                             "direction_source": tangent_source,
                             "target_overlap": diagnostics.get("target_overlap"),
                             "selected_mode_index": diagnostics.get("selected_mode_index"),
                             "amplitude": candidate["probe"]["amplitude"],
                             "transient": bool((candidate.get("crossing") or {}).get("observed")),
                             "steps": diagnostics.get("steps"),
                             "saddle_order": diagnostics.get("saddle_order"),
                             "reason": diagnostics.get("reason"),
                             "initial_rise_eV": diagnostics.get("initial_rise_eV"),
                             "max_rise_eV": diagnostics.get("max_rise_eV"),
                             "backtracks": diagnostics.get("backtracks"),
                             "ceiling_hit_step": diagnostics.get("ceiling_hit_step")}
                    if structure is not None:
                        ts_entry = register_saddle(
                            structure, diagnostics, candidate["attempt_id"], state.id, node,
                            cfg.seed + saddle_searches, found_by="min_mode_following")
                        entry["ts_candidate"] = ts_entry["id"]
                        # The channel this saddle supports is decided by descending both
                        # of its sides, which register_saddle already did -- not by which
                        # scan happened to launch the climb. A climb is free to walk to a
                        # different saddle than the one that was wanted: seeded from the
                        # malonaldehyde dividing surface it reached the hydroxyl rotamer
                        # barrier, and seeded from a probe frame in this run it reached an
                        # O-H homolysis at 2.34 eV. Attaching either to the proton transfer
                        # because that was the reason for looking would be a fabrication.
                        entry["supports"] = [c["id"] for c in network["reaction_channels"]
                                             if ts_entry["id"] in c["ts_support"]]
                        entry["requested_channel"] = requested
                        # Three different things, kept apart. requested_channel is the
                        # search's intent; supports is what descending both sides actually
                        # established; supports_requested is a hit rate for the seeding
                        # strategy and must never stand in for either.
                        #
                        # It is null, not false, when one side of the descent failed. False
                        # is a claim -- the endpoint evidence is complete and it supports
                        # some other channel -- and asserting it on incomplete evidence
                        # would turn a missing measurement into a negative result.
                        ends = ts_entry.get("connects") or []
                        complete = (len(ends) == 2
                                    and all(end.get("microstate") for end in ends))
                        entry["attribution_complete"] = complete
                        entry["supports_requested"] = (
                            (requested in entry["supports"])
                            if (complete and requested is not None) else None)
                        report["seeds"].append(entry)
                        report["found"] = ts_entry["id"]
                        return report
                    report["seeds"].append(entry)
                    if saddle_searches >= cfg.saddle_search_per_run:
                        break
                report["found"] = None
                return report

            def continue_boundary(node, state, observed):
                """Follow the delivered probe past the topology switch until it settles.

                The question this answers is the one the saddle search does not: did the
                motion that was actually delivered get to the product? Round three of the
                3-oxobutanal benchmark had three saddles supporting the proton-transfer
                channel and zero trials in which A reached B, and the two facts are
                consistent -- the probe stopped on the ridge and the channel was named by
                descending from the saddles, which is undirected structural evidence.

                Seeds are the frames that reached a boundary, transient crossings first,
                exactly as the saddle search ranks them: those are the frames nearest the
                dividing surface, which is where a continuation has something left to do.

                Endpoints in the source's own chemical state are classified read-only, not
                admitted. Admitting them would spend conformer slots and shift the
                admission order for everything after, and the question here is whether the
                path leaves the state, not how many ways it can come back to it.
                """
                nonlocal continuations
                if not cfg.boundary_continuation_enabled:
                    return {"attempted": False, "reason": "disabled"}
                if continuations >= cfg.boundary_continuation_per_run:
                    return {"attempted": False, "reason": "continuation_budget"}
                crossings = [r for r in observed
                             if (r.get("crossing") or {}).get("observed")]
                boundaries = [r for r in observed if r["status"] == "ts_candidate"]
                seeds = crossings + boundaries
                if not seeds:
                    return {"attempted": False, "reason": "no_boundary_frame"}
                report = {"attempted": True, "seeds": [], "reached": None,
                          "meaning": ("the delivered probe continued past the topology "
                                      "switch; a reached state here is a directed "
                                      "observation of the probe, not of free dynamics")}
                for candidate in seeds[:cfg.boundary_continuation_seeds]:
                    frame, origin, tangent, tangent_source = _dividing_frame(
                        output / "trials" / candidate["attempt_id"],
                        probe=candidate.get("probe"))
                    entry = {"seed_attempt": candidate["attempt_id"], "seed_frame": origin,
                             "direction_source": tangent_source,
                             "amplitude": candidate["probe"]["amplitude"],
                             "transient": bool((candidate.get("crossing") or {})
                                               .get("observed"))}
                    if frame is None or tangent is None:
                        entry["reason"] = ("no_seed_frame" if frame is None
                                           else "missing_target_direction")
                        report["seeds"].append(entry)
                        continue
                    continuations += 1
                    recorder = PathRecorder(
                        output / "paths" / f"continuation{continuations:03d}", cfg,
                        "continuation")
                    recorder.context = {"seed_attempt": candidate["attempt_id"],
                                        "direction_source": tangent_source}
                    def name_endpoint(endpoint, confirmed, checked):
                        """Pure. Deliberately so.

                        An earlier version stashed the endpoint in an enclosing dict when the
                        verdict was not a return, and the search admitted whatever was in
                        that dict. But this callback is also run on BOTH descents of a rung
                        that stalled on a saddle, and there one side can name the product
                        while the other fails -- which is precisely the case the rung reports
                        as "stopped_on_saddle, not evidence of arrival". Measured on P2 round
                        five: two of the three continuations on the directed edge came from
                        rungs whose own verdict said they had not arrived. A callback that
                        the caller runs speculatively must not be able to record anything.
                        """
                        return rung_verdict(state.structure, endpoint,
                                            node.components,
                                            registry.key_of(endpoint), cfg)

                    arrived, walk = continue_across_boundary(
                        frame, tangent, physical_factory, cfg, verdict=name_endpoint,
                        recorder=recorder)
                    entry["continuation"] = walk
                    # Only what the walk itself hands back, and only when its own outcome is
                    # an arrival. Nothing else may become an admission.
                    endpoint = arrived
                    reached = ({} if endpoint is None
                               else rung_verdict(state.structure, endpoint,
                                                 node.components,
                                                 registry.key_of(endpoint), cfg))
                    if endpoint is not None and not reached["new_state"]:
                        # A degenerate rearrangement: the bond set changed and the key did
                        # not. It is a reaction event with a self-loop channel and it opens
                        # no node, so the event is recorded from the endpoint structure
                        # directly and nothing is admitted. Admitting here would spend a
                        # conformer slot to record an event that has no new state in it.
                        channel, kind = record_channel(
                            state.structure, endpoint, node, node,
                            {"observed_by": candidate["attempt_id"]})
                        entry["transition_class"] = kind
                        entry["self_loop"] = True
                        if channel is not None:
                            entry["reaction_channel"] = channel["id"]
                        report["reached"] = node.id
                    elif endpoint is not None:
                        admission = registry.admit(endpoint, node,
                                                   candidate["attempt_id"])
                        entry["admission"] = admission.outcome
                        if admission.microstate is None:
                            entry["admission_reason"] = admission.reason
                        else:
                            target_node = admission.node
                            target_state = admission.microstate
                            path = output / structure_path(target_state)
                            if not path.exists():
                                write(str(path), target_state.structure)
                            if (admission.outcome in ("new_chemical_node", "new_microstate")
                                    and target_node.depth < cfg.max_depth):
                                queue.append((target_node, target_state))
                            channel, kind = record_channel(
                                state.structure, target_state.structure, node, target_node,
                                {"observed_by": candidate["attempt_id"]})
                            entry["transition_class"] = kind
                            entry["target"] = target_node.id
                            entry["target_microstate"] = target_state.id
                            if channel is not None:
                                entry["reaction_channel"] = channel["id"]
                            if target_node.id != node.id:
                                # Run the same classifier the trials run, on this
                                # endpoint. The label used to be hardcoded "reactive",
                                # which is a claim rather than a measurement and is simply
                                # wrong for an E/Z arrival: no bond changes there, so
                                # classify() returns a rearrangement or an elastic label.
                                label, _ = classify(state.structure, endpoint, cfg)
                                synthetic = {"classification": label,
                                             "source_microstate": state.id}
                                # The probe that produced the seed, read back off the
                                # record: the continuation carried that same probe further,
                                # so the family and amplitude on the edge are its own.
                                seed_probe = SimpleNamespace(
                                    family=candidate["probe"]["family"],
                                    amplitude=candidate["probe"]["amplitude"],
                                    unit=candidate["probe"]["unit"])
                                record_reaction(node, target_node,
                                                candidate["attempt_id"], synthetic,
                                                seed_probe,
                                                mechanism="boundary_continuation")
                                entry["reaction"] = synthetic.get("reaction")
                            report["reached"] = entry.get("target")
                    report["seeds"].append(entry)
                    if report["reached"] is not None:
                        break
                    if continuations >= cfg.boundary_continuation_per_run:
                        break
                return report

            def execute(node, state, probe, stage, parents=None):
                trial_id = f"t{len(records):06d}"
                outcome = run_trial(state.structure, probe, physical_factory, cfg,
                                    output / "trials" / trial_id, trial_id)
                record = outcome.record
                record.update(source=node.id, source_microstate=state.id, target=None,
                              target_microstate=None, stage=stage,
                              refinement_parents=parents, admission=None)
                if outcome.endpoint is not None:
                    seed = int(np.random.SeedSequence(
                        [cfg.seed, len(registry.nodes), len(records)]).generate_state(1)[0])
                    confirmed, diagnostics = confirm_minimum(
                        outcome.endpoint, physical_factory, cfg, seed)
                    record["minimum_confirmation"] = diagnostics
                    if not confirmed:
                        order = diagnostics.get("saddle_order")
                        if order == 1:
                            ts_entry = register_saddle(
                                outcome.endpoint, diagnostics, trial_id, state.id, node,
                                seed, found_by="quench_stalled",
                                quench_rounds=record.get("quench", {}).get("rounds", []))
                            record.update(status="ts_candidate", admission="ts_candidate")
                            record["ts_candidate"] = ts_entry["id"]
                        else:
                            record.update(status="not_a_minimum", admission="rejected")
                            record["failure"] = {
                                "code": "not_a_minimum",
                                "saddle_order": order,
                                "message": "stationary point is not a local minimum"}
                    elif record["status"] != "completed":
                        # A transient crossing hands over its endpoint so the saddle check
                        # can run on it, but its product claim is not trusted and must not
                        # enter the network. The endpoint being a minimum here just means
                        # the trial fell back to a basin.
                        record["admission"] = "withheld_untrusted_topology"
                    else:
                        admission = registry.admit(outcome.endpoint, node, trial_id)
                        record["admission"] = admission.outcome
                        record["admission_detail"] = admission.detail or None
                        if admission.outcome == "rejected":
                            record["network_exclusion"] = admission.reason
                        else:
                            target_node, target_state = admission.node, admission.microstate
                            record["target"] = target_node.id
                            record["target_microstate"] = target_state.id
                            path = output / structure_path(target_state)
                            if not path.exists():
                                write(str(path), target_state.structure)
                            if admission.outcome in ("new_chemical_node", "new_microstate"):
                                if target_node.depth < cfg.max_depth:
                                    queue.append((target_node, target_state))
                            channel, kind = record_channel(
                                state.structure, target_state.structure, node, target_node,
                                {"observed_by": trial_id})
                            record["transition_class"] = kind
                            if channel is not None:
                                record["reaction_channel"] = channel["id"]
                            if target_node.id != node.id:
                                record_reaction(node, target_node, trial_id, record, probe)
                registry.note_trial(node, state, record.get("target") or record["status"])
                records.append(record)
                append_json(output / "attempts.jsonl", record)
                atomic_json(output / "trials" / trial_id / "result.json", record)
                publish()
                return record

            while queue and len(records) < cfg.max_trials:
                node, state = queue.popleft()
                allowed, blocked = registry.may_search(node, state)
                if not allowed:
                    network.setdefault("budget_blocks", []).append(
                        {"node": node.id, "microstate": state.id, "reason": blocked})
                    continue
                direction_seed = int(np.random.SeedSequence(
                    [cfg.seed, node.id.encode()[-2] if node.id else 0,
                     len(node.microstates), state.trials_spent]).generate_state(1)[0])
                proposal = {}
                for probes in propose(state.structure, cfg, direction_seed,
                                      report=proposal):
                    if len(records) >= cfg.max_trials:
                        break
                    allowed, blocked = registry.may_search(node, state)
                    if not allowed:
                        break
                    coarse = []
                    for probe in probes:
                        if len(records) >= cfg.max_trials:
                            break
                        allowed, _ = registry.may_search(node, state)
                        if not allowed:
                            break
                        coarse.append(execute(node, state, probe, "coarse"))
                    if not coarse:
                        break
                    scan = {"source": node.id, "source_microstate": state.id,
                            "family": probes[0].family, "indices": list(probes[0].indices),
                            "sign": probes[0].sign, "seed": probes[0].seed,
                            "unit": probes[0].unit, "local_brackets": []}
                    observed = list(coarse)

                    # Refine only where a return to the origin microstate sits next to a
                    # departure from it. The discriminator is the microstate, not the
                    # chemical node, so a conformational transition is bracketed on the
                    # same footing as a chemical one. Rejected and unconverged trials never
                    # stand in for a return.
                    #
                    # A transient crossing counts as a departure. Elastic, transient and
                    # confirmed are ordered information about where the dividing surface
                    # lies, not two failures and a success -- and the transient one is the
                    # most informative of the three, being the amplitude closest to that
                    # surface. On malonaldehyde it was the only probe that reached the
                    # ridge at all, and treating it as a failure closed the bracket that
                    # would have found the saddle.
                    def returned(record):
                        return (record["status"] == "completed"
                                and record.get("target_microstate") == state.id)

                    def departed(record):
                        if record["status"] == "completed":
                            reached = record.get("target_microstate")
                            return reached is not None and reached != state.id
                        if record["status"] == "ts_candidate":
                            # A quench that stopped on a ridge reached SOME basin boundary.
                            # That is enough to bracket an amplitude -- returning stopped
                            # happening between here and the last amplitude that returned --
                            # and it is not enough to claim the probe crossed the coordinate
                            # it was aimed at. On 3-oxobutanal every O-H amplitude ended on
                            # a rotor barrier at -65 to -80 cm^-1 while the proton transfer
                            # sits at -3186; calling that a departure along the O-H would
                            # have been evidence for a reaction that never happened.
                            return True
                        return bool((record.get("crossing") or {}).get("observed"))

                    for low_record, high_record in zip(coarse, coarse[1:]):
                        if not (returned(low_record) and departed(high_record)):
                            continue
                        low = low_record["probe"]["amplitude"]
                        high = high_record["probe"]["amplitude"]
                        reached = high_record.get("target_microstate")
                        bracket = {"target_microstate": reached, "low": low, "high": high,
                                   "high_kind": {"completed": "confirmed",
                                                 "ts_candidate": "boundary_unattributed"}
                                                .get(high_record["status"],
                                                     "transient_crossing"),
                                   "high_kind_meaning": (
                                       "boundary_unattributed means a basin boundary was "
                                       "reached, not that the probed coordinate was "
                                       "crossed"),
                                   "meaning": "local observed window; no monotonicity assumed",
                                   "interrupted_by": None}
                        for _ in range(cfg.refinement_steps):
                            if len(records) >= cfg.max_trials:
                                break
                            allowed, _ = registry.may_search(node, state)
                            if not allowed:
                                break
                            middle = (low + high) / 2
                            refined = execute(node, state,
                                              replace(probes[0], amplitude=middle),
                                              "refine", [low_record["attempt_id"],
                                                         high_record["attempt_id"]])
                            observed.append(refined)
                            if departed(refined):
                                high, high_record = middle, refined
                                if refined["status"] != "completed":
                                    bracket["high_kind"] = (
                                        "boundary_unattributed"
                                        if refined["status"] == "ts_candidate"
                                        else "transient_crossing")
                            elif returned(refined):
                                low, low_record = middle, refined
                            else:
                                bracket["interrupted_by"] = refined["attempt_id"]
                                break
                        bracket.update(low=low, high=high)
                        scan["local_brackets"].append(bracket)

                    observed.sort(key=lambda r: r["probe"]["amplitude"])
                    scan["samples"] = [{"amplitude": r["probe"]["amplitude"],
                                        "attempt": r["attempt_id"], "status": r["status"],
                                        "admission": r["admission"], "target": r["target"],
                                        "target_microstate": r["target_microstate"]}
                                       for r in observed]
                    # A later success after an intervening return is a second window,
                    # not evidence that the first one was wrong.
                    departed = False
                    scan["nonmonotonic_return_observed"] = False
                    for entry in observed:
                        if entry["status"] != "completed" or entry["target_microstate"] is None:
                            continue
                        if entry["target_microstate"] != state.id:
                            departed = True
                        elif departed:
                            scan["nonmonotonic_return_observed"] = True
                    # A scan that produced a reaction event carries the one geometry a
                    # saddle can be reached from. Quenching to a saddle only happens when
                    # the reaction coordinate is soft -- on malonaldehyde every amplitude
                    # in the scan, including the one bracketing the dividing surface,
                    # quenched into a basin while the saddle sat at lambda_1 = -37.3.
                    # Climbing does not care how sharp the ridge is; it cares where it
                    # starts, and a Newton walk reaches the nearest stationary point. Of
                    # the frames this scan produced, only the transient crossing is nearer
                    # to the saddle than to a basin: measured, the 0.25 and 0.275 frames
                    # climbed back to a minimum and the 0.35 transient one reached the
                    # saddle in eight steps.
                    # Both interventions are named, so an improvement cannot be
                    # credited to the wrong one: the proposal priority decides which
                    # directions are offered, the pair-approach feasibility test decides
                    # which are worth offering at all.
                    scan["proposal"] = {
                        "priority": probes[0].priority,
                        "priority_reason": probes[0].priority_reason,
                        "budget": proposal}
                    scan["saddle_search"] = seek_saddle(node, state, observed)
                    # A third named intervention, kept separate from the other two for the
                    # same reason: the saddle search answers "is there a saddle on this
                    # coordinate", the continuation answers "does the delivered probe reach
                    # the product", and an improvement in either must be attributable.
                    scan["boundary_continuation"] = continue_boundary(node, state, observed)
                    network["amplitude_scans"].append(scan)
                    publish()

            network.update(
                status="completed",
                termination=("trial_budget" if len(records) >= cfg.max_trials
                             else "configured_scope_exhausted"),
                outcome_counts=dict(Counter(r["status"] for r in records)),
                classification_counts=dict(Counter(r["classification"] for r in records
                                                   if r["classification"] is not None)),
                admission_counts=dict(Counter(r["admission"] for r in records
                                              if r["admission"] is not None)),
                budgets={"chemical_nodes_used": len(registry.nodes),
                         "chemical_nodes_available": cfg.chemical_max_nodes,
                         "microstates_total": sum(len(n.microstates)
                                                  for n in registry.nodes),
                         "ts_candidates": len(ts_candidates)})
            manifest["status"] = "completed"
            publish()
            atomic_json(output / "config.json", manifest)
            return network
        except BaseException as exc:
            manifest.update(status="failed",
                            failure={"type": type(exc).__name__, "message": str(exc)})
            network["status"] = "partial" if registry.nodes else "failed"
            publish()
            atomic_json(output / "config.json", manifest)
            raise
