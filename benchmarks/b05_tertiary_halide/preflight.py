"""b05 preflight: what has to be known before criteria can be written.

Nothing here is an acceptance result. This answers "can this grid be run, and could it
answer its own question", which is a different question from "did PRRS find anything".

Gates 1-3 are b04's, re-measured on this structure rather than assumed: they are about
this grid's charge and this grid's identity layer, and a gate inherited by reference is
a gate nobody ran. Gates 6 and 7 are new and they are the ones that matter here -- a
grid whose second channel is unreachable, or whose outcomes are indistinguishable, would
report a clean null that means nothing.

    python benchmarks/b05_tertiary_halide/preflight.py            # 1-3, 6, 7: CPU, seconds
    python benchmarks/b05_tertiary_halide/preflight.py --models   # also 4, 5, 8: needs a GPU
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).parent))

from ase import Atoms  # noqa: E402

from build import (  # noqa: E402
    CENTRAL_C,
    CHLORINE,
    OXYGEN,
    build,
    elimination_products,
    substitution_products,
)

results = []


def gate(name, ok, detail, required=True):
    results.append((name, ok, required))
    mark = "PASS" if ok else ("FAIL" if required else "WARN")
    print(f"  [{mark}] {name}")
    for line in detail.strip().splitlines():
        print(f"         {line}")
    print()


def config():
    text = Path(__file__).with_name("config.json").read_text()
    return json.loads(
        "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
    )


print(f"b05 preflight  ·  tree {ROOT}")
print("=" * 78)

from prrs.chemistry import chemical_key  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402

CHARGE = dict(charge_sensitive=True, charge=-1, multiplicity=1)
BLOCKING = ("no_assignment_satisfies_valences", "degree_exceeds_valence")
cfg_json = config()

# ---------------------------------------------------------------- gate 1
hydroxide = Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.964]])
hydroxide_key = chemical_key(hydroxide, 1.2, **CHARGE)
reactants = chemical_key(build(), 1.2, **CHARGE)
vetoes = bool(cfg_json.get("closed_shell_only", True))
gate(
    "1  a bare hydroxide ion is not vetoed by the identity layer",
    (not vetoes)
    or (
        hydroxide_key.get("stereo_unresolved") not in BLOCKING
        and reactants.get("stereo_unresolved") not in BLOCKING
    ),
    f"""OH- alone            stereo_unresolved = {hydroxide_key.get("stereo_unresolved")!r}
tBuCl + OH-          stereo_unresolved = {reactants.get("stereo_unresolved")!r}
fragments            {reactants["fragments"]}
closed_shell_only    {vetoes}  (from config.json)
Same shape as b04's gate 1 and the same resolution: O's neutral valence is 2 and a
hydroxide oxygen has degree 1, so no NEUTRAL assignment exists. That is a statement
about the Lewis annotator, not about the species -- OH- is closed-shell. The refusal
would come from the model's closed_shell_only declaration, and POLAR-1 declares False
because its charge response is measured. What the identity layer says is recorded as an
annotation either way.""",
)

# ---------------------------------------------------------------- gate 2
try:
    cfg = SearchConfig(**{k: v for k, v in cfg_json.items()})
    built, why = (
        True,
        (
            f"charge_sensitive={cfg.charge_sensitive} total_charge={cfg.total_charge} "
            f"multiplicity={cfg.multiplicity}"
        ),
    )
except Exception as exc:  # noqa: BLE001 - reported below
    cfg, built, why = None, False, f"{type(exc).__name__}: {exc}"
gate(
    "2  the frozen-draft config is constructible as written",
    built,
    f"""{why}
b04's gate 2 built a SearchConfig by hand; this one builds it from config.json, so a
typo in the file that would only surface at run time surfaces here instead.""",
)

# ---------------------------------------------------------------- gate 3
from prrs.search import _stamped_with_charge  # noqa: E402

probe = build()
probe.info.clear()
stamped = _stamped_with_charge(probe, cfg or SearchConfig(**CHARGE, closed_shell_only=False))
gate(
    "3  the total charge reaches the calculator",
    stamped.info.get("charge") == -1 and stamped.info.get("spin") == 1,
    f"""stamped info: charge={stamped.info.get("charge")} spin={stamped.info.get("spin")}
Worst form of this failure is silent: the config declares -1, the model evaluates the
neutral, and nothing raises. _stamped_with_charge refuses to start when the structure
and the config disagree.""",
)

# ---------------------------------------------------------------- gate 6
# The grid's whole point is that TWO channels compete, so both have to be reachable.
# Under the 2026-09-09 direction contract this gate has to separate two reasons that an
# earlier version of it merged, because they are different findings with different fixes:
#
#   basis_uncovered          no candidate direction anywhere can move this coordinate.
#                            A permanent gap: the input or the families have to change.
#   candidate_not_selected   a candidate exists, the direction budget did not take it.
#                            A budget artifact: exhaustive mode or a bigger budget fixes it.
#
# The full basis is read by proposing in `exhaustive` mode, which by construction selects
# every symmetry-reduced candidate; the grid's own mode gives what a run would actually
# reach. Asked in terms of config.json's OWN tracked_coordinates, never particular atom
# pairs -- b04's SN2 was entered by `stretch [0,1]`, not by the nucleophile's approach --
# and tracked_coordinates are an audit input here, never a selection input (contract 1).
from dataclasses import replace as _replace  # noqa: E402

from prrs.internal import named_values  # noqa: E402
from prrs.perturbations import apply as apply_probe, propose  # noqa: E402

atoms = build()
cfg_for_propose = cfg or SearchConfig(**CHARGE, closed_shell_only=False)
tracked = [t for t in cfg_json["tracked_coordinates"] if t["kind"] == "bond"]
before = named_values(atoms.positions, tracked)


def _movers(config):
    """Which declared coordinates each proposal can actually move, and its ledger."""
    report = {}
    batches = list(propose(atoms, config, config.seed, report=report))
    moved = {t["name"]: [] for t in tracked}
    for batch in batches:
        probe = batch[0]  # smallest amplitude of this direction
        if probe.mode != "displace":
            continue  # kicks move momenta, not geometry
        candidate = atoms.copy()
        try:
            apply_probe(candidate, probe, config.bond_scale, config.active_atoms)
        except ValueError:
            continue  # a probe this geometry refuses is not reach
        after = named_values(candidate.positions, tracked)
        for name in moved:
            if abs(after[name] - before[name]) > 0.01:  # a tenth of the smallest amplitude
                moved[name].append(f"{probe.family}{list(probe.indices)}")
    return moved, report


basis_movers, basis_report = _movers(_replace(cfg_for_propose, direction_mode="exhaustive"))
run_movers, run_report = _movers(cfg_for_propose)

# Provenance, because this gate and a run do NOT look at the same thing and that has
# already misled once: the gate proposes on the INPUT structure at cfg.seed, while a run
# proposes on the RELAXED source microstate with a direction seed derived per microstate
# (SeedSequence over cfg.seed, the node id, the microstate count and trials_spent). b05's
# probe3 selected three C-H stretches that this gate, on the input, reports as
# candidate_not_selected -- both statements are true of their own structure.
import hashlib  # noqa: E402

structure_hash = hashlib.sha256(
    np.round(atoms.positions, 6).tobytes() + atoms.get_chemical_symbols().__str__().encode()
).hexdigest()[:16]

CHANNELS = {
    "substitution": ("r_C_Cl_leaving", "r_O_C_substitution"),
    "elimination": ("r_C_H_beta", "r_O_H_beta"),
}
status = {}
for channel, names in CHANNELS.items():
    if any(run_movers[n] for n in names):
        status[channel] = "covered"
    elif any(basis_movers[n] for n in names):
        status[channel] = "candidate_not_selected"
    else:
        status[channel] = "basis_uncovered"

lines = [f"{'coordinate':24} {'in basis':>9}  {'this run':>9}   verdict"]
for name in [n for names in CHANNELS.values() for n in names]:
    per_coordinate = (
        "covered"
        if run_movers[name]
        else "candidate_not_selected"
        if basis_movers[name]
        else "basis_uncovered"
    )
    lines.append(
        f"{name:24} {len(basis_movers[name]):9d}  {len(run_movers[name]):9d}   {per_coordinate}"
    )

gate(
    "6  no declared channel is basis_uncovered, and censoring is stated",
    all(v != "basis_uncovered" for v in status.values()),
    "\n".join(lines) + f"""

channel verdicts        {json.dumps(status)}
this run's mode         {run_report.get("mode")}   directions {run_report["selected_total"]}
                        of {len(run_report["candidate_ids"])} candidates
coverage_censored       {run_report.get("coverage_censored")}   by {run_report.get("censored_by")}
full basis (exhaustive) {basis_report["selected_total"]} directions, censored """
    f"""{basis_report.get("coverage_censored")}

**Scope of this gate.** structure = input.extxyz as built (sha256 of positions
{structure_hash}), proposal seed = {cfg_for_propose.seed}. That makes it an INPUT
pre-check and nothing more: a run relaxes the source first and proposes on the relaxed
microstate with a per-microstate direction seed, so its selection differs -- measured, on
b05: 152 candidates here against 179 there, and three C-H stretches selected there that
are candidate_not_selected here. **Coverage used to adjudicate a run must come from that
run's own direction_ledger_per_microstate, not from this table.**

**A PASS here does not license a negative result.** It says the grid can ask its own
question. Whether an answer of "not found" is admissible is a separate test, made per
microstate on the run's own direction_ledger (CRITERIA.md 5). Under
direction_mode='{run_report.get("mode")}' with coverage_censored=
{run_report.get("coverage_censored")}, every negative on this grid is `unanswered`.

A coordinate marked candidate_not_selected is NOT a gap in the proposal basis: the
direction exists and the budget did not take it. That is the fixable kind, and it is why
the two words are kept apart here.""",
)


# ---------------------------------------------------------------- gate 7
# Can the identity layer tell the candidate outcomes apart? Crude geometries are enough:
# chemical_key reads the bond graph, and what is being tested is the annotator, not the
# structures. Building them here is NOT a prediction of what the search will find -- it is
# the check that a report of "one channel only" would be believable.
outcomes = {
    "reactants (2 fragments)": (build(), 2),
    "substitution: tBuOH + Cl- (2)": (substitution_products(), 2),
    "elimination: alkene + H2O + Cl- (3)": (elimination_products(), 3),
}
lines, keys, shapes_ok = [], {}, True
for label, (structure, expected) in outcomes.items():
    key = chemical_key(structure, cfg_for_propose.bond_scale, **CHARGE)
    keys[label] = (
        key["hash"]
        if isinstance(key, dict) and "hash" in key
        else json.dumps(key.get("fragments"))
    )
    seen = (
        len(key["fragments"])
        if isinstance(key["fragments"], (list, tuple))
        else key["fragments"]
    )
    shapes_ok = shapes_ok and seen == expected
    lines.append(
        f"{label:38} fragments {seen} (intended {expected})  "
        f"stereo_unresolved={key.get('stereo_unresolved')!r}"
    )
distinct = len(set(json.dumps(k, sort_keys=True) for k in keys.values()))
gate(
    "7  the candidate outcomes are distinguishable by the identity layer",
    distinct == len(outcomes) and shapes_ok,
    "\n".join(lines)
    + f"""
distinct identities     {distinct} of {len(outcomes)}
The geometries here are crude on purpose and are NOT a claim about what exists on the
surface: this gate asks whether the annotator could tell these outcomes apart if the
search reached them. A mismatch between 'fragments' and 'intended' means the crude
geometry is wrong, not that the chemistry is -- fix the constructor above, not the
criteria.""",
)

# ---------------------------------------------------------------- gate 10
# The stored input must be what build.py builds. Every other gate here measures the
# in-memory build(), while a run searches the FILE, so a stale file makes every gate a
# statement about a structure nobody searched. That is not hypothetical: this grid's
# input.extxyz sat at a retracted draft's 4.30 A separation for two probe runs, because
# build.py raised ModuleNotFoundError at its own rule-2 check and exited 1 BEFORE writing,
# after printing every measurement above it. The printout was read as proof of the write.
# So the artifact gets checked, not the printout.
from ase.io import read as read_structure  # noqa: E402

stored_path = Path(__file__).with_name("input.extxyz")
fresh = build()
if stored_path.exists():
    stored = read_structure(stored_path)
    same_atoms = list(stored.get_chemical_symbols()) == list(fresh.get_chemical_symbols())
    # extxyz is text: 1e-6 A is far below any bond criterion and far above the format's
    # rounding, so a mismatch here is a different structure, not a serialisation artifact.
    deviation = (
        float(np.abs(stored.positions - fresh.positions).max())
        if same_atoms and len(stored) == len(fresh)
        else float("inf")
    )
    charge_ok = stored.info.get("charge") == fresh.info.get("charge")
else:
    stored, same_atoms, deviation, charge_ok = None, False, float("inf"), False
gate(
    "10 the stored input.extxyz is what build.py builds",
    same_atoms and deviation < 1e-6 and charge_ok,
    f"""file                    {stored_path}
present                 {stored_path.exists()}
same formula            {same_atoms}
max coordinate deviation {deviation:.3e} A   (threshold 1e-6)
charge stamped          stored={None if stored is None else stored.info.get("charge")} """
    f"""build()={fresh.info.get("charge")}
r(O...C)                stored="""
    f"""{"-" if stored is None else f"{stored.get_distance(OXYGEN, CENTRAL_C):.4f}"} """
    f"""build()={fresh.get_distance(OXYGEN, CENTRAL_C):.4f}
A FAIL means the runs and the gates are looking at different structures. Regenerate with
`python benchmarks/b05_tertiary_halide/build.py benchmarks/b05_tertiary_halide/input.extxyz`
and CHECK ITS EXIT CODE -- that script prints its measurements before it writes.""",
)

# ---------------------------------------------------------------- gate 9
# Can the identity layer even name this system's reactions? The event key is
# canonicalised under the endpoints' automorphism group, and when that group is not
# enumerated inside `automorphism_limit` the key would depend on the limit rather than on
# the reaction -- so `chemistry` refuses with IdentityUnavailable (implementation
# 5d905aaa, fail-closed by design). This gate exists because b05's first probe attempt
# discovered that by DYING at trial 1, which is the expensive way to learn it.
import time  # noqa: E402

from prrs.chemistry import automorphisms  # noqa: E402

budget = cfg_for_propose.automorphism_limit
rows, ok9 = [], True
for label, structure in (
    ("reactants", build()),
    ("substitution products", substitution_products()),
    ("elimination products", elimination_products()),
):
    started = time.perf_counter()
    permutations, info = automorphisms(structure, cfg_for_propose.bond_scale, None, budget)
    elapsed = time.perf_counter() - started
    enumerated = bool(info.get("enumerated"))
    ok9 = ok9 and enumerated
    rows.append(
        f"{label:24} enumerated={str(enumerated):5}  |Aut|={len(permutations):5d}  "
        f"{elapsed * 1000:8.1f} ms  reason={info.get('reason')}  "
        f"colour_classes={info.get('colour_class_sizes')}  "
        f"nodes={info.get('search_nodes')}"
    )
gate(
    "9  the automorphism group enumerates inside automorphism_limit",
    ok9,
    "\n".join(rows)
    + f"""
automorphism_limit      {budget}   (budget unit: search_nodes)
This substrate has |Aut| = 1296 = 3! (three equivalent methyls) x (3!)^3 (their
hydrogens). The kernel that took the cartesian product over colour-class permutations
could not enumerate it: 3! x 9! = 2,177,280 candidates -- exactly computable, and exactly
the wrong scaling -- needing a limit of 10,000,000 and 19.8 s, so b05's first probe
attempt died at trial 1 with IdentityUnavailable. Replaced 2026-09-09 by point-wise
backtracking with an incremental edge/non-edge check (user ruling): the same 1296
permutations, byte-for-byte, inside the default budget. The products enumerated before
and got faster anyway, which is the same pruning showing up where it was not needed.

A FAIL here is still NOT a knob set too low: it is a scaling statement about the
enumerator, and the answer taken in 2026-09-09 was to fix the pruning rather than to
raise the limit -- raising it would have passed the pre-check provably (2,177,280 is
exact) while still examining 2.17 million mappings per admitted structure. If this gate
ever fails again, that is the shape of the finding, not a number to bump.""",
)

# ---------------------------------------------------------------- gates 4, 5, 8
if "--models" in sys.argv:
    spec = json.loads(Path(__file__).with_name("model.json").read_text())
    try:
        from mace.calculators import MACECalculator
        from prrs import torch_guard

        torch_guard.install(trusted=[spec["path"]])
        calc = MACECalculator(
            model_paths=spec["path"],
            model_type=spec["mace_model_type"],
            device="cuda",
            default_dtype="float64",
        )
        energies = []
        for q in (0, -1, +1):
            reading = build()
            reading.info["charge"], reading.info["spin"] = q, 1
            reading.calc = calc
            calc.results.clear()
            energies.append(float(reading.get_potential_energy()))
        spread = max(energies) - min(energies)
        gate(
            "4  the model responds to total charge ON THIS STRUCTURE",
            spread > 1e-6,
            f"""E(0)={energies[0]:.6f}  E(-1)={energies[1]:.6f}  E(+1)={energies[2]:.6f}
spread = {spread:.4f} eV
b04 measured 16.67 eV on its own structure; this is the same check re-run here, because
charge response is a property of the model AND the system. A zero spread would void the
grid. A nonzero spread says the model is not blind -- it says nothing about accuracy.""",
        )

        # gate 8: where does the -1 sit on a three-fragment product set? PolarMACE returns
        # per-atom charges, so this is readable rather than inferred.
        elimination = elimination_products()
        elimination.info["charge"], elimination.info["spin"] = -1, 1
        elimination.calc = calc
        calc.results.clear()
        elimination.get_potential_energy()
        charges = calc.results.get("charges")
        if charges is None:
            gate(
                "8  the -1 localises on the leaving halide, not smeared",
                False,
                "the calculator returned no per-atom charges",
                required=False,
            )
        else:
            charges = np.asarray(charges, dtype=float)
            on_chloride = float(charges[CHLORINE])
            gate(
                "8  the -1 localises on the leaving halide, not smeared",
                on_chloride < -0.5,
                f"""sum q = {charges.sum():+.4f}   q(Cl) = {on_chloride:+.4f}
q(O) = {charges[OXYGEN]:+.4f}   q(central C) = {charges[CENTRAL_C]:+.4f}
On a dissociated three-fragment set the charge has to sit on one fragment. No model in
this project has been checked on such a state (model.json, unverified_for_this_grid).
This is a WARN, not a veto: a smeared charge would be a real finding about the model,
and it belongs in the record rather than in a pass/fail.""",
                required=False,
            )

        sys.path.insert(0, str(ROOT / "scripts"))
        from model_registry import describe  # noqa: E402

        declared = set(describe(spec["path"])["atomic_numbers"])
        present = sorted(set(int(z) for z in build().numbers))
        missing = [z for z in present if z not in declared]
        gate(
            "5  every element in the structure is one the model declares",
            not missing,
            f"""{spec["name"] if "name" in spec else spec["chosen"]}: declares {len(declared)} elements
structure needs {present} -> {"all present" if not missing else f"MISSING {missing}"}
Read from the checkpoint (scripts/model_registry.py), not from a hand-kept list.
This is the DOMAIN question only; chemistry.VALENCE answers a Lewis question the model
has no opinion on.""",
        )
    except Exception as exc:  # noqa: BLE001 - reported, not raised
        gate("4/5/8  model gates", False, f"{type(exc).__name__}: {exc}")
else:
    for skipped in (
        "4  the model responds to total charge ON THIS STRUCTURE",
        "5  every element in the structure is one the model declares",
        "8  the -1 localises on the leaving halide, not smeared",
    ):
        print(f"  [SKIP] {skipped}")
    print("         pass --models to run them (needs a GPU and loads one checkpoint)\n")

print("=" * 78)
required = [(n, ok) for n, ok, req in results if req]
passed = sum(1 for _, ok in required if ok)
print(f"required gates: {passed}/{len(required)} pass")
for name, ok, req in results:
    if not ok:
        print(f"   {'FAIL' if req else 'WARN'}  {name}")
sys.exit(0 if all(ok for _, ok in required) else 1)
