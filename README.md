# PRRS — Perturbation–Response Reaction Search

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
![Lint: ruff](https://img.shields.io/badge/lint-ruff-261230)

English | [简体中文](README.zh-CN.md)

**MLPATH project entry:** [architecture and current status (Chinese)](docs/PROJECT_OVERVIEW.md) · [documentation index](docs/README.md) · [PARS path module](pars/README.md) · [path maturation layer](mature/README.md).
MLPATH includes PRRS discovery, PARS path selection/refinement (MVP-1..5 landed), and the shooting-based maturation layer (`mature/`, M0+M1 landed). This README covers PRRS; the project overview tracks their boundaries and connect-mode progress.

**PRRS is a reaction-discovery method designed for machine-learning interatomic
potentials (MLPs).**

> Perturb a molecule in a controlled way, remove the perturbation, let the system
> respond on the MLP potential-energy surface, quench, and observe which chemical
> state it reaches.

PRRS does **not** require a predefined reaction coordinate, a predefined product,
or an initial guess of the transition state. It discovers possible reaction
channels and products first; transition states are located and verified
afterwards.

---

## How it works

```
seed state → controlled perturbation → short MLP trajectory → quench → classify
                                                                       │
     conformer of source ◄── same chemical identity ──► new chemical state
                                                                       │
                              first-order saddle (n₋ = 1) ──► TS pool
```

Each trial is a small numerical experiment in two phases: a **delivery phase**,
in which the perturbation is applied incrementally to internal coordinates, and
a **free response phase**, in which the perturbation is removed and the system
relaxes without bias. Three possible outcomes:

1. **Return to the source minimum** — elastic evidence. This direction does not
   easily lead to another basin at the tested amplitude.
2. **Relaxation to another minimum** — a basin is discovered. Conformer and new
   chemical state are distinguished explicitly: what decides is chemical
   connectivity, not energy or distance.
3. **A stationary point with exactly one negative Hessian eigenvalue** — a
   transition-state candidate, admitted to a TS pool and verified by following
   the downhill directions on both sides to see which states it connects.

### Design choices that carry the method

- **Perturbations act on primitive internal coordinates** — bond lengths, bond
  angles, torsions, collective fragment motions, intermolecular approach — not
  Cartesian noise. The reason is measured, not aesthetic: early ethanol trials
  with radial-only perturbations returned to the start in all cases while three
  real conformational basins sat 4 meV away, because pairwise radial moves never
  couple into torsions within the trajectory length. **The perturbation set must
  span the chemically relevant degrees of freedom.**
- **Delivery is incremental** (≤ 0.05 Å/rad per step, gradient re-evaluated every
  step), because a large perturbation is a sequence of locally valid small ones.
  Requested target displacement and collateral damage are reported separately;
  covalent-bond damage and atom overlaps are hard gates.
- **Symmetry is handled explicitly.** Chemically equivalent atoms never generate
  artificial "new reactions": identity is hashed by automorphism-invariant graph
  invariants, conformers are deduplicated by symmetry-aware RMSD, and stereo
  signatures use configurational chirality only. Symmetric rotors are *demoted,
  not removed* — one of them contributed a real methyl-rotation saddle point.
- **Three identity layers that are never collapsed**: the *chemical state*
  (connectivity + configurational stereochemistry), the *conformer reservoir*
  (many geometries per state, because reaction accessibility depends on the
  starting conformer), and the *reaction event* (bond-difference keys
  canonicalised under the automorphism groups of both endpoints). A degenerate
  reaction — such as an equivalent-H proton transfer — is a self-loop, not a new
  node.
- **Convergence criteria are mode-resolved.** A fixed force tolerance controls
  how small the forces are, not how far the structure is from the stationary
  point; ethanol torsions are ~150× softer than bond stretches, so soft modes
  get their own gradient-over-stiffness criterion while genuinely free modes are
  recorded, not forced to "converge".
- **Fail closed.** Anything unresolvable is rejected with a stated reason, never
  replaced by a default; existence claims must be proven by execution, not by
  records.

---

## What PRRS does — and does not — claim

| Question | Comparison | What it tests |
|---|---|---|
| Is the **search** good? | PRRS vs a reference search on the *same* MLP surface | the method |
| Is the **potential** accurate? | MLP transition state vs QM | the energy model |

PRRS is responsible for discovering channels, identifying products, locating and
verifying transition-state candidates, and building the network. It is **not**
expected to deliver final quantitative barriers: barrier heights and imaginary
frequencies are properties of the chosen surface, checked downstream by QM.
Mixing the two comparisons makes errors unattributable — in the benzoylacetone
three-model comparison every chemical-identity field agreed across models while
barrier heights differed 3.7×, which is exactly the separation this design
enforces.

---

## Status

The benchmark suite (`benchmarks/`) follows one rule: one cell = one reaction
system + acceptance criteria frozen *before* the run + a script that grades the
result mechanically against those criteria.

| Stage | System | New axis tested | Status |
|---|---|---|---|
| P0 | Ethanol | Conformer discovery; framework | Completed — 13/13 |
| P1 | Malonaldehyde | Reaction event — **degenerate** self-loop | Completed — 13/13, QM-checked |
| P2 | 3-Oxobutanal enol | Reaction-network expansion A→B | Completed — 16/16 seeds |
| P3 | Benzoylacetone | Aromaticity, chelation, competing isomers | Completed — criteria frozen and passed |
| b04 | Cl⁻ + CH₃Cl | **Total charge**, 2 fragments, heavy-atom coordinate | Completed — criteria frozen 2026-09-13, TS QM-confirmed |
| b05 | Tertiary halide + nucleophile | Substitution vs elimination competition | In progress — D2-HF v1.4 two-seed full runs done (2026-11-08); HS bracket rule frozen at v7 |
| b06 | Electrophilic addition substrate | Polar addition | Planned |
| b07 | Cyclohexene (retro-DA) · butadiene + ethylene (DA) | **Fragment number 1↔2** | Planned — n1 Diels–Alder probe scripts landed |

P0–P3 are now frozen regression fixtures: their results remain valid, and formal
acceptance passes to the P4+ independent-QM protocol. The series deliberately
moves away from repeated proton-transfer examples toward charged systems,
intermolecular reactions, competing pathways, heavy-atom bond breaking/formation,
and changes in fragment number.

### Independent QM cross-checks

Protocol: **r2SCAN-3c** for geometries and frequencies → **ωB97M-V/def2-TZVPD**
single points on the same geometry. (r2SCAN-3c alone misestimates the SN2 central
barrier by ~11 kcal — geometries are good, energies must be re-checked.)

- **Malonaldehyde (P1).** Bidirectional IRC from the PRRS saddle converges to the
  expected minima, both endpoints with zero imaginary frequencies. The MLP and QM
  saddle geometries agree to **0.0151 Å** (symmetry-aware RMSD). The MLP barrier
  (419 meV) vs QM (134 meV) differ substantially — an energy-model error, not a
  failure to locate the reaction geometry.
- **SN2: Cl⁻ + CH₃Cl (b04).** The transition state found by PRRS refined to
  2.3645/2.3645 Å, matching an independently located QM saddle to four decimal
  places. Bidirectional IRC reached the expected reactants and products, and the
  reaction-event key matches the PRRS one exactly.

> **PRRS identifies the reaction event and transition-state geometry even when
> the MLP does not reproduce QM barrier energies accurately.**

### Choice of machine-learning potential

P0–P3 are frozen on **MACE-OFF24** and are not retrospectively rerun or relabelled.
From b04 onward the default is **MACE-POLAR-1-M**, chosen on QM cross-comparisons
rather than preference: it is the only candidate tested that holds up on both
neutral proton-transfer chemistry and the ionic SN2 benchmark. MACE-OFF24 does
not respond to total molecular charge at all (bitwise-identical energies for
charge 0/−1/+1), which disqualifies it for ionic chemistry. Frozen results never
automatically serve as references for later stages.

---

## Installation

Requirements: **Python ≥ 3.10**; the core dependencies are ASE and NumPy only.

```bash
pip install -e .            # core: prrs CLI; runs on any ASE calculator
pip install -e ".[mace]"    # MACE support (mace-torch)
pip install -e ".[test]"    # pytest
```

The MACE-OFF24 medium checkpoint used by the frozen P0–P3 fixtures ships in
`models/` (sha256 prefix `e5ccf5837f685899`). Each benchmark declares its model
in `benchmarks/b0*/model.json` — **no model is ever downloaded implicitly**.

### Two environments, on purpose

The **science environment** (conda spec in `env/openmm_dev.yml`; baseline and
history in `env/BASELINE.md`) produces every published number, and its package
versions are part of each result. Developer tooling (ruff) lives in a separate
`.venv-tools` (`make tools`) so that formatting can never move a scientific
dependency. Backend parity once went red on an openmm *patch* release — never
install dev tools into the science environment.

## Quickstart

Optional local minimum / TS-candidate refinement is available as `prrs optimize`.
It supports ASE calculators and existing OpenMM/MACE contexts through geomeTRIC;
see [commands, outputs and validation boundaries](docs/GEOMETRIC_REFINEMENT.md).

```bash
# minimal self-contained demo
prrs demo --output runs/demo --seed 42

# a search run
prrs search input.extxyz --output runs/my_run \
    --config config.json --calculator mace --calculator-kwargs calc.json

# connect mode: given a product, every microstate records its distance to it
# (recorded, never routed)
prrs search input.extxyz --product product.extxyz --output runs/my_connect \
    --config config.json --calculator mace

# export an unrefined path seed for downstream NEB / dimer / string refinement;
# quench iterations are excluded — they are not dynamics
prrs path --run runs/my_run --attempt <id> --output seed.extxyz --images 9
```

Every run writes a manifest (including `implementation_sha256`, see
[Development](#development)) and a `network.json` reaction network.
`runs/README.md` records which on-disk runs count as evidence — a `network.json`
merely *existing* does not mean the run finished.

### Watching a live run

```bash
python scripts/prrs_watch.py runs/my_run
```

A read-only terminal UI: it tails each trial's `observations.jsonl`
incrementally, polls `network.json`, and draws the evolving potential-energy
landscape directly in the terminal (true colour). It never writes to `runs/`;
full-resolution PNGs go to `.scratch/prrs_watch/<run>/`.

---

## Repository layout

```
src/prrs/                 search core — 25 modules, ≈13k lines
  config.py               all protocol parameters + validation; invalid configs rejected before launch
  internal.py             primitive internals, exact rigid-fragment rotation, continuation displacements
  chemistry.py            the three identity layers, automorphism groups, symmetry-aware RMSD
  canonical.py            certified canonical labelling by individualization–refinement
  perturbations.py        six perturbation families, symmetry quotient, layered direction budget
  network.py              registry: admission, conformer reservoirs, domain gates
  state.py                active-subset graph encoding, Kabsch RMSD, response classification
  reliability.py          per-force-call fail-closed gating
  runner.py               one trial: delivery, mode-resolved quench, Hessian/min checks, saddle descent
  search.py               budget scheduling, amplitude refinement, TS pool, run manifest
  completeness.py         conditional completeness certificates for basin membership
  product.py              connect-mode product target: records distances only, never routes
  window_audit.py         matched response-horizon replay audit (never certifies a static boundary)
  geometric_backend.py    optional geomeTRIC refinement on ASE or an existing OpenMM Context
  geometric_cli.py        CLI wiring for the refinement (lazy import)
  provenance.py           calculator identity for resume verification (no execution)
  schema.py               explicit network.json schema (schema 2)
  calculators.py          calculator factory
  openmm_backend.py       tested openmm-ml/MACE primitives
  torch_guard.py          single audited wrapper around torch.load
  validation.py           independent cross-check path (ase.vibrations) — never imported by the search
tests/                    36 files, 388 test functions
benchmarks/               the benchmark suite — one directory per reaction system (criteria + grader)
docs/                     working records: handoff, acceptance docs, plans, experiment logs
scripts/                  provenance, implementation hash, live run observer
runs/                     experiment outputs — runs/README.md says which count as evidence
models/                   MACE-OFF24_medium.model (the frozen P0–P3 fixture model)
env/                      science-environment baseline: BASELINE.md, conda spec, pip freeze
```

---

## Development

```bash
make help      # what each target does
make check     # fmt-check + lint + test + digest, in that order — the full gate
make digest    # golden-sample regression gate
make hash      # the implementation_sha256 a run manifest would record
```

- **pytest asks whether the code is correct; digest asks a different question** —
  whether the code still produces the same numbers as the runs already on disk
  (332 files / 2136 frames, with every difference classified as
  `graph_invariant` / `parity_only` / `key_only`). **Zero diff is the only
  accepted evidence of a zero-behaviour change.** Intentional changes are passed
  with `--allow parity_only` etc. — but an allowance confirms the *shape* of a
  difference, not that the new values are correct; that needs chemical judgement.
- **`implementation_sha256`** (`make hash`) is the provenance anchor of every run
  manifest: a content hash of `src/prrs/*.py`, independent of VCS state. Never
  copied by hand; the fingerprint of record is tracked in `docs/HANDOFF.md` and
  `env/BASELINE.md`.
- Last recorded full gate: **537 passed, 7 skipped, digest zero-diff**
  (2026-10-01; the skips are environmental — no MACE checkpoint / CUDA on the
  gate machine).

### Project norms

These are working rules, not suggestions; results produced in violation of them
have been voided before.

- **Acceptance criteria are frozen before the data exists** and never edited in
  place afterwards; a change is an explicit revision recording the old text, the
  reason, and which data had been seen when the revision was made.
- **Interventions are named separately**, so an improvement can be attributed to
  one change at a time.
- **Incomplete ≠ failed.** The denominator is the number of *initiated* trials;
  incompletes and crashes are recorded as such, with causes
  (`scripts/scan_status.sh` distinguishes the three states mechanically).
- **Fail closed.** When something cannot be determined it is rejected with a
  stated reason — never replaced by a default value.
- **Gates are not moved because they block a wanted result.** First quantify
  what the gate blocks.
- **Measure before inferring.** Recorded refutations are assets, not
  embarrassments (`docs/HANDOFF.md` §6).

---

## Documentation

| | |
|---|---|
| [docs/README.md](docs/README.md) | document index — what is current, what is only a record |
| [docs/HANDOFF.md](docs/HANDOFF.md) | handoff overview; frozen conclusions live here |
| [docs/DECISIONS_PENDING.md](docs/DECISIONS_PENDING.md) | open items only |
| [docs/PRRS_STATUS.md](docs/PRRS_STATUS.md) | implementation status, measured claim by claim |
| [docs/onboarding/README.md](docs/onboarding/README.md) | newcomer guide |
| [benchmarks/README.md](benchmarks/README.md) | benchmark suite rules and queue |
| [env/BASELINE.md](env/BASELINE.md) | environment baseline and fingerprint history |
| [docs/HF5_SUMMARY_2026-11-10.md](docs/HF5_SUMMARY_2026-11-10.md) | current line — HF5/HF6 summary and plans |
| [docs/PLAN_2026-09-28.md](docs/PLAN_2026-09-28.md) | connect-mode two-layer bridge — G1–G4 landed, stop-at-product termination pending |

Active work as of 2026-11-11: connect mode, PARS (MVP-1..5) and the maturation
layer M0+M1 stand as before. The **D2-HF line** has run v1.1 → v1.4: same-graph
ridges sampled, HS / HC run in parallel waves, the HS bracket rule frozen at v7,
and v1.4 executed as two-seed full runs with lost-mode diagnostics; the
fp32-Hessian diagnosis and the F3 false-ridge audit are closed, and rescued
midpoint candidates are the new E2 TS. Next construction blocks: HF5
cross-system (`docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md`, summary in
`docs/HF5_SUMMARY_2026-11-10.md`) and the HF6 dimer
(`docs/PLAN_HF6_DIMER_2026-11-10.md`); an n1 Diels–Alder probe has scripts in
`benchmarks/n1_diels_alder/`.

The original method note is
[`mlp_native_perturbation_response_reaction_search.md`](mlp_native_perturbation_response_reaction_search.md).

---

## License

This project is licensed under the **GNU Affero General Public License v3.0 only**
(`SPDX: AGPL-3.0-only`). The full text is in [LICENSE](LICENSE).

In short: you may use, study, share, and modify this code, but any derivative
work — whether distributed or offered over a network — must be released under
the same license, with complete corresponding source. "Only" means the code is
licensed under AGPL-3.0 exactly, and not automatically under future GPL/AGPL
versions.
