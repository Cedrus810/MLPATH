# PARS — Pareto-Adaptive Reaction Path Search

Review and compatibility notes: [2026-10-04 fixes (Chinese)](../docs/FIXES_2026-10-04.md); 2026-10-05 electronic-state handoff, Hessian-sample uncertainty gate and verify-before-filter ingestion ([review](../docs/REREVIEW_2026-10-05.md), fix record and real b04 check in [DONE_LOG](../docs/DONE_LOG.md)). The planned maturation workflow is a separate downstream diagnostic step; the current PARS commands are `select` and `refine`.

An **independent module** next to [PRRS](../README.md): PRRS discovers
reactions; PARS organizes and selects reaction paths, then hands
representative reaction tubes to NEB/string refinement and QM validation.

```
PRRS discovers reactions;  PARS organizes and selects reaction paths;
NEB/string refines them;   QM validates them.
```

PARS treats reaction-path search as a multi-objective, multi-fidelity
selection problem over an ensemble of reaction trajectories: candidates are
partitioned into reaction channels and then mechanistically distinct tubes,
pruned within each tube by Pareto geometry under an adaptive margin,
and compressed to a small representative cover before any
expensive refinement is spent. Near-frontier candidates stay recoverable so
that MLP errors cannot irreversibly delete a path.

Design document: `../Pareto-Adaptive Reaction Path Search (PARS) — MLPATH
Integration Design.md` (this implementation covers MVP-1..5 of section 43).

---

## Independence contract

Everything under `pars/` is liftable: copy the folder anywhere and it is a
standalone project (`pyproject.toml`, `src/pars/`, `tests/`, `Makefile`,
`LICENSE`).

The rules that make this true, enforced by `tests/test_boundary.py`:

1. **No `import prrs`, anywhere.** The integration surface with PRRS is a
   *file contract*, not a Python API: PARS reads a run directory
   (`network.json`, `config.json` manifest, `attempts.jsonl`,
   `trials/*/response.traj|observations.jsonl|result.json|endpoint.extxyz`).
   The single reader is `src/pars/adapters/prrs_run.py`.
2. **Authority stays on the PRRS side.** Channel identity, symmetry and bond
   perception are consumed exactly as PRRS materialized them
   (`reaction_channels` / `event_key`, per-frame `edges`); PARS never
   re-derives them. Uncertainty is consumed as the per-frame
   `uncertainty_eV_A` the run recorded. Runs without `reaction_channels`
   (old format) are rejected fail-closed.
3. **Missing quantities are never fabricated** (design section 7): an
   objective the artifacts do not carry is `None`, excluded from
   normalization, and reported as unavailable. Dominance, margin and cover
   compare only candidates with identical available objective sets; different
   masks remain incomparable to avoid cyclic pruning and invented coverage.

   Consequence for uncertainty-aware tiering: no V1 objective carries a
   measured uncertainty (`objective_uncertainty` stays empty unless a caller
   sets it), so in the default configuration `CERTAINLY_DOMINATED` is never
   assigned; dominated paths go to the near-frontier reserve or the cold
   archive. Hard deletion needs `dominance_mode="nominal"`, a positive
   `dominance_tolerance`, or complete objective uncertainty intervals for both
   compared candidates. Partial intervals cannot establish certainty; an
   unrelated candidate's intervals cannot license deletion of a missing pair.
   Positive tolerance alone is an explicit point-estimate exclusion policy
   when neither candidate supplies intervals, not measured statistical certainty.

With `preserve_mechanism_diversity=True` (the default), tubes are established
over **all** candidates before tiering. Dominance and margins are evaluated
within a tube, so a worse cheap objective cannot erase a distinct tube before
the diversity cover runs. Pruned members remain as tube alternates.
`False` explicitly restores channel-wide pruning. Selection records state
`dominance_scope`; `refinement_order` visits one representative per tube first.
Regenerate old selections to use this policy; representative counts can grow.
Contact signatures remain an approximation: identical signatures can hide
geometrically distinct saddles. The independent CPU validation and limits are
documented in [the two-saddle benchmark](../benchmarks/pars_two_saddles/README.md).

Extraction checklist: take the folder, adjust `RUFF ?=` in the `Makefile`
(it points at the parent's `.venv-tools` by default), `pip install -e .`,
done. Dependencies are only `numpy` and `ase`.

---

## Install & run

Tests run without installation (`pythonpath=["src"]`); the interpreter is
the science environment, the same discipline as the parent repository:

```bash
cd pars
make test                      # pytest, needs openmm_dev
make check                     # fmt-check + lint + test
```

Selection — consumes a run directory, writes the Pareto selection:

```bash
pars select --run-dir ../runs/<experiment>/<seed>/<arm> --out ../runs/<experiment>/pars
```

Outputs (all JSON written atomically, NaN rejected):

```
pars_selection.json    config echo, normalization parameters (frozen),
                       per-channel tiers / tubes / representatives / audit,
                       counts incl. skipped attempts by status/classification
                       and channel members that were not ingested
candidates.jsonl       one summary per candidate: objectives, provenance
representatives/       rep-<candidate-id>.extxyz: resampled NEB seed chains
```

Refinement — refines representative seeds through a backend:

```bash
pars refine --selection ../runs/<experiment>/pars \
            --calculator-factory mypkg.calculators:make_mace
```

`--calculator-factory module:attribute` resolves at runtime and must return
an ASE calculator; PARS has no static dependency on any physics engine (a
thin wrapper around PRRS's `mace_factory` works). Refinement is staged
(endpoint relaxation and curvature checks → NEB → climbing image), monitors
pathological bands, and reports `UNRESOLVED_MLP` when an explicit uncertainty
threshold is exceeded or measurements are missing on any image. NaN/Inf or
negative uncertainty is `FAILED` even without a threshold. Endpoints that do
not relax, have negative curvature, mismatch independently relaxed basin
witnesses, or relax into the same minimum are `FAILED`. New selections carry
mapped source/quench basin witnesses; old seeds without them explicitly record
`basin_identity: unavailable` and do not establish endpoint identity.
A plain-NEB band whose top is at an endpoint or unresolved reports no barrier.
A climbing image must have exactly one negative projected Hessian eigenvalue
beyond `refinement_saddle_curvature_tolerance` (default 1e-4 eV/A^2), using
central force differences at `refinement_saddle_hessian_step` (default 0.001 A).
An unresolved or higher-order stationary point is `FAILED` without a barrier.
This adds 6N force evaluations and records the spectrum and Cartesian negative
mode in `saddle_check`. Each signed negative-mode kick must lower the saddle
energy and independently descend to a checked minimum, with strict energy
decrease and a unique pairing to the two refined endpoints. Missing connectivity,
same-basin descents, an intermediate basin, ambiguous matching or exhausted
budgets refuse CI success without a barrier. `connectivity_check` records the
signed starts, energies, landed positions, force residuals, Hessians and mapped
RMSDs. With source/quench witnesses its identity scope is `source_basin_witnesses`;
legacy seeds without them establish only `refined_endpoints_only` connectivity.
Explicit uncertainty thresholds also apply to descent trials; missing/excessive
measurements are `UNRESOLVED_MLP`, invalid values `FAILED`.

Connectivity also requires a mandatory finer replay: kick, gradient step and
displacement cap are multiplied by `refinement_descent_resolution_factor`
(default 0.25), with the same force tolerance, per-sign step limit, calculator
and hard budget. Coarse and fine endpoint labels must agree. A finite kick can
otherwise jump a narrow intermediate well while all recorded steps lower
energy. Sensitive/incomplete/uncertain finer replay refuses the barrier;
`resolution_check` retains its diagnostics. Success means
`stable_on_tested_scales`, with `continuous_connection_certified=false`.
This finite check cannot prove continuous gradient flow or exclude still
narrower hidden basins. Old single-resolution records need a new refinement.

Descent defaults are a 0.05 A Cartesian kick, gradient step 0.05 A^2/eV,
maximum trial displacement 0.05 A per atom, fmax 0.001 eV/A (capped by the
endpoint tolerance), 2000 accepted steps per sign and 20 backtracking halvings.
All `refinement_descent_*` settings are frozen in refinement conditions. Descent
and 6N force evaluations per landing use the same calculator and hard budget.
This is a numerical Cartesian witness on the supplied MLP, not mass-weighted
IRC or QM validation. Plain NEB retains a discrete band estimate with explicit
`connectivity_check.status=not_checked`. Existing refined records retain their
old conditions; rerun refinement to use this implementation and its extra cost.
NEB uses projected springs (improved
tangent) with `refinement_spring_k` = 5 eV/A^2. By default rigid motion is
removed from seeds and during NEB (`refinement_remove_rigid_motion`), which is
right for an isolated molecule and WRONG for an embedded core or any external
field: set it to false there.

Results go to `<selection>/refined/`: `refined-<id>.extxyz` and
`refined-<id>.json` per candidate, written as each finishes (the JSON is the
commit marker), and `refined_paths.json` aggregating every record in the
directory. Each record carries its `conditions`: seed sha256, actual calculator
configuration/model weights, declared model file hashes, factory/implementation
fingerprints, numerical package versions, every `refinement_*` setting and PARS
version. Factories can export `calculator_metadata()` returning JSON metadata;
alternatively pass `--calculator-metadata metadata.json`. Name external model
files and settings there; content hashes are recomputed every invocation.
Opaque calculator state is refused instead of claiming a verified identity.

Re-running `pars refine` constructs the calculator to verify its identity, then
skips candidates already refined under the same conditions. Selected candidates
with different/unknown conditions refuse before force evaluations; `--redo`
explicitly replaces them. Legacy records need a redo; re-select to add basin
witnesses to old seeds. Old source runs are read-only.

`--only --redo` may create different refinement cohorts. `refined_groups.json`
lists each cohort's conditions and candidate ids, with `homogeneous: false`
for mixed conditions; the CLI prints a notice. The flat `refined_paths.json`
keeps compatibility and records `conditions_sha256` membership. Compare results
within groups. A candidate's seed hash is excluded from its cohort key.

Endpoint curvature uses central force differences (`refinement_endpoint_hessian_step`
= 0.001 A), with rigid modes projected out only when rigid-motion removal is enabled.
The negative-curvature floor is `refinement_endpoint_curvature_tolerance`
= 1e-4 eV/A²; mapped basin witness RMSD tolerance is
`refinement_endpoint_identity_rmsd` = 0.05 A. These are explicit numerical
checks on the supplied PES and do not certify a TS or replace independent QM.
Each checked geometry adds 6N force evaluations.

---

## Module map (design sections in parentheses)

| path | responsibility |
|---|---|
| `schema.py` | `DiscoveryTrace → PathCandidate → ReactionTube → RefinedPath`, lifecycle enums (6, 27) |
| `config.py` | one frozen `ParsConfig`, JSON only, unknown keys rejected, rationale per default (30, 36) |
| `adapters/prrs_run.py` | the whole PRRS integration surface: run directory reader |
| `paths/ingest.py` | run directory → traces, classification filter (6.1) |
| `paths/candidate.py` | event-window trimming; forcing stage dropped, kept in provenance (6.2, 17.1) |
| `paths/resample.py` | arc-length reparameterization to NEB images (17.4) |
| `paths/fingerprint.py` | bond-signature sequences, discrete Fréchet path distance, recrossing (15, 8.5) |
| `paths/cluster.py` | reaction tubes by leader clustering (16) |
| `frontier/objectives.py` | cheap objectives: excursion, trigger, uncertainty, curvature, recrossing (8) |
| `frontier/normalize.py` | robust median/IQR normalization, parameters frozen into output (10) |
| `frontier/dominance.py` | nominal dominance with indistinguishability tolerance; certain dominance needs interval separation (11, 12) |
| `frontier/margin.py` | Γ(τ) directed margin, four tiers, γ rescue (13) |
| `frontier/cover.py` | greedy deterministic cover, tube-aware representative selection (14–16) |
| `frontier/archive.py` | per-channel Pareto archive, audit sample for the later QM experiment (4, 20, 28) |
| `refinement/base.py` | backend-independent refinement contract (43 MVP-5) |
| `refinement/ase_neb.py` | ASE NEB adapter: staged NEB/CI-NEB, pathology guard, uncertainty gate (18) |

CLI: `pars select`, `pars refine` (`python -m pars` works too).

---

## What V1 deliberately does not do (design section 44)

No universal Pareto library, no dynamic NEB image redistribution (V2,
section 24), no MLP retraining, no Markov bridge, no multistate/MECP handling
(`electronic_state_status` exists as a field; Stateweaver handoff is an
adapter for later), no learned path encoder. The one hypothesis V1 tests:
*can Pareto geometry + mechanism-aware covering reduce refinement cost
without destroying reaction diversity?*

Next steps on the design's MVP ladder: MVP-6 — a fixed QM-check budget
measuring ranking inversion, frontier stability, margin rescue and QM
efficiency (the audit sample in every selection record is its feeding
channel; `withheld_candidates` in `network.json` is the precedent).
