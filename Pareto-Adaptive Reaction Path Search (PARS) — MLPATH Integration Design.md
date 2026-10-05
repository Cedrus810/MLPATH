# Pareto-Adaptive Reaction Path Search (PARS)

## Design and Integration into MLPATH / PRRS

**Status:** design proposal  
**Parent project:** MLPATH / PRRS  
**Role:** reaction-path selection, compression, refinement, and multi-fidelity resource allocation  
**Core principle:**

\[
\boxed{
\text{PRRS discovers reactions;}
\quad
\text{PARS organizes and selects reaction paths;}
\quad
\text{NEB/string refines them;}
\quad
\text{QM validates them.}
}
\]

---

# 1. Motivation

PRRS currently addresses the reaction-discovery problem:

\[
A
\rightarrow
\{\text{perturbations}\}
\rightarrow
\{\text{responses}\}
\rightarrow
\{B_1,B_2,\ldots\}.
\]

This is intentionally different from conventional NEB.

NEB usually assumes that both endpoints are already known:

\[
A,\ B
\quad\Longrightarrow\quad
\text{search for a minimum-energy path}.
\]

PRRS instead asks:

\[
A
\quad\Longrightarrow\quad
\boxed{\text{What reactions are accessible at all?}}
\]

Once PRRS becomes capable of producing many successful reactive trajectories, however, a second problem appears.

Suppose one search produces

\[
10^3-10^5
\]

candidate trajectories.

Many may correspond to:

- the same product,
- the same mechanism,
- different conformers,
- different perturbation amplitudes,
- noisy variants of the same transition tube,
- high-energy but geometrically useful paths,
- low-energy paths with poor MLP reliability,
- smooth paths with high uncertainty,
- rough paths that nevertheless contain a useful saddle geometry.

The next question is therefore not merely:

> Which trajectory has the lowest MLP barrier?

Instead:

\[
\boxed{
\text{Which small subset of candidate paths preserves the chemically meaningful}
\\
\text{reaction-space diversity while being worth expensive refinement?}
}
\]

This is the role of **Pareto-Adaptive Reaction Path Search (PARS)**.

---

# 2. Position of PARS inside MLPATH

PARS is **not** a replacement for PRRS.

PARS starts after PRRS has produced candidate reactive responses.

The complete architecture becomes:

```text
Stable state / conformer
        │
        ▼
┌──────────────────────────┐
│          PRRS            │
│ perturb → respond        │
│ detect topology change   │
│ discover new basins      │
└────────────┬─────────────┘
             │
             │ raw reactive trajectories
             ▼
┌──────────────────────────┐
│          PARS            │
│                          │
│ path construction        │
│ objective extraction     │
│ Pareto archive           │
│ adaptive margin          │
│ path/mechanism cover     │
│ refinement selection     │
└────────────┬─────────────┘
             │
             │ representative reaction tubes
             ▼
┌──────────────────────────┐
│    Path refinement       │
│ NEB / CI-NEB / string    │
│ saddle refinement        │
└────────────┬─────────────┘
             │
             ▼
┌──────────────────────────┐
│      QM validation       │
│ QM SP / gradients        │
│ TS optimization          │
│ Hessian                  │
│ IRC                      │
└────────────┬─────────────┘
             │
             ▼
      Reaction network
```

The conceptual separation is:

\[
\boxed{
\text{Discovery}
\neq
\text{Path selection}
\neq
\text{Path optimization}
\neq
\text{Electronic-structure validation}.
}
\]

This separation is essential because errors in these four stages have different meanings.

---

# 3. Scientific question

PARS addresses the following problem.

Given a set of candidate reactive paths

\[
\mathcal T=
\{\tau_1,\tau_2,\ldots,\tau_M\},
\]

and a limited expensive-computation budget \(B\), choose a much smaller set

\[
\mathcal C
\subset
\mathcal T,
\qquad
|\mathcal C|\ll M,
\]

such that \(\mathcal C\) preserves:

1. chemically distinct reaction channels,
2. mechanistically distinct paths between the same endpoints,
3. promising low-barrier paths,
4. paths that are geometrically well behaved,
5. paths whose MLP predictions are sufficiently reliable,
6. paths that probe uncertain or potentially incorrect MLP regions,
7. candidates useful for discovering ranking inversions between MLP and QM.

The optimization problem is inherently multi-objective.

Therefore PARS should **not begin by collapsing everything into one scalar score**.

Instead, it constructs a path-level objective vector

\[
R(\tau)
=
\left(
R_1(\tau),
R_2(\tau),
\dots,
R_K(\tau)
\right)
\]

and operates on its Pareto geometry.

---

# 4. Important restriction: Pareto comparison is channel-conditioned

A global Pareto frontier over all reactions is chemically dangerous.

Suppose PRRS discovers:

\[
A\rightarrow B
\]

and

\[
A\rightarrow C.
\]

If the \(A\rightarrow B\) reaction has a lower barrier and lower uncertainty than every \(A\rightarrow C\) path, a naïve global Pareto algorithm could delete the entire \(A\rightarrow C\) mechanism.

That would destroy the purpose of reaction discovery.

Therefore candidates must first be partitioned into reaction channels:

\[
\mathcal T
=
\bigcup_c \mathcal T_c.
\]

A channel is defined primarily by:

- initial chemical state,
- final chemical state,
- symmetry-equivalent connectivity,
- stereochemical identity where relevant.

Pareto pruning is then performed primarily within each channel:

\[
\operatorname{Pareto}(\mathcal T_c).
\]

Thus:

\[
\boxed{
\text{A low-barrier reaction must not erase a distinct higher-barrier reaction.}
}
\]

Reaction-network diversity takes priority over global energetic dominance.

---

# 5. Same endpoints do not necessarily mean same mechanism

Even after grouping by endpoints, multiple mechanisms may exist:

\[
A\rightarrow B
\]

via

\[
M_1,\ M_2,\ M_3.
\]

For example, two paths may have identical reactant and product connectivity but differ in:

- proton-transfer ordering,
- bond-breaking/bond-forming order,
- solvent participation,
- conformational gating,
- intermediate formation,
- approach geometry,
- electronic state.

Therefore PARS introduces a second level:

```text
Reaction channel
    │
    ├── mechanism / tube 1
    ├── mechanism / tube 2
    └── mechanism / tube 3
```

A **reaction tube** is an ensemble of paths representing the same broad mechanism.

The intended hierarchy is therefore:

\[
\boxed{
\text{chemical state}
\rightarrow
\text{reaction channel}
\rightarrow
\text{reaction tube}
\rightarrow
\text{individual path}
}
\]

This prevents objective-space compression from accidentally merging physically distinct mechanisms.

---

# 6. Candidate path lifecycle

PARS should distinguish several objects rather than treating every trajectory as equivalent.

## 6.1 `DiscoveryTrace`

Raw trajectory produced by PRRS.

It may contain:

- forced perturbation frames,
- free-response dynamics,
- basin entry,
- thermal vibration,
- post-reaction relaxation.

It is **not yet an NEB path**.

---

## 6.2 `PathCandidate`

A cleaned reaction segment extracted from a `DiscoveryTrace`.

Conceptually:

\[
A
\rightarrow
\text{transition region}
\rightarrow
B.
\]

The artificial forcing stage should be recorded separately.

The candidate must preserve complete provenance:

```text
source conformer
perturbation family
perturbation parameters
random seed
trajectory ID
start state
end state
event fingerprint
MLP model
model version
```

---

## 6.3 `ReactionTube`

Cluster of mechanistically equivalent `PathCandidate`s.

A tube may contain many realizations:

\[
\mathcal U_j=
\{\tau_{j1},\tau_{j2},\ldots\}.
\]

---

## 6.4 `RefinedPath`

A path after:

- endpoint minimization,
- reparameterization,
- NEB/string optimization,
- optional climbing-image refinement.

---

## 6.5 `ValidatedPath`

A path for which sufficient electronic-structure validation exists.

Possible validation levels include:

```text
MLP-only
QM single-point checked
QM gradient checked
QM TS optimized
Hessian verified
IRC verified
```

These stages must never be conflated.

---

# 7. Path objectives

All objectives should internally be transformed so that

\[
\boxed{\text{smaller is better}}
\]

before Pareto operations.

However, different objective vectors should be used at different stages.

A raw PRRS path does not yet possess a meaningful NEB perpendicular-force residual, so missing quantities should **not** be filled with arbitrary numbers.

---

# 8. Stage I: cheap pre-refinement objectives

For a raw path candidate \(\tau\), possible objectives include:

\[
R_{\rm pre}(\tau)=
\left(
E_{\rm proxy},
A_{\rm trigger},
U_{\rm MLP},
C_{\rm path},
R_{\rm cross},
P_{\rm endpoint}
\right).
\]

## 8.1 Energy excursion

\[
E_{\rm proxy}
=
\max_t E(t)-E_A.
\]

This is **not interpreted as a validated barrier**.

It is merely a cheap path descriptor.

---

## 8.2 Trigger severity

Examples:

\[
A_{\rm trigger}
=
\text{required perturbation amplitude},
\]

or injected energy.

Smaller values may indicate easier access to the channel under the specified perturbation family.

Again,

\[
A_{\rm trigger}\neq\Delta G^\ddagger.
\]

---

## 8.3 MLP uncertainty / reliability penalty

\[
U_{\rm MLP}.
\]

Depending on available models this may include:

- ensemble disagreement,
- latent-space distance,
- committee variance,
- extrapolation score,
- force disagreement,
- known-domain checks.

The existing `reliability.py` machinery should provide these values rather than PARS implementing another uncertainty system.

---

## 8.4 Path roughness / curvature

For a path parameter \(s\),

\[
C_{\rm path}
\sim
\int
\left|
\frac{d^2R}{ds^2}
\right|^2 ds.
\]

In practice this should be computed in an alignment- and symmetry-aware representation rather than naïve Cartesian coordinates.

---

## 8.5 Recrossing penalty

A trajectory that repeatedly crosses back and forth between state definitions may be less suitable as a clean initial path.

Define, for example,

\[
R_{\rm cross}
=
N_{\rm state\ crossings}.
\]

This is a path-quality descriptor, not proof that the reaction is dynamically irrelevant.

---

## 8.6 Endpoint-confidence penalty

Paths ending in an uncertain basin assignment should be penalized:

\[
P_{\rm endpoint}.
\]

The state/canonicalization logic already present in PRRS should remain the authority for this quantity.

---

# 9. Stage II: post-refinement objectives

After NEB/string optimization, a more physically meaningful vector becomes available:

\[
R_{\rm ref}(\tau)=
\left(
\Delta E^\ddagger_{\rm MLP},
F_{\perp}^{\max},
F_{\rm NEB}^{\max},
C_{\rm path},
U_{\rm path},
S_{\rm state}
\right).
\]

Possible components are:

### MLP barrier

\[
\Delta E^\ddagger_{\rm MLP}.
\]

This is meaningful **on the MLP PES only**.

### Maximum transverse force

\[
F_{\perp}^{\max}.
\]

### NEB convergence residual

\[
F_{\rm NEB}^{\max}.
\]

### Geometric curvature / unnecessary path length

Used to identify pathological chains.

### Path-wide MLP reliability

For example,

\[
U_{\rm path}
=
\max_i U(R_i)
\]

or a high quantile rather than the mean.

### Electronic-state warning

For systems where electronic-state information exists,

\[
S_{\rm state}
\]

can encode proximity to a suspected spin/state crossing.

For ordinary universal MLPs this quantity may be unavailable and must not be fabricated.

---

# 10. Objective normalization

Raw quantities have incompatible units:

\[
{\rm eV},
\qquad
{\rm eV/\AA},
\qquad
\AA,
\qquad
\text{dimensionless uncertainty}.
\]

Therefore the frontier should operate on normalized quantities:

\[
\widetilde R_k
=
\frac{R_k-a_k}{s_k}.
\]

Robust choices for \(a_k,s_k\) include:

- median + IQR,
- predefined physically meaningful scales,
- robust quantile ranges.

The normalization parameters must be **saved in the run record**.

Changing them silently would change the frontier.

For prospective benchmarks, scaling rules should be frozen before evaluation.

---

# 11. Pareto dominance

For two candidates \(a,b\),

\[
a\prec b
\]

if

\[
R_k(a)\le R_k(b)
\quad
\forall k
\]

and strict improvement occurs for at least one objective.

Then \(b\) is dominated.

The non-dominated set is

\[
\mathcal F=
\operatorname{Pareto}(\mathcal T_c).
\]

This becomes the first compression step:

\[
M
\rightarrow
|\mathcal F|.
\]

However, strict deterministic pruning is unsafe when objectives themselves are noisy.

---

# 12. Uncertainty-aware dominance

Many quantities are estimates:

\[
R_k(\tau)
\pm
\sigma_k(\tau).
\]

Therefore PARS should distinguish:

### Nominal dominance

Based on point estimates.

Used for visualization and ranking.

### Certain dominance

Candidate \(a\) certainly dominates \(b\) only if, conservatively,

\[
U_k(a)
<
L_k(b)
\]

for all relevant objectives, modulo configured tolerances.

This means hard deletion should require much stronger evidence than merely:

\[
\hat R(a)<\hat R(b).
\]

Therefore:

\[
\boxed{
\text{uncertain dominated paths go to an archive, not to /dev/null.}
}
\]

This matters especially for MLP \(\rightarrow\) QM ranking inversion.

---

# 13. Adaptive Pareto margin

Strict Pareto membership is too brittle.

Define a practical directed distance from candidate \(\tau\) to the frontier:

\[
\Gamma(\tau)
=
\min_{\phi\in\mathcal F}
\max_k
\left[
\widetilde R_k(\tau)
-
\widetilde R_k(\phi)
\right]_+.
\]

Interpretation:

\[
\Gamma\approx0
\]

means the candidate lies on or close to the current frontier.

Larger \(\Gamma\) means it is more strongly dominated.

Candidates are divided into:

```text
ACTIVE FRONTIER
NEAR-FRONTIER RESERVE
COLD ARCHIVE
CERTAINLY DOMINATED
```

The active set for margin \(\gamma\) is

\[
\mathcal F_\gamma
=
\{\tau:\Gamma(\tau)\le\gamma\}.
\]

As expensive target information accumulates, \(\gamma\) may be increased.

Thus:

\[
\text{few QM labels}
\Rightarrow
\text{strong MLP/frontier prior},
\]

while

\[
\text{more QM evidence}
\Rightarrow
\text{larger rescue region}.
\]

This allows previously dominated MLP candidates to return if QM evidence shows that the MLP ranking is unreliable.

---

# 14. Pareto cover

Even the Pareto frontier may contain hundreds or thousands of trajectories.

Therefore frontier membership alone is insufficient.

PARS introduces an approximate **Pareto cover**:

\[
\mathcal C_t
\subset
\mathcal F_\gamma.
\]

For normalized objectives, a directed objective-space distance may be defined as

\[
d_+(a,b)
=
\max_k
[\widetilde R_k(a)-\widetilde R_k(b)]_+.
\]

The cover should retain a small number of representatives such that every relevant frontier region is represented within tolerance \(t\).

Conceptually:

\[
10000\ \text{paths}
\rightarrow
800\ \text{near-frontier paths}
\rightarrow
30\ \text{representatives}.
\]

Only the representatives need expensive refinement initially.

---

# 15. Objective-space cover alone is insufficient

Two paths can have nearly identical objective values but represent entirely different mechanisms.

For example:

\[
R(\tau_a)\approx R(\tau_b)
\]

while chemically,

```text
τa: proton transfer → bond cleavage
τb: bond cleavage → proton transfer
```

Therefore PARS must use a **two-level cover**.

## Level 1: objective-space Pareto cover

Compress redundant trade-offs.

## Level 2: reaction-path diversity

Within each covered objective region, preserve distinct mechanisms.

A path distance may contain:

\[
D_{\rm path}
=
D(
\text{event sequence},
\text{internal-coordinate trajectory},
\text{active atoms},
\text{geometry}
).
\]

Possible descriptors include:

- bond-change sequence,
- coordination-number trajectory,
- active-atom set,
- path fingerprint,
- symmetry-aware RMSD,
- path embedding,
- discrete Fréchet / Hausdorff-like path distance.

Thus the selected set is not merely

\[
\text{Pareto points},
\]

but

\[
\boxed{
\text{Pareto representatives of distinct reaction tubes}.
}
\]

---

# 16. Reaction tube construction

For channel \(c\):

\[
\mathcal T_c
=
\{\tau_1,\ldots,\tau_m\}.
\]

After frontier filtering:

\[
\mathcal F_{\gamma,c}
\subseteq
\mathcal T_c.
\]

Then cluster structurally:

\[
\mathcal F_{\gamma,c}
\rightarrow
\{
\mathcal U_{c1},
\mathcal U_{c2},
\dots
\}.
\]

Each tube \(\mathcal U_{cj}\) receives one or more representatives.

This produces:

```text
A → B
│
├── tube 1
│   ├── representative
│   └── alternate paths
│
├── tube 2
│   ├── representative
│   └── alternate paths
│
└── tube 3
```

The alternates are retained so that failed refinement of one representative does not destroy the mechanism.

---

# 17. Converting a PRRS trajectory into an NEB initial path

A raw PRRS trajectory should not be fed directly into NEB without preprocessing.

The conversion pipeline is:

\[
\text{DiscoveryTrace}
\rightarrow
\text{event window}
\rightarrow
\text{denoise}
\rightarrow
\text{align}
\rightarrow
\text{reparameterize}
\rightarrow
\text{NEB images}.
\]

## 17.1 Remove irrelevant trajectory sections

Remove:

- long pre-event oscillations,
- artificial forcing-only motion where appropriate,
- long post-product relaxation tails.

Keep provenance linking every path image back to its original frame.

---

## 17.2 Endpoint minimization

Reactant and product endpoints should be minimized on the same PES used for MLP path refinement.

---

## 17.3 Symmetry-aware alignment

The existing PRRS canonical/symmetry machinery should be used so that equivalent atoms do not create artificial path length.

---

## 17.4 Arc-length reparameterization

Define a path coordinate \(s\) and redistribute images approximately uniformly:

\[
s_i
=
\sum_{j<i}
d(R_j,R_{j+1}).
\]

The metric should preferentially use chemically meaningful internal coordinates for difficult systems.

---

## 17.5 Initial image count

V1 should use a conventional fixed or geometry-dependent image count.

Do **not** immediately couple Pareto selection to image insertion.

The first scientific question is whether Pareto geometry improves **which paths are refined**, not whether it improves the NEB integrator itself.

---

# 18. NEB refinement

For each selected representative:

\[
\tau_{\rm raw}
\rightarrow
\tau_{\rm NEB}.
\]

Recommended staged logic:

```text
endpoint relaxation
      ↓
ordinary NEB
      ↓
stable chain
      ↓
CI-NEB
      ↓
local saddle refinement
```

The calculation should monitor:

- maximum perpendicular force,
- chain continuity,
- image spacing,
- MLP reliability,
- topology consistency,
- pathological image crossing,
- convergence of the highest-energy region.

If MLP reliability becomes unacceptable in the transition region, the calculation should not silently continue and declare success.

Instead:

```text
status = UNRESOLVED_MLP
```

and the path is passed to a higher-fidelity route.

---

# 19. QM handoff

PARS does not assume that MLP path ordering is correct.

The core multi-fidelity relationship is:

\[
R_{\rm MLP}(\tau)
\quad\Longrightarrow\quad
R_{\rm QM}(\tau).
\]

A central empirical quantity is therefore the frequency of **ranking inversion**.

Example:

\[
\tau_a
\prec_{\rm MLP}
\tau_b
\]

but

\[
\tau_b
\prec_{\rm QM}
\tau_a.
\]

The workflow must explicitly measure this rather than assume it never happens.

Possible QM levels are:

```text
QM single points along selected images
QM gradients near the maximum
QM saddle refinement
QM Hessian
QM IRC
```

The amount of QM work depends on the scientific question.

---

# 20. QM batch selection

The QM budget should be allocated to several logically different purposes.

## Frontier representatives

Test paths representing important regions of the current Pareto frontier.

## Mechanism coverage

Ensure newly discovered tubes receive at least enough validation to determine whether they are physically credible.

## Near-frontier rescue

Sample candidates with small but non-zero

\[
\Gamma(\tau)
\]

to estimate how often the MLP incorrectly pushes good paths off the frontier.

## Reliability probes

Target trajectories where:

- model uncertainty is high,
- models disagree,
- forces are unusual,
- geometry leaves the training domain.

## Audit samples

Occasionally inspect candidates outside the preferred set.

Without an audit set, the algorithm cannot estimate how often Pareto pruning discards unexpected QM successes.

This is important scientifically.

---

# 21. Feedback from QM to PARS

QM results update the archive.

Suppose repeated inversions occur in one region of reaction space.

Then PARS may:

1. enlarge the adaptive margin \(\gamma\),
2. rescue archived trajectories,
3. increase sampling of the affected mechanism,
4. flag the MLP reliability model,
5. optionally generate training data for later potential improvement.

Thus the architecture becomes:

\[
\text{PRRS}
\rightarrow
\text{PARS}
\rightarrow
\text{QM}
\rightarrow
\text{PARS update}.
\]

Eventually it becomes an active multi-fidelity loop rather than a one-way pipeline.

---

# 22. Markov bridge integration

Markov bridge / conditioned trajectory sampling belongs **after endpoints or basins are known**.

It should therefore complement PRRS rather than replace it.

PRRS:

\[
A
\rightarrow
\text{discover } B.
\]

Bridge sampling:

\[
A,B
\rightarrow
\text{sample an ensemble conditioned on } A\rightarrow B.
\]

Then:

\[
\{\tau_i^{\rm bridge}\}
\rightarrow
\text{PARS frontier}
\rightarrow
\text{reaction tubes}
\rightarrow
\text{refinement}.
\]

This is especially attractive when a single deterministic NEB path hides a broad transition-path ensemble.

The long-term architecture can therefore become:

```text
PRRS
  │
  ├── discover products
  │
  ▼
A ───────────────► B
                    │
                    ▼
          conditional path sampling
          / Markov bridge / WE
                    │
                    ▼
             path ensemble
                    │
                    ▼
               PARS
                    │
                    ▼
          representative tubes
                    │
                    ▼
             NEB / string
```

---

# 23. Multi-state / spin-crossing systems

Ordinary PARS assumes that a meaningful scalar PES exists for the path being refined.

This is unsafe for systems involving:

- Fe/heme chemistry,
- near-degenerate spin states,
- redox transitions,
- strong nonadiabatic character.

Such paths require an explicit state flag:

```text
electronic_state_status:
    ordinary
    unknown
    suspected_crossing
    multistate_required
```

If a crossing is suspected:

\[
\boxed{
\text{do not force the problem into ordinary single-surface NEB}.
}
\]

Instead the path becomes a handoff candidate for Stateweaver / MECP-style treatment.

MLPATH should communicate with Stateweaver through an adapter rather than importing its implementation directly.

Thus:

```text
MLPATH/PARS
    │
    ├── ordinary path → NEB/string
    │
    └── multistate warning → Stateweaver
```

---

# 24. V2: Pareto-adaptive image placement

This is deliberately **not part of the first implementation**.

Once path-level selection is validated, Pareto ideas can be applied inside a refined path.

For image \(i\):

\[
Q_i=
\left(
E_i,
|F_{\perp,i}|,
\kappa_i,
U_i,
S_i
\right).
\]

Potential insertion signals include:

- high energy curvature,
- high perpendicular-force residual,
- high geometric curvature,
- high MLP uncertainty,
- rapidly changing electronic-state character.

Images could then be concentrated around physically difficult portions of the path.

Conceptually:

\[
\rho_{\rm image}(s)
=
f(
E,
F_\perp,
\kappa,
U,
S
).
\]

However this must be implemented carefully because NEB tangent construction and spring forces already depend on image spacing.

Therefore V2 should begin as:

```text
converge
→ inspect local path complexity
→ insert/remove images
→ restart/refine
```

rather than dynamically changing image density every optimization step.

---

# 25. Integration with the existing MLPATH codebase

The current repository already contains:

```text
src/prrs/
    canonical.py
    chemistry.py
    completeness.py
    config.py
    internal.py
    network.py
    perturbations.py
    reliability.py
    runner.py
    schema.py
    search.py
    state.py
    validation.py
    ...
```

PARS should initially be added **without refactoring the working PRRS core**.

The first implementation can remain consistent with the current relatively flat layout.

Recommended V1 additions:

```text
src/prrs/
    pathways.py
    frontier.py
    path_metrics.py
    refinement.py
```

with responsibilities:

### `pathways.py`

```text
DiscoveryTrace → PathCandidate
channel grouping
mechanism fingerprint
reaction-tube construction
trajectory trimming
path resampling
```

### `path_metrics.py`

```text
cheap path objectives
post-NEB objectives
normalization
path distance
tube distance
```

### `frontier.py`

```text
Pareto dominance
uncertainty-aware dominance
Pareto archive
frontier margin
Pareto cover
representative selection
```

### `refinement.py`

```text
generic refinement interface
NEB adapter
CI-NEB adapter
string adapter
refinement result schema
```

This avoids prematurely turning the repository into many nested packages.

---

# 26. Mature package structure

If PARS becomes large enough, the four modules can later be split without changing the public API:

```text
src/prrs/
├── canonical.py
├── chemistry.py
├── network.py
├── perturbations.py
├── reliability.py
├── search.py
├── state.py
│
├── pathways/
│   ├── candidate.py
│   ├── extract.py
│   ├── fingerprint.py
│   ├── cluster.py
│   └── resample.py
│
├── frontier/
│   ├── objectives.py
│   ├── dominance.py
│   ├── margin.py
│   ├── cover.py
│   ├── archive.py
│   └── selection.py
│
├── refinement/
│   ├── base.py
│   ├── neb.py
│   ├── string.py
│   └── saddle.py
│
└── multifidelity/
    ├── qm.py
    ├── audit.py
    └── stateweaver.py
```

The generic Pareto mathematics should remain independent of chemistry-specific code.

If another real project later reuses it, the mathematical core can then be extracted into a standalone library.

Until that happens:

\[
\boxed{
\text{keep Pareto geometry inside MLPATH.}
}
\]

---

# 27. Schema additions

A possible conceptual schema is:

```python
PathCandidate:
    id
    source_trace_id

    reactant_state_id
    product_state_id
    channel_id
    tube_id

    frames
    energies
    forces

    perturbation_metadata
    event_fingerprint

    objectives_raw
    objectives_normalized
    objective_uncertainty

    pareto_status
    pareto_margin
    cover_id

    reliability_status
    electronic_state_status

    refinement_status
    validation_status

    provenance
```

Important rule:

\[
\boxed{
\text{Every compressed/processed path remains traceable to the original PRRS trial.}
}
\]

No representative path should lose its provenance.

---

# 28. Search integration point

The existing PRRS loop conceptually becomes:

```python
results = run_prrs(state)

classified = classify_responses(results)

network.update(classified)
```

PARS extends this to:

```python
results = run_prrs(state)

classified = classify_responses(results)

network.update(classified)

candidates = build_path_candidates(classified)

for channel in group_by_channel(candidates):

    objectives = evaluate_path_objectives(channel)

    archive = build_pareto_archive(
        channel,
        objectives,
    )

    active = archive.with_margin(gamma)

    tubes = build_reaction_tubes(active)

    representatives = pareto_cover(
        tubes,
        tolerance=cover_tolerance,
    )

    refined = refine_selected_paths(representatives)

    validation_queue.add(refined)
```

The key architectural point is that `search.py` should orchestrate the operation but **not contain the Pareto implementation**.

---

# 29. Reaction network integration

The reaction network should store more than:

\[
A\leftrightarrow B.
\]

Each edge should eventually be able to contain:

```text
A → B
    channel metadata
    discovered paths
    reaction tubes
    selected representative paths
    MLP barriers
    QM barriers
    TS IDs
    validation level
    reliability flags
```

Conceptually:

\[
G_{\rm reaction}
=
(V,E)
\]

with

\[
E_{AB}
=
\{
\mathcal U_1,
\mathcal U_2,\ldots
\}.
\]

Thus different mechanisms become multiple pieces of evidence for the same network edge rather than being overwritten.

---

# 30. CLI / configuration

A possible configuration block is:

```yaml
pars:
  enabled: true

  objectives:
    pre:
      - energy_excursion
      - trigger_severity
      - mlp_uncertainty
      - path_curvature
      - recrossing
      - endpoint_penalty

    refined:
      - mlp_barrier
      - max_perpendicular_force
      - neb_residual
      - path_curvature
      - mlp_uncertainty

  normalization:
    method: robust

  dominance:
    mode: uncertainty_aware

  margin:
    enabled: true
    gamma: 0.15

  cover:
    enabled: true
    objective_tolerance: 0.10
    preserve_mechanism_diversity: true

  refinement:
    method: neb
    climbing_image: true
```

Exact defaults must be determined prospectively rather than selected after examining benchmark outcomes.

---

# 31. Tests

New unit tests should include:

```text
tests/test_frontier.py
tests/test_frontier_uncertainty.py
tests/test_pareto_margin.py
tests/test_pareto_cover.py
tests/test_path_metrics.py
tests/test_path_fingerprint.py
tests/test_reaction_tubes.py
tests/test_path_resample.py
tests/test_refinement_contract.py
```

Important invariants include:

### Permutation invariance

Candidate ordering must not alter the frontier.

### Scale invariance after normalization

Changing physical units must not alter selection.

### Symmetry invariance

Symmetry-equivalent trajectories should remain equivalent.

### Duplicate stability

Duplicating one trajectory 100 times must not change mechanism coverage.

### Channel isolation

A superior path in

\[
A\rightarrow B
\]

must not eliminate

\[
A\rightarrow C.
\]

### Archive rescue

Increasing \(\gamma\) must monotonically expand or preserve the active set.

### Deterministic cover

Given identical input and seed/configuration, representative selection must be reproducible.

---

# 32. Benchmark strategy

PARS needs its own benchmark question.

The question is not:

> Can PRRS find a reaction?

That has already been tested elsewhere.

The PARS question is:

\[
\boxed{
\text{Given many candidate paths, can PARS reduce expensive refinement cost}
\\
\text{without losing chemically important mechanisms or high-quality paths?}
}
\]

---

# 33. Baselines

PARS should be compared against:

1. random path selection,
2. lowest MLP energy excursion,
3. lowest apparent MLP barrier,
4. uncertainty-only selection,
5. path clustering only,
6. manually weighted scalar score,
7. strict Pareto frontier,
8. Pareto + adaptive margin,
9. Pareto + adaptive margin + reaction-tube cover.

The intended result is **not** merely:

> Pareto gives a prettier scatter plot.

It must demonstrate improved expensive-computation efficiency.

---

# 34. Primary evaluation metrics

Useful metrics include:

## Channel retention

\[
\text{fraction of known/discovered reaction channels retained}.
\]

## Mechanism retention

How many distinct validated reaction tubes survive compression?

## Oracle-path retention

If exhaustive reference calculations are available, does the selected subset contain the eventual best/refinable path?

## QM efficiency

\[
\frac{
\text{validated useful TS / mechanisms}
}{
\text{QM calls}
}.
\]

## Compression ratio

\[
\frac{
N_{\rm raw}
}{
N_{\rm refined}
}.
\]

## Ranking inversion rate

How often does

\[
\text{MLP preference}
\]

disagree with

\[
\text{QM preference}?
\]

## Duplicate-refinement rate

How much compute is wasted refining paths that turn out to represent the same mechanism?

---

# 35. Best initial benchmarks

The existing MLPATH benchmark ladder already contains simple reactions useful for basic correctness, but PARS particularly needs **competing paths**.

A sensible progression is:

### Malonaldehyde

Use as a sanity test.

Expected:

- many geometrically redundant proton-transfer trajectories,
- heavy compression,
- same mechanism retained.

### SN2

Tests:

- intermolecular geometry,
- ionic reaction,
- alignment,
- path representation.

### Competing substitution/elimination system

This should become the first serious PARS benchmark because there are genuinely competing channels.

The key question becomes:

\[
\boxed{
\text{Does channel-conditioned Pareto selection preserve both SN2 and E2-like behavior?}
}
\]

### Larger multi-path system

Later test:

- multiple conformers,
- multiple paths to same product,
- MLP uncertainty,
- substantial path redundancy.

---

# 36. Prospective acceptance criteria

As with the existing PRRS benchmark philosophy, criteria should be frozen before running the final evaluation.

Example questions:

- minimum required channel recall,
- maximum allowed loss of validated mechanisms,
- target compression ratio,
- maximum QM budget,
- allowable stochastic variation,
- definition of a successful TS,
- required IRC verification level.

The benchmark should not be retuned after observing which PARS settings happen to win.

---

# 37. Critical ablation: scalar score vs Pareto geometry

One of the most important experiments is:

\[
S(\tau)
=
\sum_k w_k R_k(\tau)
\]

versus

\[
R(\tau)
=
(R_1,\ldots,R_K).
\]

The scalar method forces one global trade-off.

PARS instead preserves multiple valid compromises.

The scientific hypothesis is:

> Reaction-path search contains several physically meaningful trade-offs that should not be collapsed into a single arbitrary weight vector before expensive validation.

This needs empirical demonstration.

---

# 38. Critical ablation: strict frontier vs adaptive margin

Strict Pareto pruning assumes that cheap objectives preserve enough of the true ordering.

That assumption is likely to fail occasionally when moving from MLP to QM.

Therefore compare:

```text
strict Pareto
vs
Pareto + margin
vs
full candidate set
```

Measure:

\[
P(
\text{QM-useful path was MLP-dominated}
).
\]

If this probability is non-negligible, the adaptive archive is justified empirically.

---

# 39. Critical ablation: objective cover vs mechanism-aware cover

Compare:

```text
objective-space cover only
```

against

```text
objective-space cover
+
path/mechanism diversity.
```

The expectation is that objective-only compression may merge distinct paths with similar scores.

This experiment establishes whether the concept of a **reaction tube** is necessary.

---

# 40. Relationship to the Pareto-cover theory

The motivating Pareto-cover work provides a useful conceptual framework:

\[
\text{large candidate family}
\rightarrow
\text{Pareto geometry}
\rightarrow
\text{small representative subset}
\rightarrow
\text{few expensive target labels}.
\]

However, its theoretical guarantees must **not** be copied directly to reaction-path search.

MLPATH violates several idealized assumptions:

- candidate trajectories are correlated,
- objectives are estimated rather than known exactly,
- MLP and QM ordering may invert,
- paths possess chemical/mechanistic structure beyond objective vectors,
- the candidate population evolves adaptively.

Therefore PARS should treat the theory as:

\[
\boxed{\text{design inspiration, not an inherited theorem.}}
\]

A reaction-path-specific theory can be developed later if the empirical method succeeds.

---

# 41. What PARS changes conceptually about NEB

Traditional workflow:

\[
A+B
\rightarrow
\text{guess one path}
\rightarrow
\text{NEB}.
\]

PARS workflow:

\[
A
\rightarrow
\text{discover many possible }B_i
\rightarrow
\text{generate many paths}
\rightarrow
\text{organize reaction space}
\rightarrow
\text{select representative tubes}
\rightarrow
\text{NEB}.
\]

Therefore NEB changes role.

It is no longer the reaction-discovery engine.

It becomes:

\[
\boxed{
\text{a local path-refinement operator applied to selected reaction tubes}.
}
\]

This is the central architectural shift.

---

# 42. Complete target workflow

The final intended pipeline is:

\[
\boxed{
\begin{array}{c}
\text{Conformer reservoir}\\
\downarrow\\
\text{PRRS perturbation proposals}\\
\downarrow\\
\text{short MLP responses}\\
\downarrow\\
\text{chemical-state classification}\\
\downarrow\\
\text{reaction-channel discovery}\\
\downarrow\\
\text{PathCandidate extraction}\\
\downarrow\\
\text{cheap multi-objective characterization}\\
\downarrow\\
\text{channel-conditioned Pareto archive}\\
\downarrow\\
\text{adaptive margin}\\
\downarrow\\
\text{mechanism/tube clustering}\\
\downarrow\\
\text{Pareto cover}\\
\downarrow\\
\text{representative-path refinement}\\
\downarrow\\
\text{NEB / string / saddle refinement}\\
\downarrow\\
\text{QM validation / Stateweaver handoff}\\
\downarrow\\
\text{frontier update}\\
\downarrow\\
\text{validated reaction network}
\end{array}
}
\]

---

# 43. Minimal viable implementation

The first version should remain deliberately small.

## MVP-1 — path representation

Implement:

```text
PathCandidate
trajectory trimming
channel grouping
event/path fingerprint
path resampling
```

No Pareto optimization yet.

---

## MVP-2 — static frontier

Implement:

```text
objective extraction
normalization
Pareto dominance
frontier visualization
```

Test using stored existing PRRS trajectories.

---

## MVP-3 — margin and archive

Implement:

```text
certain dominance
near-frontier margin
cold archive
rescue operation
```

---

## MVP-4 — mechanism-aware cover

Implement:

```text
path distance
reaction-tube clustering
objective-space cover
representative selection
```

This is the first point where major computational savings should become visible.

---

## MVP-5 — NEB adapter

Implement:

```text
PathCandidate
→ resampled initial chain
→ MLP NEB
→ RefinedPath
```

Keep the refinement interface backend-independent.

---

## MVP-6 — QM feedback experiment

Run a fixed number of QM checks.

Measure:

```text
ranking inversion
frontier stability
margin rescue
QM efficiency
```

This is the key experimental validation of PARS.

---

# 44. What should explicitly NOT be done in V1

Do not immediately attempt:

- a universal Pareto optimization library,
- fully dynamic NEB image redistribution,
- automatic MLP retraining,
- end-to-end differentiable reaction search,
- full Markov bridge integration,
- multistate MECP inside the same engine,
- a giant learned path encoder.

Those are plausible later extensions but would prevent identifying which part of PARS actually works.

V1 should test one clean hypothesis:

\[
\boxed{
\text{Can Pareto geometry + mechanism-aware covering reduce}
\\
\text{reaction-path refinement cost without destroying reaction diversity?}
}
\]

---

# 45. Long-term architecture

If the method succeeds, MLPATH evolves from

\[
\text{reaction discovery software}
\]

into

\[
\boxed{
\text{a hierarchical reaction-space exploration framework}.
}
\]

The levels become:

\[
\text{configuration space}
\rightarrow
\text{basins}
\rightarrow
\text{reaction channels}
\rightarrow
\text{reaction tubes}
\rightarrow
\text{MEPs / saddles}
\rightarrow
\text{electronic-state validation}.
\]

PRRS provides exploration.

PARS provides compression and decision-making.

NEB/string provides geometric refinement.

QM provides high-fidelity validation.

Stateweaver handles cases where a single adiabatic surface is no longer adequate.

---

# 46. Core methodological statement

The method can ultimately be summarized as:

> **PARS treats reaction-path search as a multi-objective, multi-fidelity selection problem over an ensemble of reaction trajectories rather than as the optimization of one preselected chain. PRRS first discovers candidate reaction channels without requiring known products. PARS then constructs channel-conditioned Pareto archives, preserves mechanistically distinct reaction tubes, and selects a compact representative cover for expensive NEB/string and QM refinement. Near-frontier candidates remain recoverable through an adaptive margin so that errors in the machine-learning potential do not cause irreversible path deletion.**

Or mathematically:

\[
\boxed{
\text{Discover}
\rightarrow
\text{Partition}
\rightarrow
\text{Frontier}
\rightarrow
\text{Cover}
\rightarrow
\text{Refine}
\rightarrow
\text{Validate}
\rightarrow
\text{Update}.
}
\]

This should be the central design logic of **Pareto-Adaptive Reaction Path Search**.