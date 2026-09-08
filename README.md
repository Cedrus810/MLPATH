# PRRS — Perturbation–Response Reaction Search

## What is PRRS?

**PRRS is a reaction-discovery method designed for machine-learning potentials (MLPs).**

The idea is simple:

> Perturb a molecule in a controlled way, remove the perturbation, let the system respond naturally, and observe which chemical state it reaches.

PRRS does **not** require:

- a predefined reaction coordinate,
- a predefined product,
- or an initial guess of the transition state.

Instead, it first discovers possible reaction channels and products. Transition states are then located and verified afterwards.

---

## Why do we need it?

Most reaction-search methods work best when we already know roughly **what reaction we are looking for**.

That is much harder when the reaction channel itself is unknown.

PRRS approaches the problem from the opposite direction:

1. start from a known molecular state,
2. systematically perturb physically meaningful internal coordinates,
3. let the system evolve freely on the MLP potential-energy surface,
4. relax the resulting structure,
5. determine whether it returned to the original state, found a new conformer, or reached a new chemical state.

Repeating this procedure allows us to gradually construct a reaction network without specifying the products in advance.

---

## Core idea

A PRRS trial can be summarized as:

**starting structure → controlled perturbation → free dynamics → relaxation → chemical-state identification**

There are three possible outcomes.

### 1. The system returns to the original minimum

The perturbation was essentially elastic.

This tells us that this direction does not easily lead to another basin at the tested amplitude.

### 2. The system reaches another minimum

PRRS has discovered either:

- another conformer of the same molecule, or
- a genuinely different chemical state.

These two cases are explicitly distinguished.

### 3. The search approaches a first-order saddle point

This becomes a transition-state candidate.

The candidate is then checked by following the downhill directions on both sides to determine which states the saddle connects.

---

## Why perturb internal coordinates?

Random Cartesian displacement is inefficient because most random motion has little chemical meaning.

PRRS instead perturbs coordinates such as:

- bond lengths,
- bond angles,
- torsions,
- collective internal motions,
- and intermolecular approach coordinates.

Early ethanol tests showed why this matters.

Using only radial perturbations repeatedly returned to the starting state even though several nearby conformational basins were known to exist. Adding torsional and collective motions allowed the search to reach the relevant parts of configuration space.

This is an important design principle:

> **The perturbation set must cover the chemically relevant degrees of freedom.**

---

## Chemical states and conformers are treated separately

A different minimum is not necessarily a different chemical species.

For example, bond rotation can create another conformer without changing chemical connectivity.

PRRS therefore uses two levels:

**Chemical state**

A distinct molecular connectivity and stereochemical state.

**Conformer reservoir**

Multiple geometries belonging to the same chemical state.

This distinction is important because reaction accessibility can depend strongly on the starting conformer.

A conformer that appears unreactive does not imply that every conformer of the same molecule is unreactive.

---

## Symmetry is handled explicitly

Chemically equivalent atoms should not generate artificial “new reactions.”

For example, transferring one of three equivalent methyl hydrogens must not create three different products simply because the atom labels are different.

PRRS therefore compares structures and reaction events while accounting for molecular symmetry.

This avoids:

- duplicated products,
- duplicated conformers,
- duplicated reaction channels,
- and artificial stereochemical distinctions caused by equivalent atoms.

---

# What has been demonstrated so far?

The method has progressed from conformational discovery to genuine reaction-network expansion.

## P0 — Ethanol

**Question:** Can the method correctly discover nearby conformational basins?

**Result:** Yes.

The test established the basic perturbation, relaxation, symmetry handling, and conformer-identification framework.

It also revealed that radial perturbations alone are insufficient and motivated the inclusion of torsional and collective perturbations.

---

## P1 — Malonaldehyde proton transfer

**Question:** Can PRRS discover a real reaction event and its transition state?

**Result:** Yes.

The proton transfer was reproducibly identified.

Because the two sides are chemically equivalent, the reaction is represented as a **degenerate self-loop** rather than as a new chemical node.

A subsequent QM calculation independently confirmed the transition-state geometry and reaction identity.

The PRRS and QM saddle geometries agree to approximately **0.015 Å RMSD**.

The MLP and QM barrier heights differ substantially, but this is an energy-model error rather than a failure to locate the reaction geometry.

---

## P2 — 3-Oxobutanal enol

**Question:** Can PRRS discover a genuinely different chemical state and expand the reaction network?

**Result:** Yes.

Across 16 independent starting seeds, the expected reaction channel was found in all 16 cases.

This case demonstrates transition from:

**single-state exploration → reaction-network expansion**

rather than simply conformational sampling.

---

## P3 — Benzoylacetone

**Question:** Does the approach still work in a more complicated system involving aromatic structure, chelation, and competing conformations?

**Result:** Yes under the predefined acceptance criteria.

The reaction channel and transition-state structure were successfully identified.

The geometry also shows the expected contraction of the donor–acceptor heavy-atom distance at the transition state.

---

## Independent QM checks

Two reaction cases have now been checked against QM calculations.

### Malonaldehyde

Bidirectional IRC calculations confirm that the PRRS transition state connects the expected minima.

The MLP and QM saddle geometries are nearly identical.

### SN2 reaction: Cl⁻ + CH₃Cl

A transition state generated from the PRRS search refined to essentially the same geometry as an independently located QM transition state.

Bidirectional IRC calculations also reached the expected products and reactants.

These tests support an important conclusion:

> **PRRS can identify the reaction event and transition-state geometry even when the MLP does not reproduce the QM barrier energy accurately.**

---

# What PRRS is — and is not — trying to predict

This distinction is central to the project.

## PRRS is responsible for

- discovering possible reaction channels,
- identifying products,
- locating candidate transition states,
- identifying which states a transition state connects,
- and building the reaction network.

## PRRS is not expected to provide final quantitative barrier energies

Absolute barrier heights and imaginary frequencies depend on the quality of the potential-energy surface.

The current workflow therefore separates two questions:

### Question 1 — Is the search algorithm good?

Compare the PRRS result with a reference search performed on the **same MLP potential-energy surface**.

This tests the search method itself.

### Question 2 — Is the MLP chemically accurate?

Compare the corresponding MLP transition state with a **QM calculation**.

This tests the potential-energy model.

Mixing these two comparisons would make it impossible to determine whether an error came from the search algorithm or from the MLP.

---

# Choice of machine-learning potential

The first development stages used **MACE-OFF24**.

This model worked adequately for the initial neutral organic systems but has an important limitation: it does not properly represent changes in total molecular charge.

This becomes unacceptable for ionic reactions such as SN2 chemistry.

For the next stage, the default model is therefore **MACE-POLAR-1-M**, which performed substantially better across both:

- neutral proton-transfer chemistry,
- and the ionic SN2 benchmark.

The earlier P0–P3 results remain frozen and are not retrospectively changed.

---

# Current limitations

The method is working, but several important problems remain.

## 1. Perturbing bond lengths can distort the rest of the molecule

The current bond-stretch perturbation can move the selected atoms while unintentionally deforming nearby bonds.

This is an artifact of how the perturbation is applied, not a necessary physical effect.

The next implementation step is therefore to move molecular fragments more rigidly when a bond is stretched.

---

## 2. Reaction-direction selection is still heuristic

PRRS currently generates several chemically motivated perturbation directions, but their priority is not yet learned or optimized.

A future version should rank perturbations according to their observed usefulness for:

- basin discovery,
- reaction discovery,
- and transition-state discovery.

---

## 3. Transition-state verification becomes expensive for larger molecules

For the current small test systems, rigorous curvature analysis is practical.

For larger systems, the same calculation will eventually become too expensive.

A scalable transition-state verification strategy will therefore be required before applying PRRS to substantially larger molecules.

---

## 4. Ionic chemistry is only entering formal validation

The SN2 system has passed the initial feasibility tests and the transition state has already been independently confirmed by QM.

However, the formal benchmark criteria and automated acceptance test are not yet frozen.

Therefore this case should currently be described as **promising but not yet a completed benchmark**.

---

# How results are evaluated

The project deliberately uses strict prospective evaluation.

Before a benchmark is run:

1. the scientific question is defined,
2. the acceptance criteria are frozen,
3. the calculation is then performed,
4. and the result is evaluated against those original criteria.

Parameters are not adjusted afterwards simply to turn a failed result into a successful one.

Incomplete calculations, numerical failures, and genuine negative results are also recorded separately.

This is intended to prevent benchmark design from gradually adapting to the desired answer.

---

# Current benchmark progression

The benchmark series is designed to add new chemical difficulty rather than repeatedly test the same reaction class.

| Stage | Chemical problem | Status |
|---|---|---|
| P0 | Conformer discovery | Completed |
| P1 | Degenerate intramolecular proton transfer | Completed + QM checked |
| P2 | Reaction-network expansion | Completed |
| P3 | More complex intramolecular proton transfer | Completed |
| b04 | Ionic intermolecular SN2 | Preflight successful + QM checked |
| b05 | Competing substitution / elimination | Next group |
| b06 | Electrophilic addition | Planned |
| b07 | Bond formation / fragmentation and Diels–Alder-type chemistry | Planned |

The next benchmarks are intentionally moving away from repeated proton-transfer examples toward:

- charged systems,
- intermolecular reactions,
- competing pathways,
- heavy-atom bond formation/breaking,
- and changes in fragment number.

---

# Next steps

The immediate priorities are:

### 1. Fix perturbation delivery

Bond stretching should move chemically connected fragments without introducing artificial internal distortion.

This must be solved before interpreting amplitude scans.

### 2. Complete the SN2 benchmark

Freeze formal acceptance criteria and evaluate the first charge-sensitive intermolecular reaction.

### 3. Expand chemical diversity

Proceed to systems with:

- substitution versus elimination competition,
- polar addition,
- bond formation and breaking,
- and changes in molecular fragment number.

### 4. Measure search efficiency

The final efficiency metric should not simply be trajectories per second.

The scientifically useful quantity is:

> **computational cost per independently QM-confirmed reaction channel discovered**

---

# Current conclusion

PRRS has already demonstrated that controlled perturbation followed by unbiased molecular response can discover:

- conformational basins,
- reaction events,
- new chemical states,
- and transition-state geometries,

without requiring a predefined reaction coordinate or product.

The strongest evidence so far is that independently performed QM calculations reproduce the reaction identity and transition-state geometry in the tested cases.

The main open question is no longer whether the basic idea works.

The next question is:

> **How broadly and efficiently can it discover chemically diverse reaction channels as system size and reaction complexity increase?**