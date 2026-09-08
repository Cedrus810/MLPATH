# MLP-Native Reaction Discovery via Active Perturbation–Response Search

## 1. Core Idea

The goal is not to build another small modification of NEB, but to design a genuinely **MLP-native reaction discovery framework**.

The central problem is defined as:

\[
\boxed{
\text{Starting from a stable structure } A,
\text{ actively perturb the system, observe its response,}
\text{ and automatically discover reaction channels } A \rightarrow B_i
}
\]

The method does **not** require:

- a full Hessian,
- a predefined product,
- a predefined reaction coordinate,
- or a predefined minimum-energy path.

A temporary name for the framework is:

> **Perturbation–Response Reaction Search (PRRS)**

The core loop is:

\[
\boxed{
\text{State}
\rightarrow
\text{Perturb}
\rightarrow
\text{Respond}
\rightarrow
\text{Classify}
\rightarrow
\text{Refine}
}
\]

More explicitly,

\[
R_0
\rightarrow
\{P_k\}
\rightarrow
\{\tau_k(t)\}
\rightarrow
\{B_1,B_2,\ldots\}
\rightarrow
\text{reaction channels}.
\]

Here,

- \(R_0\): current local minimum,
- \(P_k\): the \(k\)-th perturbation,
- \(\tau_k(t)\): a short MLP trajectory,
- \(B_i\): newly discovered basin or connectivity state.

The first target is therefore not a transition state, but:

\[
\boxed{\text{basin connectivity}}
\]

namely,

\[
A \leftrightarrow B.
\]

Once connectivity is discovered, barrier and transition-state refinement can be performed afterward.

---

# 2. Chemical-State Representation

A direct search in the full Cartesian \(3N\)-dimensional space is unnecessarily inefficient.

Instead, the current configuration should first be represented by a dynamic molecular graph:

\[
G(R)=\{V,E,W\}.
\]

The graph edges should preferably be continuous rather than purely binary:

\[
w_{ij}=f(r_{ij},Z_i,Z_j),
\]

instead of simply

\[
w_{ij}\in\{0,1\}.
\]

This makes it possible to describe:

- bond stretching,
- emerging contacts,
- coordination changes,
- proton transfer,
- ring formation,
- bond breaking.

A local environment embedding may also be extracted:

\[
z_i=\Phi(R)_i.
\]

For an MLP such as MACE, this representation could potentially reuse internal hidden representations rather than requiring a separate molecular encoder.

The complete state may be written as:

\[
S(R)=\{E,F,G,z\}.
\]

---

# 3. Perturbation Generator

The perturbation generator is the central component of the method.

The search should not rely solely on isotropic random Cartesian noise.

At least four perturbation families are useful.

---

## 3.1 Geometry Perturbation

Examples include:

\[
r_{ij}\rightarrow r_{ij}+\delta.
\]

Possible operations include:

- bond stretching,
- bond compression,
- angle distortion,
- dihedral rotation,
- nonbonded atom-pair approach,
- atom-pair separation.

Importantly, perturbations should not be restricted to existing bonds.

For a nonbonded pair \(i,j\),

\[
r_{ij}\downarrow
\]

can itself be interpreted as a possible bond-formation probe.

---

## 3.2 Force Perturbation

A temporary external force can be added:

\[
F_i'=F_i+F_i^{\rm ext}.
\]

For an atom pair,

\[
F_i^{\rm ext}=+\alpha \hat r_{ij},
\]

\[
F_j^{\rm ext}=-\alpha \hat r_{ij}.
\]

This preserves the total translational force on the pair.

The perturbation may be applied only for a short time,

\[
\tau_p \sim 5-50\ {\rm fs},
\]

and then removed.

The artificial force is therefore not responsible for generating the entire reaction path. Its purpose is only to **probe the local reactive response of the system**.

---

## 3.3 Momentum Kick

A particularly MLP-native approach is to perturb velocities rather than positions.

Inject a local amount of kinetic energy:

\[
\Delta K.
\]

For example,

\[
v_i'=v_i+\Delta v_i.
\]

The perturbation should satisfy

\[
\sum_i m_i\Delta v_i=0,
\]

and preferably also

\[
\sum_i r_i\times m_i\Delta v_i=0,
\]

to avoid introducing overall translation or rotation.

The system is then propagated under short NVE dynamics:

\[
R_0,V_0+\Delta V_p
\rightarrow
R(t).
\]

Typical trajectory lengths may be on the order of

\[
10-500\ {\rm fs}.
\]

This directly probes the nonlinear dynamical response of the potential-energy surface.

---

## 3.4 Collective Perturbation

Complex reactions generally involve multiple coupled degrees of freedom.

A collective perturbation can therefore be defined as

\[
P=\sum_k c_k q_k,
\]

where \(q_k\) may include:

- pair distances,
- bond angles,
- coordination numbers,
- local displacement modes,
- learned latent directions.

Such a perturbation can simultaneously describe events such as

\[
\text{bond breaking}
+
\text{bond formation}
+
\text{proton transfer}.
\]

This is likely essential for reactions that cannot be represented by a single pair-distance coordinate.

---

# 4. Perturbation Amplitude as a Search Variable

Each perturbation direction \(P_k\) should be explored over a range of amplitudes:

\[
a_1<a_2<\cdots<a_m.
\]

Thus,

\[
P_{k,m}=a_m P_k.
\]

The resulting response becomes

\[
\mathcal R(P,a).
\]

A typical behavior could be

\[
a<a_c:
\qquad
A\rightarrow A,
\]

while

\[
a>a_c:
\qquad
A\rightarrow B.
\]

The critical perturbation amplitude

\[
a_c
\]

is not a rigorous activation free energy,

\[
a_c\neq \Delta G^\ddagger,
\]

but it can act as a useful **reaction-accessibility or reaction-susceptibility score**.

For example,

\[
S_{\rm react}=\frac{1}{a_c}.
\]

A smaller \(a_c\) indicates that a reaction channel is easier to trigger under the chosen perturbation family.

---

# 5. Perturbation–Response Experiment

Each perturbation becomes a small numerical experiment:

\[
R_0
\xrightarrow{P_k}
R_k^*
\xrightarrow[\rm MLP]{MD}
R_k(t).
\]

The trajectory can be separated into two phases.

## Phase I: Forced Response

For

\[
0<t<t_p,
\]

the perturbation remains active.

## Phase II: Free Relaxation

For

\[
t_p<t<t_f,
\]

the perturbation is removed.

The second stage is crucial.

If the system returns to the original basin,

\[
R(t)\rightarrow A,
\]

the response is mainly elastic or conformational.

If instead

\[
R(t)\rightarrow B,
\]

a basin transition has been discovered.

---

# 6. Analyze the Full Response, Not Only the Final Structure

Each short trajectory should be represented by a response descriptor

\[
X_k=
\{
\Delta E(t),
F(t),
G(t),
r_{ij}(t),
CN_i(t),
z_i(t)
\}.
\]

Important observables include:

- energy evolution,
- force evolution,
- dynamic bond graph,
- interatomic distances,
- coordination numbers,
- latent environment embeddings.

A particularly useful event is

\[
\Delta G(t)\neq 0,
\]

indicating a chemical-topology change.

Examples include:

- C–O bond breaking,
- C–N bond formation,
- proton transfer,
- coordination change.

Each response can then be compressed into a reaction fingerprint:

\[
\mathcal F_k=
\{
E_{\rm broken},
E_{\rm formed},
\Delta CN,
\text{active atoms}
\}.
\]

Large numbers of trajectories may then be clustered into a small number of reaction channels:

\[
A\rightarrow B_1,
\]

\[
A\rightarrow B_2,
\]

\[
A\rightarrow B_3.
\]

---

# 7. Response Classification

A minimal response classification may include five classes:

\[
C=
\{
\text{elastic},
\text{rearrangement},
\text{reactive},
\text{dissociation},
\text{unphysical}
\}.
\]

## Elastic

\[
A\rightarrow A.
\]

The system returns to the original state.

## Rearrangement

The conformation changes substantially but the chemical graph remains unchanged.

## Reactive

\[
G_A\neq G_B.
\]

A genuine connectivity change occurs.

## Dissociation

One or more fragments separate.

## Unphysical

Examples include:

- excessively large forces,
- atom overlap,
- unrealistic short distances,
- obvious MLP extrapolation,
- numerical instability.

The last category is particularly important when using universal MLPs.

---

# 8. Adaptive Perturbation Search

The first search round may deliberately be broad.

For example,

\[
N=128-1024
\]

perturbations can be tested in parallel.

Subsequent rounds should use the observed response to decide which perturbations are most informative.

Define

\[
p(B|A,P),
\]

the probability that perturbation \(P\), applied to basin \(A\), reaches basin \(B\).

The next perturbation can be selected through a utility function:

\[
P^*=\arg\max_P U(P).
\]

One possible utility is

\[
U(P)
=
w_1P_{\rm reaction}
+w_2N_{\rm novel}
-w_3E_{\rm perturb}
-w_4U_{\rm MLP}.
\]

Here,

- \(P_{\rm reaction}\): estimated probability of triggering a reaction,
- \(N_{\rm novel}\): novelty of the resulting state,
- \(E_{\rm perturb}\): perturbation cost,
- \(U_{\rm MLP}\): model uncertainty or extrapolation risk.

Thus, the search naturally evolves from broad exploration to targeted probing of reactive directions.

---

# 9. Coarse-to-Fine Search in Perturbation Strength

Perturbation amplitude can be refined adaptively.

Suppose an initial scan gives

\[
a=1:\quad A,
\]

\[
a=2:\quad B.
\]

A binary or Bayesian search can then refine the transition threshold:

\[
1
\rightarrow
1.5
\rightarrow
1.75
\rightarrow
\cdots
\]

until an approximate critical amplitude is obtained:

\[
a_c.
\]

Different reaction channels can then be ranked:

\[
a_c(A\rightarrow B)
<
a_c(A\rightarrow C).
\]

This indicates that \(A\rightarrow B\) is easier to activate under the chosen perturbation family.

---

# 10. Reaction Discovery Before Reaction-Path Refinement

Once a reactive trajectory

\[
A\rightarrow B
\]

has been found, a candidate transition region can be extracted from

\[
R(t).
\]

For example,

\[
R^\dagger_{\rm candidate}
=
\arg\max_t E(R(t)).
\]

However, this structure should not immediately be identified as a true transition state.

Instead, a sequence of configurations

\[
R(t_1),\ldots,R(t_m)
\]

can be used as a **data-informed initial path**.

Only at this point should conventional local refinement methods be introduced, such as:

- string methods,
- dimer methods,
- eigenvector following,
- NEB or CI-NEB,
- direct TS optimization.

In this framework, traditional reaction-path methods are no longer the discovery engine.

They become:

\[
\boxed{\text{local refinement tools}}
\]

after connectivity has already been identified.

---

# 11. Reaction Tubes Instead of a Single MEP

A more MLP-native extension is to avoid defining a unique reaction path too early.

Suppose multiple reactive trajectories are collected:

\[
\Gamma_1,\Gamma_2,\ldots,\Gamma_M.
\]

These define a conditional ensemble

\[
p(R|A\rightarrow B).
\]

A central reaction manifold may then be extracted:

\[
\Gamma_{\rm center}.
\]

Thus,

\[
\boxed{
\text{reaction path}
=
\text{central manifold of a reactive-trajectory ensemble}
}
\]

rather than a chain of artificially spring-connected images.

This naturally generalizes the concept of a reaction path into a **reaction tube**.

---

# 12. MLP Reliability Gate

The method must prevent an MLP from exploring arbitrarily far outside its reliable domain.

Each configuration should therefore be assigned an uncertainty or out-of-distribution score:

\[
U(R).
\]

Possible measures include ensemble force disagreement:

\[
U_F=
{\rm Var}
\left[
F_1,\ldots,F_M
\right],
\]

or embedding-space distance:

\[
U_z=
d(z,z_{\rm train}).
\]

If

\[
U(R)>U_{\max},
\]

the trajectory should be stopped or escalated to a higher-level calculation.

Conceptually:

```text
normal PES
   |
 perturb
   |
reactive region
   |
MLP uncertain?
  /       \
 no       yes
 |         |
continue   QM query
```

The reliability gate prevents the reaction search from mistaking MLP extrapolation artifacts for chemistry.

---

# 13. Natural Integration with Active Learning

The framework naturally becomes a reaction-oriented active-learning method.

Conventional active learning often follows:

\[
\text{MD}
\rightarrow
\text{uncertain configuration}
\rightarrow
\text{QM}.
\]

Here the acquisition rule can instead prioritize configurations that are simultaneously:

1. reaction-relevant, and
2. uncertain.

For example,

\[
A(R)
=
P_{\rm reactive}(R)
U_{\rm MLP}(R).
\]

Thus QM calculations are concentrated on the regions that matter most for reaction discovery.

This gives:

\[
\boxed{
\text{active learning focused on reaction-relevant uncertainty}
}
\]

rather than generic uncertainty sampling.

---

# 14. Full Reaction-Discovery Loop

The full algorithm can be summarized as

\[
\boxed{
\begin{aligned}
A
&\rightarrow
\text{chemical-state encoding}\\
&\rightarrow
\text{perturbation generation}\\
&\rightarrow
\text{batched short MLP trajectories}\\
&\rightarrow
\text{response analysis}\\
&\rightarrow
\text{reaction clustering}\\
&\rightarrow
\text{adaptive perturbation}\\
&\rightarrow
\text{reaction graph expansion}.
\end{aligned}
}
\]

For every newly discovered basin \(B_i\), repeat the same process:

\[
B_i\rightarrow\{C_j\}.
\]

The final output becomes a reaction network:

\[
\boxed{
\mathcal G_{\rm rxn}
=
(V_{\rm basin},E_{\rm reaction})
}
\]

where nodes represent metastable chemical states and edges represent discovered reaction channels.

---

# 15. Minimal Viable Prototype

The first implementation should remain deliberately small.

A useful MVP could include only three perturbation families:

1. nonbonded atom-pair compression,
2. existing-bond stretching,
3. local momentum kicks.

For each perturbation direction, test several amplitudes or initial kinetic-energy levels.

A practical initial scale could be:

\[
16-32
\]

perturbation amplitudes or velocity realizations per channel.

Each candidate can be propagated using a short MACE trajectory:

\[
50-200\ {\rm fs}.
\]

The simplest detection criterion is then:

\[
\text{bond graph before}
\rightarrow
\text{bond graph after}.
\]

The MVP should answer only one question:

\[
\boxed{
\text{Can the algorithm discover a known product without being given the product?}
}
\]

Suitable benchmark reactions could include:

- simple SN2 reactions,
- proton transfer,
- small-molecule rearrangements,
- simple pericyclic reactions.

The first benchmark should avoid unnecessarily complicated catalytic chemistry.

---

# 16. Learned Perturbation Proposal

After collecting initial data,

\[
(R,P)\rightarrow Y,
\]

where

\[
Y=
\{
\text{reaction class},
B,
a_c
\},
\]

a proposal network can be introduced:

\[
\pi_\theta(P|R).
\]

Its task is:

> Given the current configuration, propose the next perturbation that is most likely to reveal a useful reaction channel.

The responsibilities are then cleanly separated:

## MLP Potential

\[
R\rightarrow E,F.
\]

It represents the physical potential-energy surface.

## Proposal Model

\[
R\rightarrow P.
\]

It controls exploration.

This creates a self-improving reaction-discovery loop without requiring the proposal network itself to predict energies.

---

# 17. Key Design Principles

## Principle 1: Do Not Rebuild an Approximate Hessian

The method should avoid collapsing back into

\[
\text{quasi-Newton}+\text{MLP}.
\]

The purpose of perturbation is not primarily to reconstruct local second derivatives.

Instead,

\[
\boxed{\text{response itself is information}}
\]

and nonlinear response should be treated as the primary observable.

---

## Principle 2: Perturbation Is a Probe, Not Merely a Bias Potential

The algorithm should not reduce perturbation to another fixed biasing method.

The key loop is:

\[
\text{perturb}
\rightarrow
\text{observe}
\rightarrow
\text{learn}
\rightarrow
\text{adapt}.
\]

How the system reacts to a perturbation determines the next search step.

---

## Principle 3: Discover Connectivity Before Demanding a TS

The first task is:

\[
\boxed{
\text{Without knowing the product, can the algorithm discover reaction connectivity?}
}
\]

Only after this succeeds should the framework worry about:

- transition states,
- barriers,
- IRC,
- free energies,
- mechanistic ranking.

This cleanly separates **reaction discovery** from **reaction characterization**.

---

# 18. Methodological Interpretation

The framework can be summarized as:

\[
\boxed{
\textbf{Reaction discovery as active perturbation–response learning on an MLP potential-energy surface.}
}
\]

The key computational shift is that traditional electronic-structure methods treat every force evaluation as expensive, whereas MLPs allow extremely large numbers of cheap and highly parallel evaluations.

Traditional algorithms therefore optimize:

\[
\boxed{\text{minimize the number of force evaluations}}
\]

whereas an MLP-native algorithm can instead exploit:

\[
\boxed{\text{many cheap dynamical experiments}}
\]

to extract reaction information.

The central methodological opportunity is therefore not to force MLPs into algorithms designed for expensive ab initio calculations, but to design algorithms around the actual computational structure of modern MLPs:

\[
\boxed{
\text{massively parallel perturbation}
+
\text{short trajectory response}
+
\text{adaptive reaction discovery}
}
\]

---

# 19. Two Main Mathematical Problems

Once the overall framework is fixed, the two most important open design questions are:

## 19.1 Perturbation-Space Parameterization

How should

\[
P
\]

be represented?

Possible choices include:

- internal coordinates,
- atom-pair forces,
- velocity-space perturbations,
- graph-localized collective coordinates,
- MLP latent-space directions,
- learned low-dimensional perturbation operators.

The perturbation space should be expressive enough to trigger multibody chemistry while remaining structured enough to search efficiently.

---

## 19.2 Response Information Score

How should the value of an observed response be quantified?

A useful score may combine:

\[
\text{reactivity}
+
\text{novelty}
+
\text{low perturbation cost}
+
\text{model reliability}.
\]

For example,

\[
I(P)
=
w_r S_{\rm reaction}
+
w_n S_{\rm novelty}
-
w_p C_{\rm perturb}
-
w_u U_{\rm MLP}.
\]

A strong definition of this information score would provide the mathematical core of the adaptive search algorithm.

---

# 20. Compact Summary

The complete conceptual pipeline is:

```text
Initial basin A
      |
      v
Chemical-state encoding
      |
      v
Generate structured perturbations
      |
      v
Batched short MLP trajectories
      |
      v
Observe nonlinear response
      |
      +----> elastic / conformational
      |
      +----> reactive
      |
      +----> dissociation
      |
      +----> MLP-OOD
      |
      v
Cluster newly discovered products
      |
      v
Estimate reaction susceptibility
      |
      v
Adapt perturbation distribution
      |
      v
Expand reaction graph
      |
      v
Refine selected channels with TS/path methods
```

The central change in viewpoint is:

\[
\boxed{
\text{Do not search directly for a path.}
\quad
\text{Perturb the system and learn what reactions it is capable of performing.}
}
\]
