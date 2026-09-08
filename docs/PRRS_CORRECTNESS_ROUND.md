# Correctness round: review findings, and what the formal framing changed

Implementation: `638c552d01341ab5` (was `158f2f3b35200290` at the start; round one
ended at `b054dcb0b763d3c1`).
Verified: `212 passed / 1 skipped`; digest regression 332 files / 2136 frames unchanged;
landscape regression 300/300 conformers, verdict fields unchanged, drift at machine
precision; `confirm_minimum` exercised on two MACE models through `sn2_probe` and on the
P3 proton transfer through `ts`.

Correction to what round one recorded here: the landscape regression line originally said
300/300 bit-identical. At the time it was written the actual output was `263/300
conformers present  MISSING [7, 15, 23, ...]` and `EXIT=1` -- shard 7 had died on CUDA OOM
and its 37 conformers had never been computed. The `263/263` line for unchanged verdict
fields was read as the conformer count. It is 300/300 now because shard 7 was recomputed,
not because it was then.

This round started from a code review listing ten defects and was then re-ordered by a
formal statement of the method. The framing turned out to predict the defect classes,
which is the main reason it is recorded here rather than only in the commit log:

| defect | which part of the framing it violates |
|---|---|
| `_projected_direction_from` removed the uniform vector, not the mass-weighted translation | the mass metric |
| `internal.displace` walked a chord where the coordinate is an arc (fixed earlier) | the perturbation map |
| "the rigid-body residual bounds the error of this whole calculation" | discovery is not certification |

A framing that only restates the code is decoration. This one located three defects of
three different kinds, so it is load-bearing.

---

## The method, as the corrections use it

Internal coordinates `q_alpha` are the probe language, and deliberately not a basis for
the whole `3N` space -- they are a chemistry-informed set. The mass metric decides how to
move one at least cost: for a coordinate with gradient `g`, the displacement achieving
`dq` with minimal `dR^T M dR` is

    dR = (dq / (g^T M^-1 g)) M^-1 g

which is what `internal.displace` computes (`direction = g/m`, `slope = sum(g*direction)`,
`increment/slope`), and the momentum kick is the same idea in velocity space with the
cross term retained so the delivered energy is exact.

Perturb-and-release is the discovery step. Quench and curvature then explain what was
found. These are different roles with different standards, and conflating them is what
produced the third defect above.

Curvature at vanishing amplitude certifies a candidate; finite-amplitude response
characterizes it. For a general probe map `Ret(av) = R + a v + a^2 c_v / 2`,

    kappa_0(v) = v^T H v + g . c_v

so the Hessian quotient is the `a -> 0` limit of perturbation response, exactly, whenever
`g . c_v` vanishes. That is a theorem about the two branches rather than an analogy, and
it has a consequence the code now states: **a Cartesian chord has `c_v = 0`, so its
`kappa_0` is `v^T H v` at any geometry, stationary or not.** The chord is the correct map
for a curvature. A curved retraction -- the arc of a torsion -- is the correct map for the
path the search travels, and its `kappa_0` differs by `g . c_v` wherever the gradient does
not vanish. The two desiderata genuinely conflict at finite amplitude; the code picks the
chord and now says so, and measures the bill.

---

## What was fixed

### The mass metric

`_projected_direction_from` (`runner.py`) subtracted an unweighted column mean, which
projects out the uniform vector instead of the mass-weighted translation
`(T_a)_ib = sqrt(m_i) delta_ab`. The two coincide only at equal masses. Measured
consequence: on water, a **pure x translation** came back as a full unit vector
`[+0.8165, -0.4082, -0.4082]` -- the output was entirely artefact, not a small residual.
Its only caller is `align_mode_with_travel`, which decides a relay's sign, so the error
fed a sign decision directly.

Fixed by removing the centre of mass in Cartesian coordinates before weighting, which is
the same projection written so its correctness is visible.

`lowest_response_direction` measured its convergence residual on the **unprojected**
gradient while stepping along the projected one, so the walk was charged for a component
it had deliberately excluded and could report a stall it was not in. Both numerator and
denominator are now projected, via a new `_without_trivial_modes` that does not normalise
-- normalising a residual discards the magnitude the test is reading.

`response_persistence` restarted from the original seed at every amplitude, so
`rotation_from_previous_deg` was the distance between independent solutions rather than a
continuation measure. It now warm-starts. The trap: `found["direction"]` is mass-weighted
and `initial` is Cartesian, so the seed needs `previous / sqrt(m)` -- invisible at equal
masses, and seeds an entirely wrong direction as soon as they differ. That is the third
defect of this same family in the project's history.

### Discovery is not certification

The rigid-body floor gate is a **necessary condition only**. An error of the form
`P_int A P_int` has exactly zero residual on every rigid-body mode while making internal
eigenvalues arbitrarily wrong, so passing the gate is evidence about invariance-breaking
error -- a wrong unit factor, a transposed block, asymmetric differencing, wrong masses --
and about nothing else. `diagnostics["floor_gate"]` now says
`necessary_condition_only`, and a test constructs exactly that error: the floor passes and
the new internal gate catches it.

`internal_mode_recheck` is that gate. Every negative eigenvalue plus the first
non-negative -- the boundary of the count -- has its Rayleigh quotient re-measured by a
direct Hessian-vector product along the eigenvector, which is different arithmetic from
the column-by-column assembly at the same bandwidth. Two design corrections were needed
and both came from measurement:

**Truncation is extrapolated away, not tolerated.** Comparing an analytic eigenvalue
against a single finite difference charges the matrix for the probe's `O(h^2)` error.
Measured on this project's MACE surface, a rotor mode whose analytic `lambda` is `-0.00318`
differences to `-0.00498` at the default step -- a 56% overestimate, and near the tolerance
that is enough to move the classification and fail a correct Hessian closed. The quotient
is now measured on a ladder `(h/2, h, 2h)` and fitted to the `h^2` law; on a well with an
analytic Hessian the single-step error is `3.4e-5` and what survives extrapolation is
`~1e-14`.

**A moving sign is reported, never gated.** On the bumpy-product test surface, whose ridge
is a gaussian narrower than the probe, the deciding curvature runs
`-33.5, -30.5, -19.4, -2.8, +16.8, +83.6` across steps `0.0025` to `0.04`. That is the
surface talking, and gating on it would reject a correct Hessian for telling the truth.

**What is deliberately not claimed:** that lying on the `h^2` line means the movement is
benign truncation. It was tempting, and wrong. For `V = -A exp(-(q/w)^2)` the symmetric
difference is `(2A/w^2)[1 - h^2/2w^2 + h^4/6w^4 - ...]`, so on `w = 0.06` the
least-squares residual is **0.3% of the spread over a fourfold ladder and 0.5% over a
twenty-fivefold one**, against a 10% criterion. Widening the ladder does not rescue it.
Reading a fitted line as sufficient would be the same mistake as reading a passing floor
gate as a certified Hessian, one level down -- the same error twice in one round, which is
why both are written down. `fits_h2_law` is used only for what it can decide (whether the
extrapolate means anything) and the only claim made about the surface is whether the
classification held.

The analytic-Hessian cross-check probed one fixed structured vector, `direction[::2] = 1`,
and a structured probe can be blind to a structured layout error. It now uses three
`config.seed`-seeded random directions, which miss an error only on a measure-zero set, at
two evaluations each. The raw asymmetry `||H - H^T||/||H||` is recorded **before**
symmetrising, because `0.5*(H + H^T)` destroys exactly the evidence that says the layout
or the units are wrong.

### The probe map

`response_curvature` declares `retraction: "cartesian_chord"` and reports
`bond_strain_by_amplitude` -- the largest fractional bond-length change at each amplitude,
pure geometry and no evaluation. `kappa` is unaffected (`c_v = 0` regardless), but
`anharmonicity` and `sign_stable` are read off finite-amplitude samples, and a sample whose
bonds are stretched is not reporting on the motion it is named after.

**Measured and closed** (`docs/experiments/c1_bond_strain_probe.py`, on the six saddle
geometries the `ts` job left on disk, MACE-OFF24):

    worst bond strain at the top of the band (a = 0.1 A*sqrt(amu)): 0.348%
    scaling: strain ~ a^1.00
    extrapolated to 5% strain: a = 1.44, i.e. 14x the band top

So the chord is not merely the correct map here, its cost is negligible -- and that is now
a number rather than an assumption. The case that motivated the arc machinery in the first
place was an O-H bond going `0.957 -> 2.346 A`, but that was `internal.displace` delivering
a full 2.09 rad torsion, two to three orders of magnitude away from this band.

The cancellation of the first-order term requires the `+a` and `-a` endpoints to carry the
same second-order term. A chord satisfies that trivially; any curved map added later must
be checked, and a test does it: injecting a 5% asymmetry makes `odd - even` a pure `1/a`,
measured as `2.13, 2.06, 2.03, 2.01` as the amplitude halves.

### Configuration

Ten integers added for the min-mode walk, the response descent and the saddle search were
in no validation bucket, and `min_mode_refresh = 0` reached `step % config.min_mode_refresh`
as a plain `ZeroDivisionError` deep inside a run. Fixing the ten is not the fix: the bucket
lists are module-level constants and a meta-test walks `dataclasses.fields`, requiring every
scalar `int`/`float` field to appear in exactly one bucket or in a named exemption set. The
meta-test was verified against negative controls -- removing a name from a bucket, and
adding a stale name -- so it is not idling.

### Internal-coordinate indices

`tuple(int(i) for i in indices)` accepted `-1`, which is a valid numpy index onto the last
atom, and `1.9`, which truncates to `1`. Both then perturb an atom nobody named and return
a real number for the wrong triple. `internal.atom_index` validates instead of coercing,
takes an optional atom count for the upper bound, and is reused by `coordinate`,
`gradient`, `displace`, `named_values`, `rotate_fragment` and the config's
`tracked_coordinates`. numpy integers pass -- indices routinely arrive from `np.where` on a
bond graph -- so the test is on the value being integral, not on `type(value) is int`.

### Stereochemistry tolerances

Both parity discriminants were scale-dependent. The tetrahedral one took
`int(np.sign(det))` of a number whose scale was the product of three bond lengths, so no
fixed threshold could serve a C-C skeleton and a C-H one at once; at a near-planar centre
that number is float noise, so thermal jiggle alone flipped the recorded configuration.
`np.sign(0.0)` is `0` as well, which was silently emitted as a third numeric parity no
configuration corresponds to. The E/Z one compared an unnormalised dot product against
`1e-6` -- a **length in Angstrom standing in for an angle**, which fires on a short
substituent nowhere near perpendicular and never on a long one that is.

Both are now dimensionless in `[-1, 1]`, and both normalisations **preserve the sign**, so
all 332 frozen files and 2136 frames are byte-identical. Below `parity_tolerance` (default
`0.01`) the configuration is recorded as `"planar"`. The tolerance is calibrated, not
picked: near the crossing `d|cos|/dangle` is about `0.017` per degree, so `0.01` is roughly
`0.6` degrees of angular slack, the scale a quench on a soft torsion leaves.

**The consequence, stated because it is real:** the band is a third macrostate. Refusing to
report a noise-determined sign requires an unresolved state, which is the same choice
`stereo_unresolved` already makes. On the critical ladder the flip points moved from
noise-determined (`|cos| = 2.18e-4`) to a deterministic threshold, at the cost of two extra
transitions. A converged minimum never lands there: 90 degrees about a locked double bond
is the torsional barrier top, not a basin.

### The `torch.load` wrapper

Two needs wrapped `torch.load` independently, and the combination was the defect, not
either one. The device-move wrapper saved whatever was installed on entry and restored it
**unconditionally** on exit, so anything installed while it was open vanished when it
closed, silently. Order decided whether a guard survived.

There is now one wrapper (`prrs/torch_guard.py`), installed at most once, never restored,
with the device move a depth counter on it rather than a second wrapper. Nothing has an
unwind path, which makes the failure mode impossible rather than unlikely. The `weights_only`
handling moved in with it, scoped to declared checkpoints rather than the process, so it
also covers `openmm_backend`.

Two coexisting guards then surfaced a second, symmetric defect measured by the other
session: each held its own trust list, the outer one imposed the restrictive default, the
inner stood down as designed, and a checkpoint the inner would have allowed failed to
unpickle. The rule is now that an outer guard overrides only where it has a positive
reason to, never to impose the restrictive default, and the trust registry is shared
through `_prrs_trusted` so a grant on either side is seen by both.

`outermost()` reports an ordering fact and the name says so: `False` has two causes --
removed from the chain, or wrapped by something else and still working -- and telling them
apart cheaply is not possible.

**Known gap:** the guard activates only where it is called, which is
`create_mace_context`. Scripts that import `mace.calculators` directly still emit the
warnings, and the `e3nn` one fires at import time so no guard installed later can help it.
Confirmed by running `sn2_probe` (three warnings, unchanged).

### Claims that were reaching too far

`search.py` now states two limits in its module docstring. The controller is a **sampler,
not an enumerator**: the proposal basis is deliberately not a basis for `3N`, so a reaction
whose coordinate it does not span is unreachable at any budget, and nothing here estimates
coverage. And the two analysis lines are **not independent**: the response probe is seeded
with the eigenvector the certification produced, so a wrong Hessian makes both wrong
together and their agreement is not corroboration.

`follow_min_mode` is labelled for what it is -- a discovery mechanism built out of
certification machinery. Refreshing a projected Hessian while climbing, far from any
stationary point, is a legitimate way to choose a direction; only `confirm_minimum` on the
converged geometry may claim an index.

That second limit is now addressable rather than only stated. `response_descent_enabled`
(default off) runs `lowest_response_direction` from a `config.seed` random direction after
a first-order saddle is confirmed and reports
`response_descent.overlap_with_unstable_mode`. Where a descent that was never told the
answer lands on the eigenvector, that agreement is corroboration. Off by default because
it costs up to `2 * response_descent_max_steps` evaluations per saddle and every recorded
result so far was produced without it.

---

## Verification

    python -m pytest tests/ -q                            205 passed, 1 skipped
    python docs/experiments/digest_regression.py          332 files / 2136 frames unchanged
    python docs/experiments/landscape_regression.py       300/300 bit-identical, EXIT 0
    python docs/experiments/parity_critical_ladder.py     flip points now deterministic
    python docs/experiments/c1_bond_strain_probe.py       0.348% at the band top
    python docs/experiments/p2_relay_overlap_margin.py    350 robust / 19 undetermined

Real-surface coverage: `internal_mode_recheck` ran through `descend_saddle` in the `ts` job
(six first-order saddles, `wrong_index` rate 50% matching the prior baseline) and through
`sn2_probe` on two independently trained MACE models, `confirm_minimum = True order = 0`
unchanged in both. `landscape` exercises `curvature_spectrum` and `floor_residual` on 300
conformers with every verdict field and every numeric field bit-identical.

Not covered on a real surface: `lowest_response_direction` and `response_persistence`
outside the new opt-in path, which nothing else in the repository calls. Unit tests are the
proportionate coverage there.

### P2's frozen relay signs

`_projected_direction_from` fed `align_mode_with_travel`, which decides a relay's sign, so
every recorded relay in P2 was reviewed. This needed no potential evaluation: the two
projections differ only by the mass-weighted centre-of-mass displacement, and the
trajectories are on disk frame by frame.

    projection error over 49982 travel vectors from 348 descent paths
      p0 0.0043   p5 0.0291   p50 0.1029   p95 0.2436   p100 0.4420

The centre of mass moves about 7% of the travel during a quench, so the artefact is not
small. Against the 369 recorded margins the result is bimodal with nothing in between:

    350 relays (94.9%)  |overlap| > 0.9359   above the largest error   robust
     19 relays (5.1%)   |overlap| < 0.0271   below the 5th percentile  undetermined
      0 relays          anywhere in between

The whole error range lies inside that empty gap, so the split does not depend on which
percentile is used; the other session confirmed it independently by finding the largest
sorted-neighbour jump, which falls exactly there with the second largest 61 times smaller.
**Undetermined means the recorded sign was decided by a quantity smaller than the artefact,
not that the sign is wrong.** Settling it needs the relay eigenvector, which is not stored,
so that step needs the potential. P2 is frozen at `0058158d`; whether to recompute those 19
is a revision-process decision, and `p2_round7` holds three of them including the smallest
margin in the repository (`1.02e-08`).

---

## Open

**`chemical_key`'s graph hash is 1-WL colour refinement, not a canonical form.** Collisions
are constructible in principle. Not fixed: strengthening it changes every digest, and
P1/P2 are frozen on the current ones, so it is a revision-process question rather than a
technical one. For C/H/O/halogen with element-typed initial colours, a collision needs a
fairly contrived graph.

**The digest itself is arguable.** Two structures are the same macrostate exactly when
their components agree, and the components are small -- a few fragments and a handful of
parities. Hashing them adds a lossy indirection: a mismatch says only "different", while
comparing components says which field. This round hit the cost directly:
`runs/e2e_ethanol` mismatched on 3 frames while **all five chemical invariants were
identical field for field**, because the component set had changed in an earlier release
(`automorphisms` removed, `stereo_unresolved` added). The digest conflated "the chemistry
changed" with "the record format changed", and only the first is a chemical fact. A short
display id for human reference plus structural comparison for identity would remove that,
at the cost of touching every frozen record's format.

**`torch_guard` covers only the OpenMM path**, as above.

**Two `implementation_sha256` values landed in one rerun batch** because source was edited
while the `ts` job was starting. `scripts/_provenance.py` refuses mixed batches, which is
why it exists; the practical rule is that source freezes before a batch starts, not during.


---

# Round two: a second review pass, and P3 assembled

Four further defects, one of which exposed an error introduced by round one's own fix.

## The recheck gate was looking at a rigid-body zero

`internal_mode_recheck` checked every negative eigenvalue plus `len(negative)`, described
in its own docstring as "the first non-negative -- the boundary of the count".
`curvature_spectrum` returns the full 3N spectrum of the projected matrix sorted ascending,
so the rigid-body modes sit as exact zeros **between** the negative internal modes and the
positive ones. `len(negative)` therefore indexes the first rigid-body mode, every time.

Measured on a C4 chain, `|P_rigid v|` being the eigenvector's projection onto the
rigid-body subspace:

```
minimum (phi=0)          negative=[]   checked=[0]
  [ 0] lam=-1.06e-15  |P_rigid v|=1.0000  <== checked
  [ 6] lam=+7.57e-02  |P_rigid v|=0.0000   <- the real first internal non-negative mode
saddle (phi=60)          negative=[0]  checked=[0,1]
  [ 0] lam=-4.78e-02  |P_rigid v|=0.0000  <== checked   (the real negative mode)
  [ 1] lam=-4.97e-16  |P_rigid v|=1.0000  <== checked   (a rigid-body zero)
  [ 7] lam=+3.34e-01  |P_rigid v|=0.0000   <- never checked
```

A curvature measured along a translation is zero, it agrees with the zero eigenvalue, and
it decides nothing. On a **minimum** -- which has no negative modes -- that was the only
mode the gate looked at, so the gate was spending six force evaluations to confirm nothing.
It could catch a reported negative mode the potential does not support; it could not catch
a missed one, and a missed one is the direction a minimum verdict depends on.

The test that was supposed to demonstrate the internal subspace did not touch it either.
`test_the_rigid_body_floor_passes_an_error_the_internal_recheck_catches` corrupted
`argmin|lambda|` on `_pair(1.2)`, which is index 1:

```
  [1] lam=+0.00e+00  |P_rigid v|=1.0000   <- what the test corrupted
  [5] lam=+1.48e+00  |P_rigid v|=0.0000   <- the only internal mode, never touched
```

It caught the corruption by measuring zero along a rigid-body direction. So the whole
`P_internal A P_internal` class had no coverage.

Fixed by choosing the boundary from the internal subspace: every mode the count read as
negative, whatever subspace it lies in, plus the softest mode with `lambda >= -tol` and
`|P_rigid v| < 0.5`. Per-mode `rigid_body_overlap` and `subspace` are recorded, because the
split is only exact where an internal mode is not degenerate with the rigid-body zeros --
eigh can return a genuinely flat one mixed with them and no threshold fixes that. Cost is
unchanged: one mode checked on a minimum before and after, two on an index-1 saddle.

## Comparing a truncated number against an untruncated one, in the other direction

The subspace fix turned three tests red. The cause was round one's own change.

Round one made the gate compare the h -> 0 extrapolate, because an analytic Hessian carries
no truncation and comparing it against a single difference charges the matrix for the
probe's error -- the measured case being a rotor whose analytic lambda is -0.00318
differencing to -0.00498 at the default step, 56%. That fix is right **for an analytic
matrix**. It is wrong for a finite-difference one, and the newly-reachable boundary mode is
where that shows:

```
bumpy-product surface, boundary mode
  assembled lambda      +4.99e-3
  ladder (h/2, h, 2h)    1.32e-3, 5.29e-3, 2.11e-2     ratios 4.00, 4.00
  h -> 0 extrapolate     5.9e-7
```

The ladder is pure `h^2` with a zero intercept: the whole of the assembled `+5e-3` **is**
the truncation of the differences it was built from. The extrapolate says flat, the matrix
says positive, and the gate failed a Hessian that is exactly what it claims to be.

So the reference has to carry the same truncation as the probe, and which comparison does
that depends on the assembly. `provenance` is now a parameter of
`internal_mode_recheck` rather than something the caller keeps to itself:

| source | reference | why |
|---|---|---|
| `analytic` | the `h -> 0` extrapolate | the matrix has no truncation to match |
| `fd` | the ladder rung at `provenance["step_A"]` | the matrix is that same difference; different arithmetic, one bandwidth |

The ladder is now centred on `provenance["step_A"]` when there is one, so the middle rung
matches the matrix's bandwidth by construction rather than by both reading the same config
field. No threshold was loosened.

This is the same error shape as round one's `fits_h2_law` mistake -- and as the D1 gate
before it -- three times in the same function: two quantities compared as though they were
the same quantity.

## Three smaller ones

**`Probe.__post_init__` truncated instead of validating.** `tuple(int(i) for i in ...)`
turns `1.9` into `1` and `-1` into a numpy index onto the last atom, so a probe acts on an
atom nobody named. `Probe("stretch", (0, 1.9), 1, 0.1, 0)` returned `indices=(0, 1)` with
no error. Now `internal.atom_index`, the same validator round one added for internal
coordinates; the upper bound is still uncheckable there because a Probe is built before it
meets a molecule.

**`rotation_from_seed_deg` was a bit-identical copy of `rotation_from_previous_deg`.**
`response_persistence` read it off `lowest_response_direction`'s `rotation_from_initial_deg`,
but after warm-starting that function's "initial" IS the previous answer:

```
amplitude   from_previous              from_seed
0.01000     None                       48.3947712917
0.02154     0.06443791094514058         0.0644379109
0.04642     0.023768806718724687        0.0237688067
0.10000     0.10912417032551684         0.1091241703
```

Only the first entry was ever the seed. Cumulative drift was recorded nowhere --
`total_rotation_deg` sums per-step rotations, which bounds the drift from above and equals
it only if the direction never turns back. Now measured against the projected seed held
across the loop. `response_persistence` has no caller outside the tests, so this is record
correctness on a function the pipeline does not currently run.

**The final spectrum in `follow_min_mode` sat outside the walk's try.** It builds a
`GuardedCalculator` and displaces along every axis, so the gate can refuse there exactly as
during the walk, and nothing above catches `GateRejected` -- `search.py:641` calls
`follow_min_mode` bare and `seek_saddle` is called bare in turn. Ten of this function's
eleven returns are failures carrying a `reason`; this one took the whole run down. Now a
return with `stage="final_spectrum"`.

## The gate on a real surface

`ts`, 8 shards, MACE-OFF24_medium. Sixteen seeds, verdicts compared field for field against
the frozen baseline in `runs/ts_verdict_baseline_6ede0854.json`:

```
16/16 seeds,  verdict field diffs: 0,  first-order saddles: 6 (baseline 6)
six sign=-1 descents all reached B, saddle_order=0
```

Every recorded gate reading (six, one per reached descent):

```
[gate sign-1] agree=True source=analytic boundary=6 lam=0.002096 ref=0.002096 against=extrapolate h2=True
[gate sign-1] agree=True source=analytic boundary=6 lam=0.001969 ref=0.001969 against=extrapolate h2=True
[gate sign-1] agree=True source=analytic boundary=6 lam=0.001771 ref=0.001771 against=extrapolate h2=True
```

`boundary=6` is an internal mode, so the gate is measuring something for the first time.
The real runs are all `source: analytic` (2372 of 2372 recorded provenance entries), which
is the extrapolate branch. `lambda` runs 0.00177 to 0.00210 -- 1.8 to 2.1 times the
`1e-3` tolerance, which is the soft-mode-near-the-boundary case that could have fired the
gate -- and the extrapolated direct product reproduces the analytic eigenvalue to six
decimals.

## Descent failures had their evidence stripped twice

`follow_unstable_mode` stored the quench diagnostics on the endpoint when it succeeded and
returned `{"reason": ...}` when it failed; `descend_saddle` then flattened them again. So
`soft_mode_moved_on_the_last_round` -- the reason for every one of P3's six failed
descents -- arrived with no coordinate, no offset and no round history. `polish_soft_modes`
set `report["moved"] = True` as a bare boolean, which made the reason unattributable by
construction. Both paths now pass the diagnostics through, and `movements` records the
torsion, the offset and whether it came from a Newton step or a scan.

Nothing is concluded from this yet: the record was added, and producing it needs a run.

---

# P3, assembled

`all shards present`, first time. Five jobs, MACE-OFF24_medium unless stated.

```
=== TS / barrier (8/8) ===
16 seeds: 6 first-order, 10 no-saddle (wrong_index 8, ceiling_blocked 2)
barrier vs A   315.9 .. 316.8 meV   (spread 0.9, n=6)
imaginary      -3086.9 .. -3085.1 cm^-1
q_PT -0.014    r_OO 2.323 A
both sides reached 0/6    connects two nodes 0/6    endpoints are A and B 0/6

=== conformer landscape (8/8) ===
300 relaxed, 229 Hessian-confirmed minima (71 had n_minus > 0)
closed chelate 80 / open form 149
chelate stabilisation 543.6 meV

=== 2D surface (8/8 rigid, 8/8 relaxed) ===
rigid    2304/2304 points, 0 NaN, E range 4541 meV, 2304 force calls
relaxed  2304/2304 points, 0 NaN, E range 1253 meV, 141217 force calls
constraint error  |dq| 0.0000 A   |dr| 0.0000 A

=== three models (3/3) ===
                      dE(B-A) meV   barrier meV      imaginary cm^-1     r_OO A
MACE-OFF24_medium        -13.117     316.3         -3086.2              2.5289
MACE-omol-0-xl-4M        -25.601      84.8    -925.7 .. -921.3          2.5195
MACE-POLAR-1-M           -19.983  95.8 .. 97.3 -1436.5 .. -1436.1       2.5118
all models agree on the (A, B) key pair: True

=== path profile ===
0 frames
```

## The conformer count was counting the sweep against itself

The landscape reported `distinct at 0.81 meV: 1  <- the defensible count` and concluded
`conformer_max_per_node must be >= 1`. All 80 confirmed closed-chelate minima span
0.0032 meV at `r_OO = 2.5236`: the 300 RDKit conformers collapse onto one structure.

But the count ran over the sweep's own rows, and the sweep does not contain A or B -- the
two structures every other P3 job is measured against. Both are closed-chelate minima of
the same substance under the same protocol, and both are separated from the sweep's lowest
by more than the 0.81 meV scale the count uses:

```
sweep lowest closed   -14637.363923 eV   r_OO 2.5236 A
endpoint A            -14637.365178 eV   r_OO 2.5287 A    -1.255 meV
endpoint B            -14637.378327 eV   r_OO 2.5181 A   -14.405 meV
distinct at 0.81 meV, sweep + [A, B]:  3
```

So the defensible count is 3, not 1, and B -- the lowest closed-chelate structure of this
substance by 14 meV -- is not among the 300 conformers the sweep generated. `p3_merge.py`
now pools the recorded endpoints into the count and prints both numbers, keeping the
sweep-only one because it is what says the 300 collapse onto a single structure. Energies
and `r_OO` are read out of the extxyz headers; nothing was recomputed.

## What P3 establishes, and what it does not

Established:

- A and B are distinct chemical states, and three independently trained models agree on the
  pair (`8710064d...` / `b1a4338a...`)
- a first-order saddle exists on the transfer coordinate, barrier 316 meV, imaginary
  wavenumber -3086 cm^-1, barrier spread 0.9 meV over six seeds
- the barrier's magnitude is model-dependent by a factor of 3.7 (316 / 85 / 96 meV); the
  sign of `dE(B-A)` is not
- 300 conformers relaxed, 229 Hessian-confirmed minima, chelate stabilisation 543.6 meV,
  three distinct closed-chelate levels once the endpoints are counted
- both surfaces fully filled, 4608 points, zero constraint error

Not established:

- **that this saddle connects A and B.** `both sides reached 0/6`. All six `sign=+1`
  descents fail with `soft_mode_moved_on_the_last_round`; all six `sign=-1` descents reach
  B. So 316 meV is the height of a first-order saddle reached by climbing from A, not the
  barrier of A -> B. This is P3's one structural gap and it was visible in the `ts` output
  before this round began.
- why those descents fail. The record now carries the polish rounds that would say, but no
  run has produced it. A cause was proposed during this round and is **not** recorded here,
  because it was disputed and never verified.
- POLAR-1's per-atom charges: `Q(shared H)` is `None` in all five recorded entries, in the
  section that exists to record them. The dipole comes through; the charge does not.
- the path profile is empty (`PATHFILE` unset), so no path is drawn on the surface.
