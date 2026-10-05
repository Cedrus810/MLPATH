# Two saddles, identical mapped endpoints

This is a CPU analytic engineering benchmark, **not molecular/QM validation**.
The frozen [protocol](protocol.json) defines the potential, pool, seeds, budget,
refinement, stationary-point checks and acceptance before the first execution.

One mobile particle in an external field has
`E=(x²+y²−1)²+0.8y²+0.08y³+z²`. Two harmonically anchored observers make
contact histories depend on geometry. Two exact stationary solutions at
`x=0, y≈−0.8051774/+0.7451774` have one negative Hessian eigenvalue each.
Independent downhill steps from ± the negative mode reach the same two
positive-Hessian minima (`x=−1/+1`). Contacts use the usual distance/covalent
radius rule; the selector receives no saddle identities or refined barriers.
Rigid-motion removal is disabled because the potential is an external field.

Twelve path seeds are generated from fixed geometry amplitudes. Selection
plans are saved **before** the reference/oracle or any refinement is measured.
PARS, channel-wide pruning, cheap energy excursion and 20 random selections
share the pool, refinement parameters, a cap of two candidates and 6000 force
computations per arm. Each arm actually refines its own seeds with a cold
calculator; outcomes are not replayed from the oracle cache. The counter
includes endpoint Hessians, witness checks, NEB/CI and failed calculations.
The all-pool oracle and shared preparation/certification costs are reported
separately. Equal candidate counts alone are not called equal compute budgets.

Run from the repository root in the declared science environment:

```sh
/home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_two_saddles/run.py --output runs/my_two_saddles
```

The output directory must be new. `manifest.json` freezes code, environment,
protocol, candidate input digest and all selection plans; `reference.json`
records Hessians/descent evidence; `oracle.json` covers every candidate;
`arms.json`/`report.json` record cold executions, costs and coverage. Known
reference saddles, including a saddle missed by every selector, remain in the
denominator. A failed or unmatched refinement contributes zero coverage and
still consumes budget. A budget exhaustion is a failure, never a free retry.

A pass requires independently certified reference saddles, both in the pool,
and both retained by PARS inside the budget. It does **not** require speedup
or superiority over random. Current channel-wide pruning can also retain
both on this pool: the benchmark is a working multi-saddle comparison, not
proof that this pool isolates the pruning bug. The separate dominated-tube
regression checks that bug. Contact signatures can still merge distinct
geometric saddles; this test does not establish universal descriptor fidelity.

The two early runs did not reach acceptance: one exposed missing provenance
constructor fields, another showed LBFGS returning to the saddle from a
small negative-mode displacement. Reference descent now uses independently
checked downhill steps. These were runner corrections; potential, pool,
budget and acceptance thresholds in the frozen protocol were not retuned.

Results and remaining molecular gates are recorded in
[the repair record](../../docs/MAJOR_FIXES_2026-10-04.md).
