# PARS downhill connectivity counterexample

The analytic external field `E=x²(x²−1)²+(y²+z²)/2` has minima at -1,0,+1
and index-one saddles at ±1/√3. A band between the outer minima cannot use
one of these saddles as proof of direct connectivity: each saddle descends
to the intermediate minimum and only one outer minimum.

```sh
PYTHONPATH=pars/src /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_connectivity/run.py --output runs/new_connectivity_counterexample.json
```

The output must be new. It records the exact potential, effective configuration,
runner/module hashes, numerical environment, band and complete refinement report.
The script reports the actual outcome, without silently replacing it by an expected value.

`results/2026-10-04/baseline-14d8325.json` was generated using actual exported
PARS sources from commit `14d8325`, with this same calculator/script and environment.
It returns `ci_neb_converged` and ~0.148148 eV. `current.json` retains the two
minimum witnesses, identifies the intermediate basin and refuses the barrier.
`mutations.json` records two regressions caught in isolated module copies.

These are engineering checks on an analytic external field, not molecular or
QM evidence, and Cartesian downhill descent is not mass-weighted IRC.
See [implementation, limits and budgets](../../docs/CONNECTIVITY_FIX_2026-10-04.md).

An additional narrow-well counterexample is available in `resolution.py`:
the 0.05 A signed kick skips a 0.002 A wide intermediate well. Actual exported
`7d299d5` sources accept the apparent outer-minimum connection; the finer replay
finds the intermediate and refuses it. Records are in `results/2026-10-05/`.
See [resolution check and costs](../../docs/CONNECTIVITY_RESOLUTION_2026-10-05.md).
Two tested resolutions are numerical evidence, not continuous IRC proof.
