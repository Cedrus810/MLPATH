# P3 heavy jobs — node bundle

Everything GPU-bound for P3 Track A (benzoylacetone chelate proton transfer), packaged to
run on the node. Nothing here is a search run and nothing here evaluates acceptance
criteria: **the P3 criteria are not frozen yet.** These produce the numbers the criteria
need to be written from, plus the visualisation the user asked for.

## TF32 — read this first

**TF32 does not apply to this workload and must not be enabled.** TF32 accelerates fp32
matmuls. The whole PRRS stack runs `default_dtype="float64"`, which TF32 leaves untouched,
so switching it on buys exactly nothing. Dropping to fp32 in order to reach it would be
actively destructive: the frozen results carry rigid-mode floor residuals near `1e-16`, the
analytic Hessian is relied on to machine precision, and the truncation analysis in
`PRRS_STATUS.md §2.6f` is written for float64. An fp32 run would look completely normal and
be worthless.

`_env_guard.py` forces `allow_tf32 = False` on both matmul and cudnn (note torch ships
`cudnn.allow_tf32 = True` by **default**, so this is a real change), sets
`float32_matmul_precision("highest")`, asserts a float64 CUDA tensor, and **verifies the
model file hash** against the frozen `e5ccf5837f685899`, failing closed. If TF32 was on at
import it says so once, because whether it was on decides whether a run is comparable.

The node's actual win is throughput: **23.8 ms/point vs 52.2 ms** on the 2080 Ti, ~2.2×.

## The `torch.load` warning, fixed at the cause

The noise was:

```
mace/calculators/mace.py:226: UserWarning: Environment variable
TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD detected, since the `weights_only` argument was not
explicitly passed to `torch.load`, forcing weights_only=False.
```

It is not filtered, it is removed. torch 2.6 flipped `weights_only` to `True`, which breaks
two loads: e3nn's `constants.pt` holds a pickled `slice` (so `import e3nn.o3` raises
`UnpicklingError` outright — the var is a hard dependency, not a nuisance), and
`mace/calculators/mace.py:226` loads a full pickled model object. The common workaround,
`TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1`, is exactly what produces the warning:
`torch/serialization.py:1499` warns on **every** `torch.load` where the variable is set and
the callsite did not pass the argument itself. At 8 workers that is twice per worker start.

The guard instead:

- registers `slice` via `torch.serialization.add_safe_globals`, so e3nn's file loads under
  `weights_only=True` for real — it genuinely is a weights-only load;
- wraps `torch.load` to pass `weights_only=False` **explicitly**, and only for the one file
  whose sha256 it just verified. Because the callsite now sets the argument, torch's warning
  condition (`weights_only_not_set`) is false and nothing is emitted;
- **removes the variable from `os.environ`** itself. This last part matters: torch reads it
  with `os.getenv` *inside* `torch.load`, so a set variable short-circuits ahead of the
  allowlist. An earlier version tested clean only because the test used `env -u`; under the
  pool, where the variable was inherited from the launching shell, all eight workers printed
  the warning again.

Net result is narrower than the environment variable, not merely quieter: arbitrary pickles
keep torch's safe default, and the single exemption is bound to the checkpoint hash.
Verified with `-W error::UserWarning` armed and the variable deliberately set.

## Running

Every job is **sharded across a worker pool**, `PER_GPU=8` as in
`scripts/run_p2_temperature_scan.sh`. Shards are a stride, not a contiguous block, so each
worker gets a mix of easy and hard items and they finish together.

```bash
# from the repo root
PY=/path/to/env/bin/python ./scripts/p3_node_bundle/run_all.sh all
PY=... ./scripts/p3_node_bundle/run_all.sh ts                  # 16 seeds over 8 shards
PY=... ./scripts/p3_node_bundle/run_all.sh pes 64              # 64 rows over 8 shards
PY=... PER_GPU=4 ./scripts/p3_node_bundle/run_all.sh landscape # 300 conformers
PY=... ./scripts/p3_node_bundle/run_all.sh merge               # assemble + global numbers
PY=... ./scripts/p3_node_bundle/run_all.sh status              # running / done / tracebacks
```

`PATHFILE=<extxyz>` adds the second panel to the `pes` job. Logs land in
`runs/p3_node_logs/<tag>_shardNN.out`. Every worker checkpoints its own JSON after each item,
so a kill leaves partial results rather than nothing — unlike the search scans, which
`PRRS_STATUS.md §9` records as not resumable.

**`merge` is not optional.** Each worker sees only its own stride, so the questions that
matter — how many distinct confirmed minima, the seed distribution, the assembled surface —
are not answerable inside a shard, and a worker printing its own count would print a number
that reads like the answer and is not. `p3_merge.py` also **reports missing shards and exits
nonzero**, because a merge over 7 of 8 looks identical to a complete one otherwise.

`all` runs TS first on purpose: cheapest and most decisive. If there is no first-order saddle
on the transfer coordinate, the surface and the landscape are describing the wrong thing.

### Concurrency: what is measured, and the two separate ceilings

Throughput (`scripts/_saturation_probe.py`): the RTX 5080 saturates at **1.7×** from 1→24
processes; the 2080 Ti reaches **3.2×**. **The faster card has the lower ceiling**, which is
what rules out GPU compute as the contended resource. One thread + one pinned core per
process moved 1.70× to 1.72× — i.e. nothing — ruling out CPU-thread contention too.
**Do not enable CUDA MPS: measured at 0.15×, 6.8× slower.** The remaining suspect is CUDA
context switching between processes sharing one card, still unmeasured.

Memory is a *separate* ceiling and on 11 GB it binds first. Measured here: one worker holds
**~1350 MiB** for this 22-atom float64 system, so `PER_GPU=8` sat at **10780 / 11264 MiB —
96% of a 2080 Ti**, and a ninth process OOMed immediately. `run_all.sh` now prints the
per-GPU fit and warns when `PER_GPU` exceeds it (`P3_WORKER_MIB` to override). A 16 GB card
takes 8 comfortably (~67%); past ~11 it will not.

PER_GPU=8 is kept even though aggregate throughput stops rising near 7: the total time is the
same and the early shards land sooner.

## The three jobs

| script | what it produces | shard unit | cost at 8 shards |
| --- | --- | --- | --- |
| `p3_ts_barrier.py` | first-order saddle, barrier, imaginary frequency, two-sided descent endpoints | 16 seeds (8 amplitudes × 2 climb seeds) | 2 seeds/worker |
| `p3_pes2d.py` | rigid + relaxed `(q_PT, r_OO)` surface, plus the surface along a supplied path | grid rows | ~6 rows/worker at 48² |
| `p3_landscape_full.py` | all 300 conformers relaxed **and Hessian-confirmed** | conformers | ~38/worker |
| `p3_merge.py` | the global numbers, and what is missing | — | seconds |

### `p3_ts_barrier.py`

Drives `follow_min_mode` and `descend_saddle` out of `runner.py`, the same entry points the
search uses, so the barrier is comparable with P1's 419.14 meV and P2's +361 meV rather
than being a separate construction.

**A target direction is mandatory here.** Without one, min-mode following climbs the
globally softest mode. On 3-oxobutanal that returned eight rotor saddles at −65…−80 cm⁻¹
while the transfer sat at −3186 (`PRRS_STATUS.md §8.16` item 1). Benzoylacetone is worse:
a phenyl twist on top of the same methyl rotor, and the source structure's softest internal
mode is only +0.0025 eV/Å². The target used is the shared proton along the O···O line.

The seed grid is 16 (8 amplitudes × 2 climb seeds), not 5, for reasons independent of
filling workers: `PRRS_STATUS.md §8.15.4` measured that **the seed decides which saddle the
climb reaches**, so a five-point scan is a five-sample estimate of a distribution rather
than a search; and P2's residual failures were mode-following not converging within 200
steps, a per-seed outcome. Sixteen makes the converged/unconverged split reportable instead
of anecdotal. Both climb seeds use the **same** target direction — only `follow_min_mode`'s
own `seed` differs, so this varies the walk and not the question.

Saddle order, two-sided descent, and the barrier are reported **separately** — the barrier
is only meaningful once the first two hold. `cm⁻¹` is printed only for a confirmed
first-order saddle, per `PRRS_CURVATURE_ROUTES.md §E`. A climb that returns no geometry is
recorded as `no_saddle` **with its reason**, which is not the same thing as the reaction
failing.

### `p3_pes2d.py` — the plotting trap, stated up front

Coordinates are the pair proposed in `HANDOFF.md §5`: `q_PT = r(O_d–H) − r(O_a–H)`, zero at
the symmetric point and signed toward each well, and `r_OO`, the heavy-atom gate.

**The projection is honest about positions and dishonest about energy if misused.** A PRRS
path lives in 3N dimensions. The relaxed surface value at a projected point is the minimum
over the other 3N−2 coordinates, so **the path generally lies above the surface at its own
projection.** Plotting path energy as if it were surface energy is wrong.

The correct figure is a pair: the path drawn on the contours, and beside it the path's own
energy against arc length **with the surface energy along the same projected line
overlaid**. The gap between those two curves *is the diagnostic* — it answers whether these
two coordinates actually opened the reaction up. The script emits both, so pass a path file
as the 4th argument and the two-panel figure needs no re-run:

```bash
$PY -u scripts/p3_node_bundle/p3_pes2d.py 48 "$MODEL" runs/p3_pes2d runs/p3_ts_barrier/ts_amp0.20.extxyz
```

The relaxed pass warm-starts each point from the neighbour already done, which is what
keeps it near ~100 calls/point instead of ~400.

### `p3_landscape_full.py`

The local probe relaxed 18 of 300 conformers; they collapsed onto 4 distinct energies
within 1.5 meV, all one `chemical_key`, with the open form ~543 meV above the chelate. That
is a sample, not a count, and **none of those 18 had a Hessian run on it** — so none is a
confirmed minimum. `conformer_max_per_node` for P3 needs the real number of distinct
confirmed closed-chelate minima. This runs relax + analytic Hessian on all 300.

## Inputs it expects

`runs/p3_A_source.extxyz` and `runs/p3_B_source.extxyz`, written by
`docs/experiments/p3_source_prep.py` (falls back to `p3_{A,B}_chelated.extxyz`). Both are
confirmed minima, `n_minus = 0`.

**The substance matters and is easy to get wrong.** The chelating enol is the **Z** isomer,
`CC(=O)/C=C(\O)c1ccccc1`. The E isomer points the OH away from the carbonyl: r(O···O) =
4.29 Å across *every one* of 300 conformers, and it cannot undergo intramolecular transfer
at all. A probe built on E returned entirely self-consistent-looking numbers — real minima,
distinct keys, differing graph hashes, a clean single-proton event key — for a substance
that cannot do the reaction. Nothing in the identity layer catches that, which is why a
geometric `r(O···O)` condition belongs in the frozen preflight.

## Eigenvector convention

`curvature_spectrum` returns `eigenvalues[order], vectors[:, order]` — eigenvectors are
**columns**, in the **mass-weighted** basis. `vec[i]` is row `i`, a slice across all modes,
not mode `i`. Cartesian displacement needs `dx = v / sqrt(m)`. Getting this wrong produces
a plausible-looking per-atom amplitude distribution that means nothing; it already happened
once in this bundle's development.
