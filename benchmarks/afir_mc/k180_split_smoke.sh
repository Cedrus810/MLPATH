#!/bin/bash
# 2026-11-20 smoke of the trial rules --sc-split / --skip-start-min / --confirm-dedup
# (docs/PLAN_AFIR_PAIR_SPLIT_2026-11-20.md §7.2 A, C, D; user "可以动手试试") on b04 sA:
# same seed / shell / poses as PBS 29729 b04_sA, budget capped at 40 000 (PBS 29729 first
# SN2 hit at 23 518). Not a recall measurement.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
D=benchmarks/afir_mc/results/2026-11-20/split_smoke180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $D
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b04 --mode sc --sc-poses 32 --wave 32 \
  --out $D/b04_sA --seed 20261023 --shell 2.0 --budget 40000 \
  --sc-split 1 --skip-start-min 1 --confirm-dedup 1 > $D/b04_sA.log 2>&1 || echo "b04 FAILED"
echo "all done $(date '+%F %T')"
