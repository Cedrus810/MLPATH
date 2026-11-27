#!/bin/bash
# 2026-11-25 smoke of --locator window (explore.afir_window; user "可以") on b04 sA, the
# PBS 29729 settings otherwise, budget capped at 40 000; then afir_validate.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
D=benchmarks/afir_mc/results/2026-11-25/window_smoke180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $D
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b04 --mode sc --sc-poses 32 --wave 32 \
  --out $D/b04_sA --seed 20261023 --shell 2.0 --budget 40000 --locator window > $D/b04_sA.log 2>&1 || echo "b04 FAILED"
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json $D/b04_sA:b04:sA > $D/validate.log 2>&1 || echo "validate FAILED"
echo "all done $(date '+%F %T')"
