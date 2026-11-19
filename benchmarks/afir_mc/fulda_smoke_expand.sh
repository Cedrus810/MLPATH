#!/bin/bash
# Smoke 2026-11-15 (user: 先去掉“两个极大”停止，再做括区逐步展开): defaults now
# --stop-after-maxima 0, --sc-coord frozen, brackets grown from the peak (expand_brackets).
# b05 SC (as b05_smoke_scan2000: 2 orientations, budget 300 000) then N1 SC (as
# n1_smoke_frozen: 2 orientations, budget 80 000).
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b05 --mode sc --sc-poses 2 \
  --out $D/b05_smoke_expand --seed 20261023 --shell 2.4 --budget 300000 --wave 32 > $D/b05_smoke_expand.log 2>&1 || echo "b05 FAILED"
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system n1 --mode sc --sc-poses 2 \
  --out $D/n1_smoke_expand --seed 20261023 --shell 2.1 --budget 80000 --wave 32 > $D/n1_smoke_expand.log 2>&1 || echo "n1 FAILED"
echo "all done $(date '+%F %T')"
