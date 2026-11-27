#!/bin/bash
# 2026-11-26 smoke of the AFIR-W1 candidate rules (PLAN_FREEZE_AND_HELDOUT §2; user "做起来"):
# window locator (default) + A/C/D + --double-xx first + --double-order class + --reflow, on
# N1 sA with 4 orientations, 2 reflow orientations per new conformer, budget 80 000.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
D=benchmarks/afir_mc/results/2026-11-26/w1_smoke180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $D
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system n1 --mode sc --sc-poses 4 --wave 32 \
  --out $D/n1_sA --seed 20261023 --shell 2.1 --budget 80000 \
  --sc-split 1 --skip-start-min 1 --confirm-dedup 1 --double-xx first --double-order class --reflow 2 \
  > $D/n1_sA.log 2>&1 || echo "n1 FAILED"
echo "all done $(date '+%F %T')"
