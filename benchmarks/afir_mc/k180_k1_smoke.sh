#!/bin/bash
# 2026-11-26 unimolecular code-path smoke (PLAN_FREEZE_AND_HELDOUT §3.0): s01 malonaldehyde
# enol, SC-AFIR on the molecule itself, default locator (window), budget 30 000.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
D=benchmarks/afir_mc/results/2026-11-26/k1_smoke180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $D
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system s01 --mode sc --wave 32 \
  --out $D/s01 --seed 20261126 --shell 2.0 --budget 30000 > $D/s01.log 2>&1 || echo "s01 FAILED"
echo "all done $(date '+%F %T')"
