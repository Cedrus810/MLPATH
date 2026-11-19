#!/bin/bash
# Smoke 2026-11-15 (user: SC 只推被选中的原子对，推力坐标与 TS 定位坐标解耦): defaults
# --sc-coord pair, --ts-coord frozen-peak, no two-maxima stop, brackets grown from the peak.
# b05 SC and N1 SC, 2 orientations each (as the earlier smokes).
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b05 --mode sc --sc-poses 2 \
  --out $D/b05_smoke_pair --seed 20261023 --shell 2.4 --budget 300000 --wave 32 > $D/b05_smoke_pair.log 2>&1 || echo "b05 FAILED"
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system n1 --mode sc --sc-poses 2 \
  --out $D/n1_smoke_pair --seed 20261023 --shell 2.1 --budget 300000 --wave 32 > $D/n1_smoke_pair.log 2>&1 || echo "n1 FAILED"
echo "all done $(date '+%F %T')"
