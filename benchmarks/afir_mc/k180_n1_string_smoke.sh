#!/bin/bash
# Smoke 2026-11-19 on 180: simplified string + climbing image as TS locator on the same ONE
# N1 DA path as the failed CI-NEB smoke (PBS 29729 sA n6|0-12+3-11; kernel brackets hysteresis).
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
D=benchmarks/afir_mc/results/2026-11-19/n1_string_smoke180
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
mkdir -p $D
$P -u benchmarks/afir_mc/n1_string_check.py --run $R/n1_sA --seed 20261023 --shell 2.1 --sc-poses 32 \
  --labels 'n6|0-12+3-11' --out $D/sA
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json $D/sA:n1:sA
echo "all done $(date '+%F %T')"
