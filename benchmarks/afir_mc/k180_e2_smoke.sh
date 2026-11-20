#!/bin/bash
# 2026-11-19 smoke of e2_string_climb.py on 180: b04 sA (11 paths of PBS 29729), k = 0, 2.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
R=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
D=benchmarks/afir_mc/results/2026-11-19/e2_smoke180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
mkdir -p $D
$P -u benchmarks/afir_mc/e2_string_climb.py --run $R/b04_sA --system b04 --out $D/b04_sA > $D/b04_sA.log 2>&1
for k in 0 2; do
  $P -u benchmarks/afir_mc/afir_validate.py $D/validate_k$k.json "$D/b04_sA/k$k:b04:sA" > $D/validate_k$k.log 2>&1 || echo "validate k$k FAILED"
done
echo "all done $(date '+%F %T')"
