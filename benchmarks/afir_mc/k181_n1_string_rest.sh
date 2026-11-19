#!/bin/bash
# 2026-11-19 on 181 (user: "181 直接跑吧 然后disown"): simplified string + climbing image on the
# other 4 N1 DA paths of PBS 29729 (sA n25|0-11+3-12, n4|0-11+3-12; sB n25|0-11+3-12,
# n6|0-11+3-12), then afir_validate.py.
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
D=benchmarks/afir_mc/results/2026-11-19/n1_string_rest181
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
mkdir -p $D
$P -u benchmarks/afir_mc/n1_string_check.py --run $R/n1_sA --seed 20261023 --shell 2.1 --sc-poses 32 \
  --labels 'n25|0-11+3-12,n4|0-11+3-12' --out $D/sA > $D/sA.log 2>&1 || echo "sA FAILED"
$P -u benchmarks/afir_mc/n1_string_check.py --run $R/n1_sB --seed 20261024 --shell 2.1 --sc-poses 32 \
  --labels 'n25|0-11+3-12,n6|0-11+3-12' --out $D/sB > $D/sB.log 2>&1 || echo "sB FAILED"
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json $D/sA:n1:sA $D/sB:n1:sB > $D/validate.log 2>&1 || echo "validate FAILED"
echo "all done $(date '+%F %T')"
