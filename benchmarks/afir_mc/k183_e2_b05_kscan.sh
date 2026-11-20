#!/bin/bash
# 2026-11-20: b05 sB, the 19 paths where the HF kernel matched a reference and E2 (k=2) did
# not, with the same start (|dE0| <= 0.005 eV): string snapshots k = 0,2,10,50,200,1000,
# climbing image from each. Frames copied from e2_pbs/b05_sB/frames (no replay).
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
R=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
D=benchmarks/afir_mc/results/2026-11-20/e2_b05_kscan
L="n1|2-12+13-14,n1|3-14,n2|2-12+13-14,n2|3-14,n4|2-7+11-14,n5|1-7+3-14,n5|3-10+5-14,n5|3-14,n5|3-14+6-11,n5|3-14+7-8,n6|2-8+7-14,n6|5-14,n7|5-14,n8|0-11+13-14,n8|1-11+13-14,n8|1-2+11-14,n8|11-14,n8|2-12+11-14,n8|3-10+11-14"
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
$P -u benchmarks/afir_mc/e2_string_climb.py --run $R/b05_sB --system b05 --out $D/b05_sB \
  --labels "$L" --checkpoints 0,2,10,50,200,1000 > $D/b05_sB.log 2>&1
for k in 0 2 10 50 200 1000; do
  $P -u benchmarks/afir_mc/afir_validate.py $D/validate_k${k}.json "$D/b05_sB/k${k}:b05:sB" > $D/validate_k${k}.log 2>&1 || echo "validate k${k} FAILED"
done
echo "all done $(date '+%F %T')"
