#!/bin/bash
# 2026-11-25 generality test: bond-change peak windows on the 19 b05 sB k-scan chains (user
# "你要看183 空 你吃 183"). k = 0, 10, 200; then afir_validate per k.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
O=benchmarks/afir_mc/results/2026-11-25/b05_peak_window
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 AFIR_MAX_ATOMS_PER_CALL=256
mkdir -p $O
$P -u benchmarks/afir_mc/peak_window_refine.py benchmarks/afir_mc/results/2026-11-23/e2_b05_kscan/b05_sB/chains.npz b05 $O \
  --labels "n1|2-12+13-14,n1|3-14,n2|2-12+13-14,n2|3-14,n4|2-7+11-14,n5|1-7+3-14,n5|3-10+5-14,n5|3-14,n5|3-14+6-11,n5|3-14+7-8,n6|2-8+7-14,n6|5-14,n7|5-14,n8|0-11+13-14,n8|1-11+13-14,n8|1-2+11-14,n8|11-14,n8|2-12+11-14,n8|3-10+11-14" --ks 0,10,200 > $O/refine.log 2>&1 || echo "refine FAILED"
for k in 0 10 200; do
  $P -u benchmarks/afir_mc/afir_validate.py $O/validate_k${k}.json "$O/k${k}:b05:sB" > $O/validate_k${k}.log 2>&1 || echo "validate k${k} FAILED"
done
echo "all done $(date '+%F %T')"
