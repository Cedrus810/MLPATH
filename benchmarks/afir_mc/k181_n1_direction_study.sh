#!/bin/bash
# 2026-11-19: when is the path direction good enough for the climbing image? 5 cached N1 DA
# paths; string snapshots k = 0,10,25,50,100,200,400; climb from each; afir_validate per k.
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/afir_mc/results/2026-11-19
D=$R/n1_direction_study
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
mkdir -p $D
$P -u benchmarks/afir_mc/n1_direction_study.py --out $D \
  --case "sA:n6|0-12+3-11:$R/n1_string_smoke180/sA/frames_n6_0-12_3-11.npz" \
  --case "sA:n25|0-11+3-12:$R/n1_string_rest181/sA/frames_n25_0-11_3-12.npz" \
  --case "sA:n4|0-11+3-12:$R/n1_string_rest181/sA/frames_n4_0-11_3-12.npz" \
  --case "sB:n25|0-11+3-12:$R/n1_string_rest181/sB/frames_n25_0-11_3-12.npz" \
  --case "sB:n6|0-11+3-12:$R/n1_string_rest181/sB/frames_n6_0-11_3-12.npz" > $D/study.log 2>&1
for k in 0 10 25 50 100 200 400; do
  specs=()
  for t in sA sB; do [ -d $D/k$k/$t ] && specs+=("$D/k$k/$t:n1:$t"); done
  if [ ${#specs[@]} -gt 0 ]; then
    $P -u benchmarks/afir_mc/afir_validate.py $D/validate_k$k.json "${specs[@]}" > $D/validate_k$k.log 2>&1 || echo "validate k$k FAILED"
  else echo "k$k: no converged climb"; fi
done
echo "all done $(date '+%F %T')"
