#!/bin/bash
# 2026-11-19 (E1 of docs/TS_LOCATOR_DESIGN_2026-11-19.md): climbing-image start point and
# fixed tangent: node / span2 / hermite, string snapshots k = 0,2,5,10; same 5 cached N1 DA
# paths as n1_direction_study; afir_validate per arm and k.
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/afir_mc/results/2026-11-19
D=$R/n1_tangent_variants
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
mkdir -p $D
$P -u benchmarks/afir_mc/n1_tangent_variants.py --out $D \
  --case "sA:n6|0-12+3-11:$R/n1_string_smoke180/sA/frames_n6_0-12_3-11.npz" \
  --case "sA:n25|0-11+3-12:$R/n1_string_rest181/sA/frames_n25_0-11_3-12.npz" \
  --case "sA:n4|0-11+3-12:$R/n1_string_rest181/sA/frames_n4_0-11_3-12.npz" \
  --case "sB:n25|0-11+3-12:$R/n1_string_rest181/sB/frames_n25_0-11_3-12.npz" \
  --case "sB:n6|0-11+3-12:$R/n1_string_rest181/sB/frames_n6_0-11_3-12.npz" > $D/study.log 2>&1
for arm in node span2 hermite; do
  for k in 0 2 5 10; do
    specs=()
    for t in sA sB; do [ -d $D/${arm}_k$k/$t ] && specs+=("$D/${arm}_k$k/$t:n1:$t"); done
    if [ ${#specs[@]} -gt 0 ]; then
      $P -u benchmarks/afir_mc/afir_validate.py $D/validate_${arm}_k$k.json "${specs[@]}" > $D/validate_${arm}_k$k.log 2>&1 || echo "validate ${arm} k$k FAILED"
    else echo "${arm} k$k: no converged climb"; fi
  done
done
echo "all done $(date '+%F %T')"
