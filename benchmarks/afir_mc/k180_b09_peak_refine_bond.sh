#!/bin/bash
# 2026-11-25 peak-window refinement with bond-change windows (user "可以 顺着往下做").
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
O=benchmarks/afir_mc/results/2026-11-25/b09_peak_refine_bond
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $O
$P -u benchmarks/afir_mc/b09_peak_refine.py benchmarks/afir_mc/results/2026-11-24/e2_b09_sB/b09_sB $O \
  --labels "n25|2-15,n28|2-15,n3|0-2" --ks 0,10,200 --window bond --tag _bond > $O/refine.log 2>&1 || echo "refine FAILED"
echo "all done $(date '+%F %T')"
