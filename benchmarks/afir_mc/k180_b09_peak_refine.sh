#!/bin/bash
# 2026-11-24 discriminating test: peak-window refinement + climb on the b09 sB C:-Cb chains.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
O=benchmarks/afir_mc/results/2026-11-24/b09_peak_refine
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $O
$P -u benchmarks/afir_mc/b09_peak_refine.py benchmarks/afir_mc/results/2026-11-24/e2_b09_sB/b09_sB $O > $O/refine.log 2>&1 || echo "refine FAILED"
echo "all done $(date '+%F %T')"
