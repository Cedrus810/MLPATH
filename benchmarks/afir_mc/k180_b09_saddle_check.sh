#!/bin/bash
# 2026-11-25 post-hoc order + two sides of the b09 saddle from bond-window refinement.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
O=benchmarks/afir_mc/results/2026-11-25/b09_saddle_check
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $O
$P -u benchmarks/afir_mc/b09_saddle_check.py benchmarks/afir_mc/results/2026-11-25/b09_peak_refine_bond/refine_bond.json "n28|2-15" 10 8 11 $O/n28_k10.json > $O/check.log 2>&1 || echo "check FAILED"
echo "all done $(date '+%F %T')"
