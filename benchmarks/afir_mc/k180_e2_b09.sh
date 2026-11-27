#!/bin/bash
# 2026-11-24: b09 sB (PBS b09_phase1, 2026-11-22) bond-changing paths whose HF kernel gave no
# saddle -> short String + fixed-direction climb (user 2026-11-23 「可以」). Checkpoints
# k = 0,2,10,50,200. Identification afterwards (classify_b09 with symmetric RMSD).
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
R=benchmarks/afir_mc/results/2026-11-22/b09_phase1
D=benchmarks/afir_mc/results/2026-11-24/e2_b09_sB
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $D
$P -u benchmarks/afir_mc/e2_string_climb.py --run $R/b09_sB --system b09 --out $D/b09_sB \
  --labels "n25|2-15,n28|2-15,n3|0-2" --checkpoints 0,2,10,50,200 > $D/b09_sB.log 2>&1 || echo "e2 FAILED"
echo "all done $(date '+%F %T')"
