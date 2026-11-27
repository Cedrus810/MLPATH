#!/bin/bash
# 2026-11-24 b09 chain diagnosis (geometry only, CPU) for the two C:-Cb paths of b09 sB.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
D=benchmarks/afir_mc/results/2026-11-24/e2_b09_sB/b09_sB
O=benchmarks/afir_mc/results/2026-11-24/b09_chain_diag
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
for L in "n25|2-15" "n28|2-15"; do
  $P -u benchmarks/afir_mc/b09_chain_diag.py $D "$L" $O/${L/|/_}.json > $O/${L/|/_}.log 2>&1 || echo "$L FAILED"
done
echo "all done $(date '+%F %T')"
