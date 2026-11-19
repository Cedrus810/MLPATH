#!/bin/bash
# 2026-11-15 (user 「做」): bias-force breakdown, rho vs frozen q, at the biased minimum.
# b05: E2-position paths of b05_smoke_expand_rho (2 orientations); N1: the three known
# failure paths of the PBS N1 SC sA run (32 orientations).
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/push_force_breakdown.py --run $D/b05_smoke_expand_rho --system b05 \
  --seed 20261023 --shell 2.4 --sc-poses 2 --labels 'n1|2-14,n0|2-14,n1|5-14' \
  --out $D/push_breakdown/b05
$P -u benchmarks/afir_mc/push_force_breakdown.py --run $D/pbs/n1_sc_sA --system n1 \
  --seed 20261023 --shell 2.1 --sc-poses 32 --labels 'n28|0-12,n24|0-11,n6|0-12' \
  --out $D/push_breakdown/n1
echo "all done $(date '+%F %T')"
