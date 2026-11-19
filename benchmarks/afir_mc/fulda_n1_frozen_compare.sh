#!/bin/bash
# 2026-11-15 (user: 冻结权重坐标，保留 γ 与 SC 分组，用三条已知失败路径对照): N1 SC sA paths
# n28|0-12, n24|0-11, n6|0-12, each with rho_AB / frozen q, with / without the two-maxima stop.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/n1_frozen_compare.py --run $D/pbs/n1_sc_sA --seed 20261023 --shell 2.1 \
  --sc-poses 32 --labels 'n28|0-12,n24|0-11,n6|0-12' --out $D/n1_frozen_compare
