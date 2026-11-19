#!/bin/bash
# Replay 2026-11-15 (user 「一起重放」): three N1 SC sA paths that pushed the DA positions to
# rho below the reference saddles' rho without forming a bond; every frame kept.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/n1_da_replay.py --run $D/pbs/n1_sc_sA --seed 20261023 --shell 2.1 \
  --sc-poses 32 --labels 'n28|0-12,n24|0-11,n6|0-12' --out $D/n1_da_replay
