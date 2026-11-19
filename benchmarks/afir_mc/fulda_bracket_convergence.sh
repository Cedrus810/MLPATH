#!/bin/bash
# Measurement 2026-11-15 (user 「跑」): convergence of the fixed-rho transverse relaxations of the
# bracket scan on b05 SC sA (PBS 28756), 3 bond-changing paths, up to 2000 FIRE steps.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/bracket_convergence.py --run $D/pbs/b05_sc_sA --system b05 \
  --seed 20261023 --shell 2.4 --sc-poses 32 --ks 9,39,93 --max-steps 2000 \
  --out $D/bracket_convergence
