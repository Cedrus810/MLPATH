#!/bin/bash
# Smoke 2026-11-15 (user: 双对协同探索, alpha (r_kl + r_mn)/sqrt 2, TS pair-set union):
# N1 SC, 2 orientations; single-pair base layer then the two-pair layer, total budget 60 000
# (single layer ~29k in n1_smoke_pair, so a few two-pair waves).
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system n1 --mode sc --sc-poses 2 \
  --out $D/n1_smoke_double --seed 20261023 --shell 2.1 --budget 60000 --wave 32
