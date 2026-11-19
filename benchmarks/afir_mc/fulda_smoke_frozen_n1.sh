#!/bin/bash
# Smoke of the frozen-weight SC coordinate (design §2.1, user 2026-11-15) on fulda:
# N1 SC, 2 random orientations (nodes n0-n2 include DA-position pairs), budget 80 000.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system n1 --mode sc --sc-poses 2 --sc-coord frozen \
  --out $D/n1_smoke_frozen --seed 20261023 --shell 2.1 --budget 80000 --wave 32
