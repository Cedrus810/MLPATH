#!/bin/bash
# Smoke 2026-11-15: expand_brackets was not reached by fulda_smoke_expand.sh (frozen q gave no
# bond change on b05 / N1). b05 SC with --sc-coord rho (bond changes known), otherwise as
# b05_smoke_scan2000 (2 orientations, budget 300 000) for the bracket-cost comparison.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b05 --mode sc --sc-poses 2 --sc-coord rho \
  --out $D/b05_smoke_expand_rho --seed 20261023 --shell 2.4 --budget 300000 --wave 32
