#!/bin/bash
# Smoke of the 2026-11-15 changes (event-key fallback, early discard off, atoms per call 512)
# on fulda: b05 SC, 2 random orientations, budget 40 000.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
mkdir -p $D
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b05 --mode sc --sc-poses 2 --out $D/b05_smoke_fix --seed 20261023 --shell 2.4 --budget 40000 --wave 32
