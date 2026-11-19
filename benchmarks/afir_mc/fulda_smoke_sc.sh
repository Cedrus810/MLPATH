#!/bin/bash
# Smoke of SC-AFIR mode (batched runner) on fulda: b04, reactant MINs from the input + 8 poses.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-14
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b04 --mode sc --sc-poses 8 --out $D/b04_sc_smoke --seed 20261023 --shell 2.0 --budget 60000 --wave 32
