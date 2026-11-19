#!/bin/bash
# SC-AFIR on b04 (with batched bracket narrowing), seed sA, budget = D2-HF-v1.5 b04 sA ledger (32 255 forwards).
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-14
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b04 --mode sc --sc-poses 8 --out $D/b04_sc_sA_narrow --seed 20261023 --shell 2.0 --budget 32255 --wave 32
