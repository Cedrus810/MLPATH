#!/bin/bash
# Smoke of SCAN_MAX_STEPS = 2000 (bracket scan / narrowing step cap; DONE_LOG 2026-11-15) on fulda:
# b05 SC, 2 random orientations (as b05_smoke_fix), budget 300 000 so the bracket stage is not cut.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-15
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b05 --mode sc --sc-poses 2 --out $D/b05_smoke_scan2000 --seed 20261023 --shell 2.4 --budget 300000 --wave 32
