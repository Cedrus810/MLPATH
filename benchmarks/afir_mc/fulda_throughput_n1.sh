#!/bin/bash
# N1 part of fulda_throughput.sh (the first run stopped at the b05 wave-256 CUDA OOM).
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-14
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/throughput.py $D/throughput_n1.json --system n1 --shell 2.1 --seed 20261023 --waves 8 32 128
