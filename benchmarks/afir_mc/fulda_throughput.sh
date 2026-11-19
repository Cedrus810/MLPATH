#!/bin/bash
# Throughput of batched AFIR paths vs wave size on fulda (RTX 5090); b05 then N1.
set -euo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-14
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
$P -u benchmarks/afir_mc/throughput.py $D/throughput_b05.json --system b05 --shell 2.4 --seed 20261023
$P -u benchmarks/afir_mc/throughput.py $D/throughput_n1.json --system n1 --shell 2.1 --seed 20261023
