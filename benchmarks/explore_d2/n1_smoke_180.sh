#!/usr/bin/env bash
# HF-5 N1 smoke on 180 (user 2026-11-09 「去180 冒烟」): measure-noise, then a small v1.2 run with
# that noise_eV (b04-smoke-sized L0 and budgets). Covers the code path of n1_hf5_chain_pbs.zsh.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
R=benchmarks/explore_d2/results/2026-11-09
N=$R/n1_noise_180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
[ -f $N/noise.json ] || $P -u benchmarks/explore_d2/run_d2.py --system n1 --measure-noise $N \
  --seed 20261023 --shell 2.1 --n-samples 8 --n-far 4
NOISE=$($P -c "import json; print(json.load(open('$N/noise.json'))['noise_eV_rule'])")
echo "noise_eV $NOISE"
$P -u benchmarks/explore_d2/run_d2.py --system n1 --run $R/n1_hf5_v12_smoke --hf --start-cluster \
  --same-graph-ridge --seed 20261023 --shell 2.1 --n-samples 2 --n-far 2 --noise-eV $NOISE \
  --per-pose-quota 100 --search-budget 20000 --certify-budget 20000 --scan-cap 400 \
  --reverse-cap 100 --release-cap 300 --certify-cap 800 --audit-fraction 0.1 --max-active 4 \
  --wave-size 16
echo "=== done $(date '+%F %T') ==="
