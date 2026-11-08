#!/usr/bin/env bash
# Smoke of midpoint_sign_resolution.py (S and F, 1 point x 1 repeat) and
# midpoint_coarse_validate.py (1 failed + 1 control item, provisional bound 0.0105).
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/explore_d2/results
D=$R/2026-11-02/midpoint_coarse_smoke
RUNS=(
  $R/2026-10-26/b05_hf_sA $R/2026-10-26/b05_hf_sB
  $R/2026-10-28/b05_hf_v11_sA $R/2026-10-28/b05_hf_v11_sB
  $R/2026-10-28/b05_hf_v11r_sA $R/2026-10-28/b05_hf_v11r_sB
  $R/2026-10-29/b05_hf_v4_sA $R/2026-10-29/b05_hf_v4_sB
)
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
date
$P benchmarks/explore_d2/midpoint_coarse_validate.py $D/validate 0.0105 1 20261102 1 "${RUNS[@]}"
date
$P benchmarks/explore_d2/midpoint_sign_resolution.py $D/sign S 1 1 20261102 1 "${RUNS[@]}"
$P benchmarks/explore_d2/midpoint_sign_resolution.py $D/sign F 1 1 20261102 1 "${RUNS[@]}"
date
echo SMOKE_DONE
