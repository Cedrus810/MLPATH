#!/usr/bin/env bash
# Launch one group of midpoint_sign_resolution.py on the full inputs
# (docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md §2). Rerunning resumes.
#   bash midpoint_sign_resolution_run.sh GROUP    (GROUP = S | F)
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/explore_d2/results
OUT=$R/2026-11-02/midpoint_sign_resolution
mkdir -p "$OUT"
RUNS=(
  $R/2026-10-26/b05_hf_sA $R/2026-10-26/b05_hf_sB
  $R/2026-10-28/b05_hf_v11_sA $R/2026-10-28/b05_hf_v11_sB
  $R/2026-10-28/b05_hf_v11r_sA $R/2026-10-28/b05_hf_v11r_sB
  $R/2026-10-29/b05_hf_v4_sA $R/2026-10-29/b05_hf_v4_sB
)
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
exec /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_d2/midpoint_sign_resolution.py \
  "$OUT" "$1" 8 8 20261102 0 "${RUNS[@]}"
