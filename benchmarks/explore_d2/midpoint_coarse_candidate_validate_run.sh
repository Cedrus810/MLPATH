#!/usr/bin/env bash
# v5 vs v6 (coarse-stage candidate) acceptance on frozen brackets
# (docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md). Rerunning the same command resumes.
#   bash midpoint_coarse_candidate_validate_run.sh OUT_SUBDIR LIMIT REPEATS
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/explore_d2/results
OUT=$R/2026-11-07/midpoint_coarse_candidate_validate/$1
mkdir -p "$OUT"
RUNS=(
  $R/2026-10-26/b05_hf_sA $R/2026-10-26/b05_hf_sB
  $R/2026-10-28/b05_hf_v11_sA $R/2026-10-28/b05_hf_v11_sB
  $R/2026-10-28/b05_hf_v11r_sA $R/2026-10-28/b05_hf_v11r_sB
  $R/2026-10-29/b05_hf_v4_sA $R/2026-10-29/b05_hf_v4_sB
)
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
exec /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_d2/midpoint_coarse_candidate_validate.py \
  "$OUT" 0.025 22 20261102 "$2" "$3" "${RUNS[@]}"
