#!/usr/bin/env bash
# Launch one mode of midpoint_cold_replay.py on the full inputs (docs/PLAN_HF_MIDPOINT_COLD_2026-11-01.md).
# Rerunning the same command resumes: finished JSONs are skipped.
#   bash midpoint_cold_run.sh MODE    (MODE = AB | CD)
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/explore_d2/results
OUT=$R/2026-11-01/midpoint_cold/full
mkdir -p "$OUT"
RUNS=(
  $R/2026-10-26/b05_hf_sA $R/2026-10-26/b05_hf_sB
  $R/2026-10-28/b05_hf_v11_sA $R/2026-10-28/b05_hf_v11_sB
  $R/2026-10-28/b05_hf_v11r_sA $R/2026-10-28/b05_hf_v11r_sB
  $R/2026-10-29/b05_hf_v4_sA $R/2026-10-29/b05_hf_v4_sB
)
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
exec /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_d2/midpoint_cold_replay.py \
  "$OUT" "$1" 3 30 20261101 0 "${RUNS[@]}"
