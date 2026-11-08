#!/usr/bin/env bash
# TODO item 9c: diagnose v5 fine-stage stops (midpoint_fine_diag.py). Rerunning resumes.
#   bash midpoint_fine_diag_run.sh OUT_SUBDIR REPS ITEMS
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/explore_d2/results
OUT=$R/2026-11-06/midpoint_fine_diag/$1
mkdir -p "$OUT"
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
exec /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_d2/midpoint_fine_diag.py \
  "$OUT" 0.025 $R/2026-11-04/midpoint_fire_validate/full/items.json 20261102 "$2" "$3" $R
