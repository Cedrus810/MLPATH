#!/usr/bin/env bash
# O5 optimizer comparison (docs/PLAN_HF_MIDPOINT_O5_2026-11-02.md). One process, run
# slowly beside other work (user 2026-11-02: "慢慢跑"). Rerunning the same command resumes.
#   bash midpoint_optimizer_compare_run.sh OUT_SUBDIR LIMIT REPS_N
set -euo pipefail
cd /home/ruigengji/MLPATH
R=benchmarks/explore_d2/results
OUT=$R/2026-11-02/midpoint_optimizer_compare/$1
mkdir -p "$OUT"
RUNS=(
  $R/2026-10-26/b05_hf_sA $R/2026-10-26/b05_hf_sB
  $R/2026-10-28/b05_hf_v11_sA $R/2026-10-28/b05_hf_v11_sB
  $R/2026-10-28/b05_hf_v11r_sA $R/2026-10-28/b05_hf_v11r_sB
  $R/2026-10-29/b05_hf_v4_sA $R/2026-10-29/b05_hf_v4_sB
)
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
exec /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_d2/midpoint_optimizer_compare.py \
  "$OUT" NFS "$3" 30 20261101 "$2" "${RUNS[@]}"
