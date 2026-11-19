#!/bin/bash
# Overnight queue on fulda (RTX 5090), 2026-11-14 (user: "这个下班 跑 所以你是不是 多准备一点计算").
# All code paths were smoked on fulda the same day. Each item resumes (wave files); a failing
# item does not stop the queue. Rerun this script to resume the whole queue.
#  1. SC-AFIR, reactant MINs = input + 32 random orientations, budget 5 x B* (B* = D2-HF-v1.5
#     ledger of the same system / seed; R at B* is read back from the per-wave spend);
#  2. MC-AFIR (force between whole fragments) as the no-pair-grouping control, budget B*;
#  3. validation (order + reference match + distinct TS / events) of everything above.
set -uo pipefail
cd /home/ruigengji/MLPATH
D=benchmarks/afir_mc/results/2026-11-14/overnight
mkdir -p $D
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
R=benchmarks/afir_mc/run_mc_afir_batch.py
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
run() {  # name, then runner args
  local name=$1; shift
  echo "=== $(date '+%F %T') start $name"
  $P -u $R "$@" --out $D/$name > $D/$name.log 2>&1 || echo "=== $name FAILED (exit $?)"
  echo "=== $(date '+%F %T') end $name"
}
SC=(--mode sc --sc-poses 32 --wave 32)
run b04_sc_sA --system b04 --shell 2.0 --seed 20261023 --budget 161275  "${SC[@]}"
run b04_sc_sB --system b04 --shell 2.0 --seed 20261024 --budget 165930  "${SC[@]}"
run b05_sc_sA --system b05 --shell 2.4 --seed 20261023 --budget 2084515 "${SC[@]}"
run b05_sc_sB --system b05 --shell 2.4 --seed 20261024 --budget 2075845 "${SC[@]}"
run n1_sc_sA  --system n1  --shell 2.1 --seed 20261023 --budget 1676490 "${SC[@]}"
run n1_sc_sB  --system n1  --shell 2.1 --seed 20261024 --budget 1761340 "${SC[@]}"
run b05_mc_sA --system b05 --shell 2.4 --seed 20261023 --budget 416903 --mode mc --wave 32
run n1_mc_sA  --system n1  --shell 2.1 --seed 20261023 --budget 335298 --mode mc --wave 32
echo "=== $(date '+%F %T') start validate"
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json \
  $D/b04_sc_sA:b04:sA $D/b04_sc_sB:b04:sB $D/b05_sc_sA:b05:sA $D/b05_sc_sB:b05:sB \
  $D/n1_sc_sA:n1:sA $D/n1_sc_sB:n1:sB $D/b05_mc_sA:b05:sA $D/n1_mc_sA:n1:sA \
  > $D/validate.log 2>&1 || echo "=== validate FAILED (exit $?)"
echo "=== $(date '+%F %T') done"
