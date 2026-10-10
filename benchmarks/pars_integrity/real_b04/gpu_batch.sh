#!/usr/bin/env bash
# Real POLAR-1 check of the 2026-10-05 electronic-state fix on b04 (charge -1, singlet).
#
# 1. counterfactual.py: refine the probe4_ledger representatives the way the pre-fix
#    code did (electronic_state stripped), recording what POLAR-1 actually received.
# 2. For each run: pars select, then pars refine through logged_polar1, which counts
#    the (charge, spin) of every structure POLAR-1 evaluates into <out>/rep_<run>/seen.json.
#
# Resumable: pars refine skips representatives already refined under the same conditions.
# Originally run 2026-10-05 on kasuga01 with outputs in a session scratchpad (since moved to
# results/2026-10-05/real_b04/); paths below are the in-repo equivalents.
#
# Usage: benchmarks/pars_integrity/real_b04/gpu_batch.sh <out dir> [run ...]
#   default runs: b04_probe b04_probe2_openceiling b04_probe3_rigidbond b04_frozen1 b04_frozen2
set -euo pipefail
cd /home/ruigengji/MLPATH
OUT=${1:?usage: gpu_batch.sh <out dir> [run ...]}
shift
RUNS=("$@")
[ ${#RUNS[@]} -gt 0 ] || RUNS=(b04_probe b04_probe2_openceiling b04_probe3_rigidbond b04_frozen1 b04_frozen2)
HERE=benchmarks/pars_integrity/real_b04
export OPENBLAS_NUM_THREADS=1 PYTHONPATH=pars/src:src:.:$HERE
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1   # a hard crash leaves a traceback in the log
PY=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
# TRACE_EVERY=<seconds>: run refine through traced_pars.py (periodic stack dumps into the log)
PARS=($PY -m pars)
[ -z "${TRACE_EVERY:-}" ] || PARS=($PY "$HERE/traced_pars.py")
set -x

if [ ! -f "$OUT/counterfactual.json" ]; then
  [ -d "$OUT/b04_real/sel" ] || $PY -m pars select --run-dir runs/b04_probe4_ledger --out "$OUT/b04_real/sel"
  $PY "$HERE/counterfactual.py" "$OUT/b04_real/sel" "$OUT/counterfactual.json"
fi
for r in "${RUNS[@]}"; do
  [ -d "$OUT/rep_$r/sel" ] || $PY -m pars select --run-dir "runs/$r" --out "$OUT/rep_$r/sel"
  # One run failing must not stop the rest; record its exit code (2026-10-08: the batch ended
  # inside probe3 t4 after ~67 min with no output, and nothing after it ran).
  status=0
  SEEN_OUT="$OUT/rep_$r/seen.json" "${PARS[@]}" refine --selection "$OUT/rep_$r/sel" \
    --calculator-factory logged_polar1:make || status=$?
  echo "RUN_EXIT $r $status $(date '+%F %T')"
done
echo BATCH_DONE
