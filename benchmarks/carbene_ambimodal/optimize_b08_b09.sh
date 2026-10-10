#!/usr/bin/env bash
# Re-optimise the b08 / b09 stationary points of Wang & Xue (JACS 2026, 148, 8723) on MACE-POLAR-1-M
# (gas phase) with geomeTRIC, starting from the SI's M06-2X/def2-TZVP-SMD geometries.
#   b08: :CCl2 + 1,2-bis(methylene)cycloheptane (7C) -- TS-CCl2 (ambimodal), TS-CCl2-C4
#   b09: dimethoxycarbene + 1-cyclohexene-1-carboxaldehyde (6O) -- TS-2OMe (ambimodal),
#        TS-2OMe-cis, TS-2OMe-C2, TS-2OMe-2 (the (2+1) -> (4+1) rearrangement)
# The question is only whether each DFT saddle / minimum has a POLAR-1 counterpart nearby;
# result.json never certifies a stationary point or a connection (docs/GEOMETRIC_REFINEMENT.md).
#
# Resumable: a structure whose output directory exists is skipped (prrs optimize refuses to
# overwrite). One failing structure does not stop the rest; its exit code is logged.
#
# Usage: benchmarks/carbene_ambimodal/optimize_b08_b09.sh <out dir>
set -uo pipefail
cd /home/ruigengji/MLPATH
OUT=${1:?usage: optimize_b08_b09.sh <out dir>}
S=benchmarks/carbene_ambimodal/structures
export PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1
PY=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
mkdir -p "$OUT"

TS=(ts/TS-CCl2 ts/TS-CCl2-C4 ts/TS-2OMe ts/TS-2OMe-cis ts/TS-2OMe-C2 ts/TS-2OMe-2)
MIN=(olefins/1_2-bis_methylene_cycloheptane_7C carbenes/CCl2
     olefins/1-Cyclohexene-1-carboxaldehyde_6O carbenes/Dimethoxycarbene carbenes/Dimethoxycarbene_2
     products/6O-_2_+_1 products/6O-_4_+_1)

run() {  # run <structure> [--transition]
  local name=${1#*/} dest="$OUT/${1#*/}"
  [ -e "$dest" ] && { echo "SKIP $name (exists)"; return; }
  local status=0
  $PY -m prrs optimize "$S/$1.xyz" --output "$dest" \
    --calculator benchmarks.b04_sn2_chloride.calculator:polar1 \
    --convergence-set GAU --maxiter 300 "${@:2}" > "$OUT/$name.log" 2>&1 || status=$?
  echo "OPT_EXIT $name $status $(date '+%F %T')"
}
for s in "${TS[@]}"; do run "$s" --transition; done
for s in "${MIN[@]}"; do run "$s"; done
echo ALL_DONE
