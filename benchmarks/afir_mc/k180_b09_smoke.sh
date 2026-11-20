#!/bin/bash
# 2026-11-20 b09 smoke on 180 (docs/B09_DISTANCE_2026-11-20.md §5 A-B; user "180跑").
# 1) unbiased input (SI reactant conformers, TS-2OMe atom order); 2) SC-AFIR with the
# 2026-11-15 defaults, 2 random orientations, small budget -- only to see the pipeline run on
# b09 and measure per-step cost. Not a measurement of recall.
set -euo pipefail
cd /home/ruigengji/MLPATH
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
IN=benchmarks/carbene_ambimodal/results/2026-11-20/b09_afir_input
D=benchmarks/afir_mc/results/2026-11-20/b09_smoke180
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 AFIR_MAX_ATOMS_PER_CALL=512
mkdir -p $IN $D
echo "=== $(date '+%F %T') build input"
if [ ! -s $IN/input.extxyz ]; then
  $P -u benchmarks/carbene_ambimodal/build_b09_afir_input.py $IN > $IN/build.log 2>&1
fi
SHELL_R=$($P -c "import json;print(json.load(open('$IN/build.json'))['complex']['suggested_shell'])")
echo "=== $(date '+%F %T') shell $SHELL_R; AFIR smoke"
$P -u benchmarks/afir_mc/run_mc_afir_batch.py --system b09 --mode sc --sc-poses 2 \
  --out $D/b09_sA --seed 20261023 --shell $SHELL_R --budget 60000 --wave 32 > $D/b09_sA.log 2>&1 || echo "b09 FAILED"
echo "all done $(date '+%F %T')"
