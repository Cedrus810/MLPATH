#!/bin/bash
# b05 QM chain for one TS (protocol of docs/experiments/b04_qm): IRC both ways ->
# OPT+Freq of both IRC ends -> wB97M-V/def2-TZVPD single points on TS and both minima.
# All r2SCAN-3c / TightSCF, total charge -1. IRC starts from the TS Freq Hessian (.hess)
# instead of recomputing it (same level, same geometry). Each step is skipped when its
# output already ends normally, so rerunning resumes.
#   usage: run_chain.sh NAME NPROCS      (NAME: new_10-12 | ctrl_6-9)
set -euo pipefail
N=$1; P=$2
O=/home/korakuen/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export PATH=$O:$PATH
export LD_LIBRARY_PATH=$O${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
cd "$(dirname "$0")"

done_ok() { [ -f "$1" ] && grep -q "ORCA TERMINATED NORMALLY" "$1"; }
run() {  # run INPUT_BASENAME
  if done_ok "$1.out"; then echo "skip $1"; return; fi
  echo "=== $(date '+%F %T') start $1"
  "$O/orca" "$1.inp" > "$1.out" 2>&1
  done_ok "$1.out" || { echo "FAILED $1"; exit 1; }
  echo "=== $(date '+%F %T') done $1"
}

head_r2scan() { printf '! r2SCAN-3c TightSCF %s\n%%pal nprocs %s end\n%%maxcore 2100\n' "$1" "$P"; }
head_wb()     { printf '! wB97M-V def2-TZVPD DEFGRID2 TightSCF\n%%pal nprocs %s end\n%%maxcore 2100\n' "$P"; }

# 1. IRC
{ head_r2scan IRC
  printf '%%irc MaxIter 100 InitHess read Hess_Filename "ts_%s.hess" end\n' "$N"
  printf '* xyzfile -1 1 ts_%s.xyz\n' "$N"; } > "irc_$N.inp"
run "irc_$N"

# 2. OPT+Freq of both ends
for side in F B; do
  { head_r2scan "OPT Freq"; printf '* xyzfile -1 1 irc_%s_IRC_%s.xyz\n' "$N" "$side"; } > "ep_${N}_$side.inp"
  run "ep_${N}_$side"
done

# 3. wB97M-V single points: TS and both optimised minima
{ head_wb; printf '* xyzfile -1 1 ts_%s.xyz\n' "$N"; } > "wb_ts_$N.inp"
run "wb_ts_$N"
for side in F B; do
  { head_wb; printf '* xyzfile -1 1 ep_%s_%s.xyz\n' "$N" "$side"; } > "wb_ep_${N}_$side.inp"
  run "wb_ep_${N}_$side"
done
echo "=== $(date '+%F %T') chain $N complete"
