#!/bin/bash
# b05 QM chain, second form (2026-11-01, user "可以"): the product-side IRCs were stopped
# after their geometry showed E2 products (Cl- gone, water formed) on a long flat descent,
# so product minima are optimised from the last IRC frame. Same levels as run_chain.sh /
# docs/experiments/b04_qm: r2SCAN-3c TightSCF, charge -1; wB97M-V/def2-TZVPD DEFGRID2
# single points. Each step is skipped when its output ends normally; rerun resumes.
#   usage: run_chain2.sh NAME NPROCS      (NAME: new_10-12 | ctrl_6-9)
set -euo pipefail
N=$1; P=$2
O=/home/korakuen/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2
export PATH=$O:$PATH
export LD_LIBRARY_PATH=$O${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
cd "$(dirname "$0")"

done_ok() { [ -f "$1" ] && grep -q "ORCA TERMINATED NORMALLY" "$1"; }
run() {
  if done_ok "$1.out"; then echo "skip $1"; return; fi
  echo "=== $(date '+%F %T') start $1"
  "$O/orca" "$1.inp" > "$1.out" 2>&1
  done_ok "$1.out" || { echo "FAILED $1"; exit 1; }
  echo "=== $(date '+%F %T') done $1"
}
r2() { printf '! r2SCAN-3c TightSCF %s\n%%pal nprocs %s end\n%%maxcore 2100\n' "$1" "$P"; }
wb() { printf '! wB97M-V def2-TZVPD DEFGRID2 TightSCF\n%%pal nprocs %s end\n%%maxcore 2100\n' "$P"; }
opt() {  # opt NAME START_XYZ
  { r2 "OPT Freq"; printf '* xyzfile -1 1 %s\n' "$2"; } > "$1.inp"; run "$1"
}
sp() {  # sp NAME XYZ
  { wb; printf '* xyzfile -1 1 %s\n' "$2"; } > "$1.inp"; run "$1"
}

case "$N" in
  new_10-12)
    # reactant side not walked yet: backward IRC only
    { r2 IRC
      printf '%%irc Direction backward MaxIter 100 InitHess read Hess_Filename "ts_%s.hess" end\n' "$N"
      printf '* xyzfile -1 1 ts_%s.xyz\n' "$N"; } > "ircB_$N.inp"
    run "ircB_$N"
    opt "ep_${N}_prod" "prod_${N}_start.xyz"
    opt "ep_${N}_react" "ircB_${N}_IRC_B.xyz"
    ;;
  ctrl_6-9)
    opt "ep_${N}_react" "react_${N}_start.xyz"
    opt "ep_${N}_prod" "prod_${N}_start.xyz"
    ;;
  *) echo "unknown $N"; exit 2 ;;
esac
sp "wb_ts_$N" "ts_$N.xyz"
sp "wb_ep_${N}_react" "ep_${N}_react.xyz"
sp "wb_ep_${N}_prod" "ep_${N}_prod.xyz"
echo "=== $(date '+%F %T') chain2 $N complete"
