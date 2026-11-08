#!/bin/bash
# Re-optimise one b05 QM endpoint that stopped at ORCA's 50-cycle limit in run_chain3
# (2026-11-01, user "可以"): start from that run's last geometry (ep_<TAG>.xyz), allow 300
# cycles, Freq at the converged minimum, then the wB97M-V/def2-TZVPD single point.
# Unlike run_chain*.sh, an OPT step counts as done only if the optimisation converged
# ("OPTIMIZATION RUN DONE"); ORCA ending normally without convergence is a failure.
# New basenames (ep2_*, wb_ep2_*) so the earlier outputs stay as they are. Rerun resumes.
#   usage: run_endpoint2.sh TAG NPROCS   (TAG: new_10-12_react | new_10-12_prod | ctrl_6-9_prod)
set -euo pipefail
T=$1; P=$2
O=$(dirname "$(command -v orca)")   # ORCA needs its full path for parallel runs
cd "$(dirname "$0")"

ended() { [ -f "$1" ] && grep -q "ORCA TERMINATED NORMALLY" "$1"; }
converged() { ended "$1" && grep -q "OPTIMIZATION RUN DONE" "$1"; }
run() {  # run BASENAME CHECK
  if "$2" "$1.out"; then echo "skip $1"; return; fi
  echo "=== $(date '+%F %T') start $1"
  "$O/orca" "$1.inp" > "$1.out" 2>&1 || true
  "$2" "$1.out" || { echo "FAILED $1 ($2)"; exit 1; }
  echo "=== $(date '+%F %T') done $1"
}

{ printf '! r2SCAN-3c TightSCF OPT Freq\n%%pal nprocs %s end\n%%maxcore 2100\n' "$P"
  printf '%%geom MaxIter 300 end\n'
  printf '* xyzfile -1 1 ep_%s.xyz\n' "$T"; } > "ep2_$T.inp"
run "ep2_$T" converged

{ printf '! wB97M-V def2-TZVPD DEFGRID2 TightSCF\n%%pal nprocs %s end\n%%maxcore 2100\n' "$P"
  printf '* xyzfile -1 1 ep2_%s.xyz\n' "$T"; } > "wb_ep2_$T.inp"
run "wb_ep2_$T" ended
echo "=== $(date '+%F %T') endpoint2 $T complete"
