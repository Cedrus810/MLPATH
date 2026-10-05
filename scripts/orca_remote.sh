#!/usr/bin/env bash
# Run one ORCA job on kasuga02 with settings that match the machine.
#
# The first few jobs of 2026-09-04 were run with PAL32, no binding and
# %maxcore 3000 -- that asks for 32 x 3000 MB = 96 GB on a 93 GB box, and it
# spreads ranks across two NUMA nodes at random. Small jobs survived it by
# luck. This script exists so nothing has to remember the right numbers.
#
#   kasuga02: 2 sockets x 20 cores = 40 physical (80 with HT, not used)
#             NUMA node0 = 0-19,40-59   node1 = 20-39,60-79
#             93 GB RAM
#
# Usage:  scripts/orca_remote.sh <local .inp>  [nprocs]
# The input's %pal / %maxcore lines are rewritten to match `nprocs`.
set -euo pipefail

HOST=${ORCA_HOST:-korakuen@192.168.0.183}
ORCA=${ORCA_DIR:-/home/korakuen/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2}
PHYSICAL=40          # physical cores; hyperthread siblings are deliberately unused
RAM_MB=93000
FRACTION=75          # percent of RAM ORCA may claim, leaving room for the OS

INP=${1:?usage: orca_remote.sh <input.inp> [nprocs]}
N=${2:-32}
(( N <= PHYSICAL )) || { echo "nprocs $N exceeds $PHYSICAL physical cores" >&2; exit 2; }
MAXCORE=$(( RAM_MB * FRACTION / 100 / N ))

BASE=$(basename "$INP" .inp)
REMOTE=/tmp/orca_$BASE.$$

# Strip any PAL/maxcore the caller wrote and set them from the machine instead.
# The PALn shortcut only accepts 2/3/4/8/16/32/64 -- PAL40 is an input error -- so
# an arbitrary core count has to go through the %pal block.
BODY=$(sed -E 's/[[:space:]]*PAL[0-9]+//I; /^%maxcore/Id; /^%pal/,/^ *end/Id' "$INP")

ssh -o BatchMode=yes "$HOST" "mkdir -p $REMOTE"
{
  printf '%s\n' "$BODY" | head -1
  printf '%%pal nprocs %d end\n' "$N"
  printf '%%maxcore %d\n' "$MAXCORE"
  printf '%s\n' "$BODY" | tail -n +2
} | ssh -o BatchMode=yes "$HOST" "cat > $REMOTE/$BASE.inp"

ssh -o BatchMode=yes "$HOST" "
  export LD_LIBRARY_PATH=$ORCA:\${LD_LIBRARY_PATH:-}
  # One rank per physical core, bound, spread evenly over both NUMA nodes so
  # memory is local to the rank using it. OMP_NUM_THREADS=1 stops OpenBLAS from
  # opening a second layer of threads underneath MPI and oversubscribing.
  export OMPI_MCA_hwloc_base_binding_policy=core
  export OMPI_MCA_rmaps_base_mapping_policy=numa
  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  cd $REMOTE
  t0=\$(date +%s); $ORCA/orca $BASE.inp > $BASE.out 2>&1; t1=\$(date +%s)
  echo \"[orca] $BASE  PAL$N  maxcore ${MAXCORE}MB  \$((t1-t0))s  \$(grep -c 'TERMINATED NORMALLY' $BASE.out | sed s/1/OK/ | sed s/0/FAILED/)\"
  echo \"[orca] remote dir $REMOTE\"
"
