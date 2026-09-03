#!/usr/bin/env bash
# P3 heavy jobs. Sharded worker pool, PER_GPU=8, following scripts/run_p2_temperature_scan.sh.
#
# Concurrency numbers are measured, not guessed (scripts/_saturation_probe.py):
#   RTX 2080 Ti  ceiling 3.2x   (7 workers 3.7x, 20 workers 3.2x -- no gain past ~7)
#   RTX 5080     ceiling 1.7x   (1 -> 24 workers: per-call 23.8 -> 328 ms, aggregate 1.72x)
# The faster card has the LOWER ceiling, which is what rules out GPU compute as the
# contended resource. One thread and one pinned core per process removed CPU-thread
# contention as a cause too (1.70x -> 1.72x, i.e. nothing). Do NOT enable CUDA MPS: it
# measured 0.15x, 6.8x SLOWER. The leverage is more cards, not more processes per card.
#
# PER_GPU=8 is P2's default and is kept: past ~7 the aggregate does not rise, but the first
# shards still finish sooner, and every job here checkpoints per shard.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY="${PY:-python}"
MODEL="${MODEL:-models/MACE-OFF24_medium.model}"
PER_GPU="${PER_GPU:-8}"
B=scripts/p3_node_bundle
LOG="${LOG:-runs/p3_node_logs}"; mkdir -p "$LOG"

usage() {
  cat <<USE
usage: $0 {ts|paths|pes|figure|quench-ab|models|landscape|all|merge|status} [grid]

  ts         first-order saddle, barrier, imaginary freq, two-sided descent (8 amplitudes)
  paths      record the delivered path (climb + both descents) with named coordinates,
             so it can be projected onto the surface. Failed descents keep their frames.
  figure     build the two-panel figure from whatever is on disk (CPU only)
  quench-ab  2x2 over (quench_fmax, quench_steps). Tolerance alone is NOT one variable:
             at a fixed step budget it also decides whether the criterion is reachable.
  models     the same reaction under MACE-OFF24, omol-0-xl-4M and POLAR-1-M. Relative
             quantities only -- absolute energies across models are meaningless.
  pes        rigid + relaxed (q_PT, r_OO) surface                            (grid rows)
  landscape  all 300 conformers relaxed AND Hessian-confirmed
  all        ts, then pes, then landscape
  merge      assemble shards and report the global numbers
  status     what is running, what is done, what is missing

env: PY=<python>  MODEL=<path>  PER_GPU=$PER_GPU  PATHFILE=<extxyz for the pes 2nd panel>
     P3_QUENCH=<fmax, default 2e-3>   P3_QUENCH_STEPS=<default 500>
     P3_MACE_DIR=<dir holding the three checkpoints>
     P3_AMPS=<comma list, models job only>   P3_WORKER_MIB=<footprint override>
USE
  exit 2
}
[ $# -ge 1 ] || usage
GRID="${2:-48}"

gpus() {
  if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then echo "$CUDA_VISIBLE_DEVICES" | tr ',' ' '
  else nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' '; fi
}
GPU_LIST=($(gpus)); NGPU=${#GPU_LIST[@]}; SLOTS=$(( NGPU * PER_GPU ))

# Measured footprint: one worker on this 22-atom system holds ~1350 MiB of a float64
# MACE-OFF24 context. At PER_GPU=8 that is 10780 of 11264 MiB on an RTX 2080 Ti -- 96% of
# the card, with a ninth process OOMing outright. So the slot count has a memory ceiling
# as well as the throughput ceiling, and it is the memory one that bites first on 11 GB.
# Checked rather than assumed, and only warned about: the operator may know better than a
# hard-coded footprint, and refusing to start would be worse than saying so.
capacity_check() {
  # Read the override HERE, not at the top of the script. Assigning P3_WORKER_MIB as a
  # command prefix cannot affect a variable that was already expanded during startup --
  # the `models` job set 4096 that way and the check still printed 1350.
  local worker_mib="${P3_WORKER_MIB:-1350}"
  local total per_gpu_fit
  total=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1) || return 0
  [ -n "$total" ] || return 0
  per_gpu_fit=$(( total / worker_mib ))
  echo "   gpu memory ${total} MiB, ~${worker_mib} MiB/worker -> about ${per_gpu_fit} fit per gpu"
  if [ "$PER_GPU" -gt "$per_gpu_fit" ]; then
    echo "   WARNING: PER_GPU=$PER_GPU exceeds that. Expect CUDA OOM on the later shards."
    echo "            Either PER_GPU=$per_gpu_fit, or set P3_WORKER_MIB if the footprint differs."
  fi
}

# One thread and one dedicated core per process. A 22-atom model gets nothing from intra-op
# parallelism, and without this N processes each open a thread pool and fight over the
# cores: on an RTX 5080 per-call time grew 13.9x from 1 to 24 workers while aggregate
# throughput rose only 1.70x, so most of each call was contention, not GPU work.
# Override with P3_THREADS; drop the pinning with P3_PIN=0.
CORE=0
NCORES=$(nproc)
PATTERN='p3_(ts_barrier|pes2d|landscape_full|three_model|paths)\.py'
# Live-worker count from the PIDs we launched, NOT from pgrep.
# Two reasons pgrep was wrong here. First, `pgrep -fc` prints "0" AND exits 1 when nothing
# matches, so `pgrep -fc ... || echo 0` emitted TWO lines and every arithmetic test
# downstream failed with "integer expected". Second and worse, -f matches any command line
# containing the script name -- including the operator's own monitoring shells -- so the
# count read 11 while 8 workers were running, and a slot limiter that over-counts stalls
# instead of launching. Tracking the PIDs is exact and cannot collide with anything.
PIDS=()
running() {
  local alive=0 pid
  for pid in "${PIDS[@]:-}"; do
    [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null && alive=$(( alive + 1 ))
  done
  echo "$alive"
}

launch() {        # launch <script> <tag> <shard> <nshard> <gpu> [extra args...]
  local script=$1 tag=$2 shard=$3 nshard=$4 gpu=$5; shift 5
  local pin=()
  if [ "${P3_PIN:-1}" != "0" ] && command -v taskset >/dev/null; then
    pin=(taskset -c "$(( CORE % NCORES ))"); CORE=$(( CORE + 1 ))
  fi
  # mkdir here, not only at the top: the quench grid overrides LOG per group with a
  # prefix assignment, and a missing directory makes the redirect fail, which means the
  # command never runs at all -- eight "launched" lines and zero workers.
  mkdir -p "$LOG"
  CUDA_VISIBLE_DEVICES="$gpu" PYTHONPATH=src \
  OMP_NUM_THREADS=${P3_THREADS:-1} MKL_NUM_THREADS=${P3_THREADS:-1} \
  OPENBLAS_NUM_THREADS=${P3_THREADS:-1} NUMEXPR_NUM_THREADS=${P3_THREADS:-1} \
  nohup "${pin[@]}" "$PY" -u "$B/$script" "$@" --shard "$shard/$nshard" \
    > "$LOG/${tag}_shard$(printf %02d "$shard").out" 2>&1 &
  PIDS+=("$!")
  echo "  $tag shard $shard/$nshard -> gpu $gpu  pid $!"
}

pool() {          # pool <script> <tag> <nshard> [extra args...]
  local script=$1 tag=$2 nshard=$3; shift 3
  echo "== $tag: $nshard shards over ${NGPU} gpu(s) [${GPU_LIST[*]}], $PER_GPU per gpu =="
  capacity_check
  local i
  for (( i=0; i<nshard; i++ )); do
    # Slot limiting rather than launching all at once: past the saturation point extra
    # processes only add context-switch cost, and on a 16 GB card 32 concurrent OOMed.
    while [ "$(running)" -ge "$SLOTS" ]; do sleep 15; done
    launch "$script" "$tag" "$i" "$nshard" "${GPU_LIST[$(( i % NGPU ))]}" "$@"
  done
  wait
  echo "== $tag: all $nshard shards exited  $(date -Is) =="
}

case "$1" in
  ts)        pool p3_ts_barrier.py     ts        "$SLOTS" "$MODEL" runs/p3_ts_barrier ;;
  landscape) pool p3_landscape_full.py landscape "$SLOTS" "$MODEL" runs/p3_landscape_full ;;
  pes)
    # Shards are rows, so more shards than rows would leave workers idle.
    n=$(( SLOTS < GRID ? SLOTS : GRID ))
    # Default the second panel to the recorded descent if `paths` has already run: the
    # panel is what makes the projection honest about energy, so it should not be opt-in
    # by accident.
    pf="${PATHFILE:-}"
    if [ -z "$pf" ] && [ -f runs/p3_paths/MACE-OFF24_medium/descent-1/descent.extxyz ]; then
      pf=runs/p3_paths/MACE-OFF24_medium/descent-1/descent.extxyz
      echo "   PATHFILE defaulted to $pf"
    fi
    pool p3_pes2d.py pes "$n" "$GRID" "$MODEL" runs/p3_pes2d "$pf"
    ;;
  all)
    # TS first on purpose: cheapest and most decisive. If there is no first-order saddle on
    # the transfer coordinate then the surface and the landscape are describing the wrong
    # thing, and the operator should see that before spending hours on them.
    pool p3_ts_barrier.py ts "$SLOTS" "$MODEL" runs/p3_ts_barrier
    "$PY" "$B/p3_merge.py" runs || echo "  (ts shards incomplete -- see above)"
    n=$(( SLOTS < GRID ? SLOTS : GRID ))
    pool p3_pes2d.py pes "$n" "$GRID" "$MODEL" runs/p3_pes2d "${PATHFILE:-}"
    pool p3_landscape_full.py landscape "$SLOTS" "$MODEL" runs/p3_landscape_full
    "$PY" "$B/p3_merge.py" runs
    ;;
  paths)
    # Records the delivered path -- climb plus both descents -- so there is something to
    # draw on the surface. A FAILED descent still writes its frames, and those frames are
    # the point: `side +1` has never reached A, and where that motion stalls in
    # (q_PT, r_OO) is visible nowhere else. quench_steps defaults to 5000 here so the
    # descent is not budget-bound (see the quench-ab note).
    export P3_WORKER_MIB="${P3_WORKER_MIB:-4096}"
    pool p3_paths.py paths 3 runs/p3_paths
    ;;
  figure)
    # Pure CPU, reads what is on disk. Safe to run anywhere.
    exec "$PY" "$B/p3_figure.py" runs
    ;;
  models)
    # Sharded by model, three shards: each relaxes its own endpoints once and then walks
    # every amplitude, so endpoint work is not repeated per amplitude. omol-0-xl-4M is a
    # 414 MB checkpoint and its resident footprint is well above OFF24's ~1350 MiB, so
    # three concurrent models want a roomy card -- this is the job to give the 16-24 GB
    # node rather than an 11 GB desktop.
    export P3_WORKER_MIB="${P3_WORKER_MIB:-4096}"
    pool p3_three_model.py models 3 runs/p3_three_model
    ;;
  quench-ab)
    # A 2x2 over (quench_fmax, quench_steps), not a scan of the tolerance alone.
    #
    # The first version of this varied only the tolerance and was wrong: quench_fmax sets
    # the convergence criterion AND, at a fixed 500-step FIRE budget, whether that
    # criterion is reachable. Going 2e-3 -> 2e-4 changed the two-sided-descent failure from
    # "soft_mode_moved_on_the_last_round" to "fmax_not_reached" for ALL THREE models, i.e.
    # it measured the budget. The cell (2e-4, 500) is that already-known result and is
    # included so the grid is complete and the confound is visible rather than argued about.
    #
    # Reading the grid:
    #   (2e-3, 500)   the protocol as it stands -- the baseline
    #   (2e-3, 5000)  budget alone. If this changes anything, the baseline was budget-bound.
    #   (2e-4, 500)   known: fmax_not_reached. Kept as the confounded cell.
    #   (2e-4, 5000)  tolerance with the budget to reach it -- the actual question.
    for cell in "0.002 500" "0.002 5000" "0.0002 500" "0.0002 5000"; do
      set -- $cell; q=$1; st=$2
      tag="q${q}_s${st}"
      echo "######## quench_fmax = $q  quench_steps = $st  ($(date -Is)) ########"
      rm -rf runs/p3_ts_barrier
      P3_QUENCH="$q" P3_QUENCH_STEPS="$st" LOG="runs/p3_node_logs_$tag" \
        pool p3_ts_barrier.py "ts_$tag" "$SLOTS" "$MODEL" runs/p3_ts_barrier
      if [ -d runs/p3_ts_barrier ]; then
        rm -rf "runs/p3_ts_barrier_$tag"
        mv runs/p3_ts_barrier "runs/p3_ts_barrier_$tag"
      else
        echo "   WARNING: no runs/p3_ts_barrier produced -- every worker died before writing."
        echo "            Check runs/p3_node_logs_$tag/ before reading anything else."
      fi
      echo "######## cell $tag done ($(date -Is)) ########"
    done
    echo
    for tag in q0.002_s500 q0.002_s5000 q0.0002_s500 q0.0002_s5000; do
      [ -d "runs/p3_ts_barrier_$tag" ] || { echo "=== $tag MISSING ==="; continue; }
      echo "=== $tag ==="
      rm -rf runs/p3_ts_barrier
      cp -r "runs/p3_ts_barrier_$tag" runs/p3_ts_barrier
      "$PY" "$B/p3_merge.py" runs | sed -n '/TS . barrier/,/merged/p'
      rm -rf runs/p3_ts_barrier
    done
    ;;
  merge)  exec "$PY" "$B/p3_merge.py" runs ;;
  status)
    echo "== processes =="
    # -f over-matches (it sees monitoring shells too); fine for a human listing,
    # which is why the slot limiter uses tracked PIDs instead.
    pgrep -af -- "$PATTERN" | sed 's/^/  /' || echo "  none running"
    echo "== shard files =="
    for d in runs/p3_ts_barrier runs/p3_pes2d runs/p3_landscape_full; do
      [ -d "$d" ] || continue
      c=$(ls "$d" 2>/dev/null | grep -c shard) || c=0
      printf "  %-28s %s shard artefacts\n" "$d" "${c:-0}"
    done
    echo "== logs with a traceback =="
    grep -l "Traceback\|SystemExit\|CUDA out of memory" "$LOG"/*.out 2>/dev/null | sed 's/^/  /' \
      || echo "  none"
    ;;
  *) usage ;;
esac
