#!/usr/bin/env bash
# Run one prrs command on kasuga02, detached, and say where its log is.
#
# All computation goes to kasuga02 by standing decision (user, 2026-09-09), so this box
# keeps its GPU free. Nothing is copied: /home/ruigengji is the same NFS mount on both
# machines, which is also why `runs/<name>` means the same directory from either side.
#
# Two things do NOT carry over, and both were mistakes before they were flags:
#   - ssh does not cd. Without it the working directory is the remote login's home.
#   - ~/.cache/mace is per-machine and we log in as a different user there, so the
#     checkpoint is named explicitly. Same file: sha256 e5ccf583...bc69.
#
# Liveness is not claimed from here: this host shares the filesystem with kasuga02 but
# not its process table (see scan_status.sh, which had to learn that the hard way).
#
# Usage:  scripts/prrs_remote.sh search input.extxyz --output runs/foo --config c.json
#         scripts/prrs_remote.sh --wait search ...      block until it exits
set -euo pipefail

HOST=${PRRS_HOST:-korakuen@192.168.0.183}
DIR=${PRRS_REMOTE_DIR:-/home/ruigengji/MLPATH}
PY=${PRRS_REMOTE_PYTHON:-/home/ruigengji/miniforge3/envs/openmm_dev/bin/python}
MODEL=${PRRS_MACE_MODEL:-/home/ruigengji/MLP/mace/MACE-OFF24_medium.model}

WAIT=no
[ "${1:-}" = "--wait" ] && { WAIT=yes; shift; }
[ $# -gt 0 ] || { echo "usage: prrs_remote.sh [--wait] <prrs args...>" >&2; exit 2; }

LOG=.scratch/prrs_remote_$(date +%Y%m%d_%H%M%S).log
printf -v ARGS '%q ' "$@"

ssh -o BatchMode=yes "$HOST" "
  cd $DIR
  mkdir -p .scratch
  export PYTHONPATH=$DIR:$DIR/src
  export PRRS_MACE_MODEL=$MODEL
  setsid nohup $PY -m prrs $ARGS > $LOG 2>&1 < /dev/null &
  sleep 3
  pgrep -f 'python -m prrs' > /dev/null || { echo 'it exited immediately:'; cat $LOG; exit 1; }
"
echo "[prrs] launched on $HOST, log $DIR/$LOG (same path from here)"

if [ "$WAIT" = yes ]; then
  while ssh -o BatchMode=yes "$HOST" "pgrep -f 'python -m prrs' > /dev/null"; do sleep 60; done
  echo "[prrs] finished; tail of $LOG:"
  tail -3 "$LOG"
fi
