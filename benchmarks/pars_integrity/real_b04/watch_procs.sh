#!/usr/bin/env bash
# Independent observer for the 2026-10-08 silent deaths on kasuga02: every 5 s, log korakuen's
# python/bash processes (state, elapsed, RSS) and GPU memory, in its own session so whatever
# takes the job down does not also take this down. Stops by itself after $1 seconds.
set -uo pipefail
end=$(( $(date +%s) + ${1:?seconds} ))
while [ "$(date +%s)" -lt "$end" ]; do
  echo "== $(date +%T)"
  ps -u korakuen -o pid,ppid,pgid,sid,stat,etime,rss,cmd | grep -E 'pars|gpu_batch' | grep -v -E 'grep|watch_procs' | cut -c1-160
  nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader 2>&1
  sleep 5
done
