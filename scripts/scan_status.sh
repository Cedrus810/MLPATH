#!/usr/bin/env bash
# What state is each seed actually in?
#
# network.json is NOT a completion marker: publish() rewrites it throughout a run with
# status "running", so its mere existence cannot tell a finished run from one that died
# mid-flight. The status field plus a liveness check can.
cd "$(dirname "$0")/.."
PY=${PY:-/home/ruigengji/miniforge3/envs/openmm_dev/bin/python}
printf "%-6s %-26s %-8s %-7s %s\n" seed status alive trials note
for d in runs/p2_T300_seed*; do
  [ -d "$d" ] || continue
  s=${d##*p2_T300_seed}
  if [ -f "$d/network.json" ]; then
    st=$("$PY" -c "import json;print(json.load(open('$d/network.json')).get('status','?'))" 2>/dev/null || echo unreadable)
  else
    st=no_network_json
  fi
  # pgrep sees the LOCAL process table. These runs share a filesystem with the machine
  # that launched them but not a process table, so on any other host "no process" means
  # "not visible from here" and nothing more. Reading it as "not running" once produced a
  # complete and entirely wrong story -- a 48-second kill window, a PBS walltime, eight
  # seeds declared interrupted -- and a recommendation to rm -rf seven directories whose
  # jobs were still running. So the host is checked first and liveness is simply not
  # claimed when it cannot be known.
  ran_on=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['versions']['platform'])" "$d/config.json" 2>/dev/null)
  here=$("$PY" -c "import platform;print(platform.platform())" 2>/dev/null)
  if [ -n "$ran_on" ] && [ "$ran_on" != "$here" ]; then
    alive="?"
  else
    alive=no
    pgrep -f "[-]-output $d " >/dev/null 2>&1 && alive=yes
    pgrep -f "[-]-output $d\$" >/dev/null 2>&1 && alive=yes
  fi
  n=$(wc -l < "$d/attempts.jsonl" 2>/dev/null || echo 0)
  # "no process and status still running" is TWO different states and calling both DEAD is
  # the same finished-versus-failed conflation this script exists to remove, one level down:
  #
  #   failed       something raised. The CLI prints "prrs: <Error>: ..." and/or a trial
  #                left an error.txt. There is a defect to find.
  #   interrupted  no error anywhere. The run was cut short from outside -- a job walltime,
  #                a qdel, a session ending. Nothing is broken and nothing needs fixing;
  #                it simply has not answered yet.
  #
  # The difference is mechanical, so it is decided mechanically rather than by tone.
  log="docs/experiments/$(basename "$d").out"
  note=""
  if [ "$st" = "running" ] && [ "$alive" = "?" ]; then
    # Not knowable from here. File freshness is the only signal available, and it is
    # labelled as the inference it is rather than dressed up as a liveness check.
    fresh=$(( ( $(date +%s) - $(stat -c %Y "$d/network.json" 2>/dev/null || echo 0) ) / 60 ))
    note="launched on another host ($(echo "$ran_on" | cut -c1-24)...); last wrote ${fresh} min ago"
    if [ "$fresh" -le 15 ]; then
      note="$note -- PROBABLY STILL RUNNING, do not delete"
    else
      note="$note -- check on that host before concluding anything"
    fi
  elif [ "$st" = "running" ] && [ "$alive" = "no" ]; then
    err=$(grep -m1 -E '^prrs: [A-Za-z]+Error' "$log" 2>/dev/null)
    if [ -n "$err" ] || [ -n "$(find "$d" -name error.txt -print -quit 2>/dev/null)" ]; then
      note="FAILED: ${err:-a trial wrote error.txt} -- find the defect"
    else
      note="INTERRUPTED: no process, no error. Cut short from outside; not a failure."
    fi
  fi
  [ "$st" = "no_network_json" ] && note="never published; log: $(tail -1 "$log" 2>/dev/null | cut -c1-70)"
  printf "%-6s %-26s %-8s %-7s %s\n" "$s" "$st" "$alive" "$n" "$note"
done
echo
echo "completed: $(grep -c completed <<<"$(for d in runs/p2_T300_seed*; do "$PY" -c "import json;print(json.load(open('$d/network.json')).get('status','?'))" 2>/dev/null; done)")"
