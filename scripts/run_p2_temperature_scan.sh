#!/usr/bin/env bash
# P2 temperature scan -- see docs/P2_TEMPERATURE_SCAN_PROTOCOL.md (frozen before running).
#
#   ./scripts/run_p2_temperature_scan.sh check          verify the knob first (REQUIRED)
#   ./scripts/run_p2_temperature_scan.sh run [PER_GPU]  run all 16 seeds
#   ./scripts/run_p2_temperature_scan.sh report         aggregate whatever has finished
#
# Concurrency saturates early and the ceiling is machine-specific, so run `probe` first.
# Measured:
#   RTX 2080 Ti  ceiling 3.2x   (7 workers 3.7x, 20 workers 3.2x -- no gain past ~7)
#   RTX 5080     ceiling 1.7x   (1 -> 24 workers: per-call 23.8 -> 328 ms, aggregate 1.72x)
#                               32 workers OOMs on 16 GB (20 of 32 failed)
# The FASTER card has the LOWER ceiling, which rules out GPU compute as the contended
# resource. Pinning one core and one thread per process changed nothing (1.70 -> 1.72), so
# it is not CPU thread contention either; the remaining suspect is CUDA context switching
# between processes, whose cost does not shrink when the card gets faster. CUDA MPS is the
# mechanism meant to remove exactly that, but on the 2080 Ti it measured 0.15x on
# launch-bound work -- do not enable it without re-measuring on the card in front of you.
#
# There are only 16 seeds, so concurrency above 16 does nothing here. On a single GPU
# prefer PER_GPU=8: total time is the same (throughput is saturated either way) but the
# first eight seeds -- the pre-registered n=8 prefix -- finish at half time instead of
# everything landing at the end.
#
# Seeds are spread across the GPUs named in CUDA_VISIBLE_DEVICES (default: all visible).
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
PY=${PY:-/home/ruigengji/miniforge3/envs/openmm_dev/bin/python}
SEEDS=(17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79)
PER_GPU=${2:-3}
CONF=docs/experiments/configs
INPUT=runs/p2_oxobutanal/p2_input.extxyz

# The checkpoint now travels with the repo (models/), so nothing here depends on whose home
# directory it happens to sit in. It was originally only under /home/kasuga/.cache, a home
# with mode 700, which meant any job running as a different user could not read it at all.
# Override for a machine that keeps it elsewhere:
#   MACE_MODEL=/abs/path/to/MACE-OFF24_medium.model ./scripts/run_p2_temperature_scan.sh run
MACE_MODEL=${MACE_MODEL:-$ROOT/models/MACE-OFF24_medium.model}
MACE_DEVICE=${MACE_DEVICE:-cuda}
# Written into the repo, not into TMPDIR, and never deleted. The previous version put it in
# $TMPDIR and removed it from an exit handler -- which fired as soon as the launcher finished
# handing out seeds, i.e. while the LAST seed's python was still starting and had not yet
# opened the file. Measured: seed 79 died with
#   FileNotFoundError: '/var/tmp/pbs.3899.yayoi/p2_mace_kwargs.1470241.json'
# On PBS the job's TMPDIR is also wiped at job end, so the run could not be reproduced
# afterwards either. Keeping the file is provenance, not litter: it records what ran.
KW=$CONF/mace_kwargs.generated.json

# pgrep -c prints a count AND exits 1 when the count is zero, so `|| echo 0` would append a
# second line and every numeric test on it fails with "integer expression expected".
# The bracket in '[p]rrs search' is load-bearing: the regex matches the text "prrs search"
# while a shell whose own argv contains the literal "[p]rrs search" does not, so this cannot
# count itself. Without it the helper reported 3 running searches on an idle machine.
running_searches() { pgrep -f '[p]rrs search' 2>/dev/null | wc -l | tr -d ' '; }

preflight() {
  local bad=0
  [ -x "$PY" ] || { echo "python not executable: $PY   (set PY=/path/to/python)"; bad=1; }
  [ -f "$INPUT" ] || { echo "missing input structure: $ROOT/$INPUT"; bad=1; }
  [ -f "$CONF/p2_round6.json" ] || { echo "missing base config: $CONF/p2_round6.json"; bad=1; }
  if [ ! -f "$MACE_MODEL" ]; then
    echo "MACE checkpoint not found: $MACE_MODEL"
    echo "  it is not in this repo. Point MACE_MODEL at the absolute path on this machine:"
    echo "    MACE_MODEL=/abs/path/MACE-OFF24_medium.model $0 ${1:-run}"
    bad=1
  fi
  command -v nvidia-smi >/dev/null || { echo "nvidia-smi not on PATH"; bad=1; }
  [ "$bad" -eq 0 ] || exit 5
  printf '{"model_paths": "%s", "device": "%s", "default_dtype": "float64"}\n' \
    "$MACE_MODEL" "$MACE_DEVICE" > "$KW"
  echo "  model: $MACE_MODEL"
  local got; got=$(sha256sum "$MACE_MODEL" | cut -c1-16)
  echo "  sha256: $got..."
  if [ "$got" != "e5ccf5837f685899" ]; then
    echo "  WARNING: this is not the checkpoint rounds 3-6 used (e5ccf5837f685899...)."
    echo "  Results will not be comparable with those runs. Continuing anyway."
  fi
}

emit_config() {   # emit_config <seed> <max_trials|-> <path>
  "$PY" - "$1" "$2" "$3" <<'PYEOF'
import json, sys
sys.path.insert(0, "src")
from prrs.config import SearchConfig
seed, trials, out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
cfg = json.load(open("docs/experiments/configs/p2_round6.json"))
cfg["temperature_K"] = 300.0
cfg["seed"] = seed
if trials != "-":
    cfg["max_trials"] = int(trials)
SearchConfig.from_dict(cfg)          # validate before anything is launched
json.dump(cfg, open(out, "w"), indent=1, sort_keys=True)
PYEOF
}

# One thread and one dedicated core per process. A 12-atom model gets nothing from intra-op
# parallelism, and without this N processes each open a thread pool and fight over the cores:
# measured on an RTX 5080, per-call time grew 13.9x from 1 to 24 workers while aggregate
# throughput rose only 1.70x, so most of each call was contended rather than GPU-bound.
# Override the thread count with SCAN_THREADS; drop the affinity with SCAN_PIN=0.
CORE=0
NCORES=$(nproc)
launch() {        # launch <seed> <config> <outdir> <gpu>
  local pin=()
  if [ "${SCAN_PIN:-1}" != "0" ] && command -v taskset >/dev/null; then
    pin=(taskset -c "$(( CORE % NCORES ))"); CORE=$(( CORE + 1 ))
  fi
  CUDA_VISIBLE_DEVICES="$4" PYTHONPATH=src \
  OMP_NUM_THREADS=${SCAN_THREADS:-1} MKL_NUM_THREADS=${SCAN_THREADS:-1} \
  OPENBLAS_NUM_THREADS=${SCAN_THREADS:-1} NUMEXPR_NUM_THREADS=${SCAN_THREADS:-1} \
  nohup "${pin[@]}" "$PY" -m prrs search "$INPUT" \
    --output "$3" --config "$2" --calculator prrs.calculators:mace_factory \
    --calculator-kwargs "$KW" > "docs/experiments/$(basename "$3").out" 2>&1 &
  echo "  seed $1 -> $3 (gpu $4) pid $!"
}

gpus() {
  if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then echo "$CUDA_VISIBLE_DEVICES" | tr ',' ' '
  else nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' '; fi
}

case "${1:-}" in
check)
  echo "== step 0: does temperature actually perturb the run? =="
  # Two check runs at once destroy each other: the second one's rm -rf deletes the
  # directories the first one is still writing into, and then neither produces trials.
  # Measured the hard way, so it is an interlock rather than a comment.
  preflight check
  if [ "$(running_searches)" -gt 0 ]; then
    echo "refusing: a prrs search is already running. Wait for it or kill it first:"
    pgrep -af 'prrs search' | sed 's/^/  /'
    exit 3
  fi
  for s in 17 19; do
    rm -rf "runs/knob_T300_s$s"
    emit_config "$s" 3 "/tmp/knob_T300_s$s.json"
    launch "$s" "/tmp/knob_T300_s$s.json" "runs/knob_T300_s$s" "$(gpus | awk '{print $1}')"
  done
  wait
  for s in 17 19; do
    if [ ! -s "runs/knob_T300_s$s/attempts.jsonl" ]; then
      echo "seed $s produced no trials; its log:"; tail -5 "docs/experiments/knob_T300_s$s.out"
      exit 4
    fi
  done
  "$PY" - <<'PYEOF'
import json, sys
from pathlib import Path
rows = {}
for seed in (17, 19):
    p = Path(f"runs/knob_T300_s{seed}/attempts.jsonl")
    rows[seed] = [json.loads(l) for l in p.read_text().splitlines()] if p.exists() else []
for seed, rs in rows.items():
    print(f"-- seed {seed} --")
    for r in rs:
        print(f"   {r['attempt_id']} {r['probe']['family']}{r['probe']['indices']} "
              f"amp={r['probe']['amplitude']} "
              f"injected_KE={r.get('injected_kinetic_eV'):.6f} "
              f"drift={r.get('max_nve_drift_eV_atom'):.2e} "
              f"status={r['status']} target={r.get('target_microstate')}")
# The right observable is the kinetic energy at the "perturbed" frame, which is where
# thermal_momenta has just been applied. injected_kinetic_eV is the KICK injection and is
# identically zero for a displace probe however hot the run is -- reading it as "is the
# thermostat on" reports NOT VERIFIED on a run whose thermostat is plainly working.
def perturbed_ke(seed):
    out = []
    for trial in sorted(Path(f"runs/knob_T300_s{seed}/trials").iterdir()):
        obs = trial / "observations.jsonl"
        if not obs.exists():
            continue
        for line in obs.read_text().splitlines():
            row = json.loads(line)
            if row["phase"] == "perturbed":
                out.append(row["kinetic_eV"])
                break
    return out

ke = {s: perturbed_ke(s) for s in rows}
res = {s: [(r["status"], r.get("target_microstate")) for r in rs] for s, rs in rows.items()}
print(f"\nperturbed-frame kinetic energy  seed 17: {ke.get(17)}")
print(f"                                seed 19: {ke.get(19)}")
ke_differs = bool(ke.get(17)) and ke.get(17) != ke.get(19) and any(ke[17])
out_differs = res.get(17) != res.get(19)
print(f"\nthermal momenta present and seed-dependent: {ke_differs}")
print(f"at least one trial's outcome differs:       {out_differs}")
print("  (outcomes may legitimately agree on only three trials; the kinetic-energy")
print("   criterion is the one that decides whether the knob moves the system)")
if ke_differs:
    print("\nKNOB VERIFIED -- the scan can measure something. Run:  ./scripts/run_p2_temperature_scan.sh run")
else:
    print("\nKNOB NOT VERIFIED -- temperature does not move the tested mechanism either.")
    print("STOP. Running 16 seeds would buy an answer that is already determined.")
    sys.exit(1)
PYEOF
  ;;
probe)
  # Find this card's saturation point instead of guessing it. Aggregate throughput stops
  # rising once the GPU-resident fraction of one process is filled; on a 2080 Ti that was
  # 3.2x at both 7 and 20 workers. A faster card does the GPU part in less wall time, so
  # that fraction is smaller and the ceiling higher -- by how much is a measurement.
  preflight probe
  "$PY" scripts/_saturation_probe.py "$KW"
  ;;
run)
  echo "== launching 16 seeds, $PER_GPU per GPU =="
  preflight run
  GPU_LIST=($(gpus)); n=${#GPU_LIST[@]}; slots=$(( n * PER_GPU ))
  echo "   gpus: ${GPU_LIST[*]}   concurrency: $slots"
  i=0
  for s in "${SEEDS[@]}"; do
    out="runs/p2_T300_seed$s"
    # An existing directory means one of two different things and they must not be
    # conflated: a finished run (skip it) or a run that died before writing network.json
    # (needs redoing, and the CLI will refuse to reuse a non-empty directory).
    if [ -d "$out" ]; then
      # network.json is NOT a completion marker: publish() rewrites it throughout a run
      # with status "running". Testing only for its existence reported a run that died
      # mid-flight as "already run" -- the same finished/failed conflation this branch was
      # added to remove, just moved from the directory to the file.
      st=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1])).get('status','?'))" \
           "$out/network.json" 2>/dev/null || echo unreadable)
      if [ "$st" = "completed" ]; then
        echo "  seed $s: completed, skipping"
      elif pgrep -f "[-]-output $out " >/dev/null 2>&1; then
        echo "  seed $s: still running (status $st), leaving it alone"
      elif [ -f "$out/network.json" ]; then
        echo "  seed $s: DEAD -- status '$st' with no process. $(wc -l < "$out/attempts.jsonl" 2>/dev/null || echo 0) trials done."
        echo "      to redo it:  rm -rf $out && $0 run $PER_GPU"
      else
        echo "  seed $s: FAILED EARLIER, no network.json. Its log ends with:"
        tail -2 "docs/experiments/$(basename "$out").out" 2>/dev/null | sed 's/^/      /'
        echo "      to redo it:  rm -rf $out && $0 run $PER_GPU"
      fi
      continue
    fi
    while [ "$(running_searches)" -ge "$slots" ]; do sleep 20; done
    cfg="$CONF/p2_T300_seed$s.json"; emit_config "$s" - "$cfg"
    launch "$s" "$cfg" "$out" "${GPU_LIST[$(( i % n ))]}"
    i=$(( i + 1 ))
  done
  echo "all launched; wait for them, then: $0 report"
  echo "then check they share one implementation:  $0 provenance"
  ;;
provenance)
  # Refuses to let a mixed-implementation batch be quoted as one measurement.
  "$PY" scripts/_provenance.py runs/p2_T300_seed "${SEEDS[@]}"
  ;;
report)
  args=(); for s in "${SEEDS[@]}"; do args+=("$s:runs/p2_T300_seed$s"); done
  PYTHONPATH=src "$PY" docs/experiments/p2_multiseed.py "${args[@]}"
  echo; echo "per-seed criteria:"
  for s in "${SEEDS[@]}"; do
    [ -f "runs/p2_T300_seed$s/network.json" ] || continue
    echo "-- seed $s --"
    PYTHONPATH=src "$PY" docs/experiments/p2_check.py "runs/p2_T300_seed$s/network.json" | tail -3
  done
  ;;
*)
  sed -n '2,12p' "$0"
  echo "  probe                      measure this card's saturation point"
  exit 2;;
esac
