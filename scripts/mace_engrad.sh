#!/usr/bin/env bash
# ORCA execs ProgExt with its own PATH; conda is not on it. Name the interpreter.
exec "${PRRS_ENGRAD_PYTHON:-/home/ruigengji/miniforge3/envs/openmm_dev/bin/python}" \
     "$(dirname "$(readlink -f "$0")")/mace_engrad.py" "$@"
