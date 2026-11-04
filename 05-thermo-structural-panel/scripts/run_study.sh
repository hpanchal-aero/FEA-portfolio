#!/bin/bash
# Project 05 study driver. Runs a named group of cases from the study matrix
# (docs/specification.md, Amendments 12 and 14) in order and stops at the first
# case with a failed check (exit status 1) or a triggered stop rule (exit status 2).
# Usage: bash scripts/run_study.sh imperfection | thickness | mesh
cd "$(dirname "$0")/.." || exit 1
PY=/usr/bin/python3
case "$1" in
  imperfection)
    CASES=("--h 2 --bc rest --imp 0.01" "--h 2 --bc rest --imp 0.5") ;;
  thickness)
    CASES=("--h 3 --bc free --imp 0.1" "--h 3 --bc clamp --imp 0.1" "--h 1 --bc free --imp 0.1" "--h 1 --bc clamp --imp 0.1") ;;
  mesh)
    CASES=("--h 2 --bc rest --imp 0.1 --n 40") ;;
  *)
    echo "usage: bash $0 imperfection | thickness | mesh"
    exit 1 ;;
esac
for c in "${CASES[@]}"; do
  echo "=== $PY scripts/run_case.py $c"
  if ! $PY scripts/run_case.py $c; then
    echo "STOPPED: case failed or stop rule triggered: $c"
    exit 1
  fi
done
echo "GROUP $1 COMPLETE"
