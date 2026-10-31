#!/bin/bash
# Records which snappyHexMesh variant produced the mesh of each steady CFD run and re-checks the
# production meshes with checkMesh.
# Usage: scripts/record_production_meshes.sh        (writes results/production_meshes.txt)
# A run's mesh "source" is every feas* case whose constant/polyMesh/points has the same md5 as the run.
ROOT=~/FEA-portfolio/04-crm-wing-cfd-structural
BASE=$ROOT/simulation/cfd_v1
OUT=$ROOT/results/production_meshes.txt
RUNS="u_L8_p0p36 u_p0p36 u_L10_p0p36 u_defl_p0p36 u_defl_a2p07 u_p0p5 u_p0p3 a_m0p5 a_m1p5 b_m0p8 t_m0p5"
CHECK="u_L8_p0p36 u_p0p36 u_L10_p0p36 u_defl_p0p36"
SRCS="feas_L8_novolq feas_L9_novolq feas_L10_novolq feas_L9_novolq_defl"
KEYS='level|maxLocalCells|maxGlobalCells|minVol |maxNonOrtho|BoundarySkewness|InternalSkewness'

declare -A SUM
for c in "$BASE"/feas*/; do
  cf=${c}constant/polyMesh/points
  [ -f "$cf" ] && SUM[$(basename "$c")]=$(md5sum < "$cf" | cut -c1-32)
done

{
echo "# Production mesh record, generated $(date -u +%Y-%m-%dT%H:%MZ) by scripts/record_production_meshes.sh"
echo "# source = feas* case folder whose constant/polyMesh/points has the same md5 as the run"
echo
echo "## Run -> mesh"
for r in $RUNS; do
  f=$BASE/runs/$r/constant/polyMesh/points
  if [ ! -f "$f" ]; then echo "$r  no points file"; continue; fi
  m=$(md5sum < "$f" | cut -c1-32)
  n=$(grep -m1 -o 'nCells:[0-9]*' "$BASE/runs/$r/constant/polyMesh/owner")
  src=""
  for k in "${!SUM[@]}"; do [ "${SUM[$k]}" = "$m" ] && src="$src $k"; done
  echo "$r  $n  md5 ${m:0:12}  source:${src:- none}"
done
echo
echo "## snappyHexMeshDict settings of the source cases"
for s in $SRCS; do
  echo "--- $s"
  grep -n -E "$KEYS" "$BASE/$s/system/snappyHexMeshDict"
done
echo
echo "## checkMesh (constant mesh of the run)"
for r in $CHECK; do
  echo "--- $r"
  checkMesh -case "$BASE/runs/$r" -constant </dev/null 2>&1 | grep -E "cells:|Max aspect|Min volume|non-orthogonality|Max skewness|Mesh OK|Failed|\*\*\*|Concave|Zero|Negative|command not found"
done
} > "$OUT" 2>&1
echo "written $OUT ($(wc -l < "$OUT") lines)"
