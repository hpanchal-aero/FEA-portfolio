#!/bin/bash
# Surface max-level sweep for the snappyHexMesh feasibility mesh.
# usage: [QTIGHT=1] [SEDX='sed -E expressions'] [TAG=_name] snappy_level_sweep.sh L [L ...]
# For each L a clean case simulation/cfd_v1/feas_L<L><TAG> is built from simulation/cfd_v1/feas
# (system, triSurface, feature edges) with surface level (7 L), features level L,
# maxLocalCells 2M, maxGlobalCells 3M.
# QTIGHT=1 : maxBoundarySkewness 4, maxInternalSkewness 3.5, relaxed maxNonOrtho 70.
# SEDX     : extra sed -E expressions applied to system/snappyHexMeshDict after the standard edits.
# Summary is appended to results/cfd_snappy_level_sweep<TAG>.txt.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$ROOT/simulation/cfd_v1/${BASECASE:-feas}
TAG=${TAG:-}
OUT=$ROOT/results/cfd_snappy_level_sweep$TAG.txt
PY=/usr/bin/python3
if [ "$#" -lt 1 ]; then echo "usage: [QTIGHT=1] [SEDX=...] [TAG=_name] $0 L [L ...]"; exit 1; fi
mkdir -p "$ROOT/results"
for L in "$@"; do
  D=$ROOT/simulation/cfd_v1/feas_L$L$TAG
  echo "=== surface max level $L tag '$TAG' QTIGHT=${QTIGHT:-0} SEDX='${SEDX:-}' : start $(date '+%F %T') ===" >> "$OUT"
  rm -rf "$D"
  mkdir -p "$D/constant" "$D/0"
  cp -r "$SRC/system" "$D/system"
  cp -r "$SRC/constant/triSurface" "$D/constant/triSurface"
  cp -r "$SRC/constant/extendedFeatureEdgeMesh" "$D/constant/extendedFeatureEdgeMesh"
  sed -i -E "s/level   \([0-9]+ [0-9]+\);/level   (7 $L);/; s/level   [0-9]+;/level   $L;/; s/maxLocalCells +[0-9]+;/maxLocalCells       2000000;/; s/maxGlobalCells +[0-9]+;/maxGlobalCells      3000000;/" "$D/system/snappyHexMeshDict"
  if [ "${QTIGHT:-0}" = "1" ]; then
    sed -i -E "s/maxBoundarySkewness +[0-9.]+;/maxBoundarySkewness 4;/; s/maxInternalSkewness +[0-9.]+;/maxInternalSkewness 3.5;/; s/maxNonOrtho +75;/maxNonOrtho     70;/" "$D/system/snappyHexMeshDict"
  fi
  if [ -n "${SEDX:-}" ]; then
    sed -i -E "$SEDX" "$D/system/snappyHexMeshDict"
  fi
  {
    echo "dictionary lines:"
    grep -n -E "level|maxLocalCells|maxGlobalCells|maxNonOrtho|Skewness|minVol " "$D/system/snappyHexMeshDict"
  } >> "$OUT"
  (
    cd "$D" || exit 1
    ulimit -v 10000000
    blockMesh </dev/null > log.blockMesh 2>&1 || { echo "blockMesh FAILED"; exit 1; }
    timeout 3600 /usr/bin/time -v snappyHexMesh -overwrite </dev/null > log.snappyHexMesh 2> log.snappy.time
    echo "snappyHexMesh exit code $?"
    checkMesh </dev/null > log.checkMesh 2>&1
    echo "checkMesh exit code $?"
    postProcess -func writeCellCentres -time 0 </dev/null > log.pp1 2>&1
    postProcess -func writeCellVolumes -time 0 </dev/null > log.pp2 2>&1
  ) >> "$OUT" 2>&1
  {
    echo "--- time/memory ---"
    grep -E "Elapsed|Maximum resident" "$D/log.snappy.time"
    echo "--- cells per level (snappy log) ---"
    grep -A 14 "Cells per refinement level" "$D/log.snappyHexMesh" | tail -n 16
    echo "--- snap / errors ---"
    grep -E "Snapped mesh|FOAM FATAL|Not refining|exceed|bad_alloc|Finished meshing" "$D/log.snappyHexMesh" "$D/log.snappy.time"
    echo "--- snappy final check ---"
    grep -A 12 "Checking final mesh" "$D/log.snappyHexMesh"
    echo "--- checkMesh ---"
    grep -E "cells:|Max aspect|Min volume|Mesh non-orthogonality|severely|Max skewness|Error|Failed|Mesh OK" "$D/log.checkMesh"
    echo "--- cell size at feature points ---"
    $PY "$ROOT/scripts/check_te_cells.py" "$D" 0 "$D/constant/extendedFeatureEdgeMesh/wing_edgeMesh.obj" </dev/null
  } >> "$OUT" 2>&1
  rm -f "$D/0/C" "$D/0/Ccx" "$D/0/Ccy" "$D/0/Ccz" "$D/0/Vc"
  echo "=== surface max level $L tag '$TAG' : end $(date '+%F %T') ===" >> "$OUT"
done
echo "SWEEP DONE $TAG" >> "$OUT"
