#!/bin/bash
# Layer-stage test: add boundary layers to an existing snapped snappyHexMesh mesh.
# usage: [SEDX='sed -E expressions'] snappy_layers_test.sh <source_case> <tag> <nLayers> <firstThickness_m> <expansion>
# Copies simulation/cfd_v1/<source_case> to simulation/cfd_v1/<source_case>_<tag>, switches
# snappyHexMeshDict to addLayers only (absolute thicknesses), applies optional extra sed edits (SEDX),
# runs snappyHexMesh and checkMesh, then measures layer coverage and surface fidelity.
# Summary is appended to results/cfd_snappy_layers_<tag>.txt.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
if [ "$#" -ne 5 ]; then echo "usage: $0 <source_case> <tag> <nLayers> <firstThickness_m> <expansion>"; exit 1; fi
SRCNAME=$1
TAG=$2
NL=$3
T1=$4
ER=$5
SRC=$ROOT/simulation/cfd_v1/$SRCNAME
D=$ROOT/simulation/cfd_v1/${SRCNAME}_$TAG
OUT=$ROOT/results/cfd_snappy_layers_$TAG.txt
PY=/usr/bin/python3
PATCHES='(inlet outlet symmetry farfieldY bottom top)'
if [ ! -d "$SRC/constant/polyMesh" ]; then echo "source mesh not found: $SRC"; exit 1; fi
mkdir -p "$ROOT/results"
echo "=== layers on $SRCNAME : nLayers $NL first $T1 m expansion $ER SEDX='${SEDX:-}' : start $(date '+%F %T') ===" >> "$OUT"
rm -rf "$D"
mkdir -p "$D/0"
cp -r "$SRC/constant" "$D/constant"
cp -r "$SRC/system" "$D/system"
sed -i -E "s/castellatedMesh +true;/castellatedMesh false;/; s/snap +true;/snap            false;/; s/addLayers +false;/addLayers       true;/; s/relativeSizes +true;/relativeSizes       false;/; s/layers +\{\}/layers { wing { nSurfaceLayers $NL; } }/; s/expansionRatio +[0-9.]+;/expansionRatio      $ER;/; s/finalLayerThickness +[0-9.]+;/firstLayerThickness  $T1;/; s/minThickness +[0-9.]+;/minThickness        0.0001;/" "$D/system/snappyHexMeshDict"
if [ -n "${SEDX:-}" ]; then
  sed -i -E "$SEDX" "$D/system/snappyHexMeshDict"
fi
{
  echo "dictionary lines:"
  grep -n -E "castellatedMesh|snap  |addLayers |relativeSizes|layers |expansionRatio|LayerThickness|minThickness|minVol |maxFaceThicknessRatio|maxThicknessToMedialRatio|minMedialAxisAngle|nSmooth|nLayerIter|featureAngle" "$D/system/snappyHexMeshDict"
} >> "$OUT"
(
  cd "$D" || exit 1
  ulimit -v 10000000
  timeout 5400 /usr/bin/time -v snappyHexMesh -overwrite </dev/null > log.snappyHexMesh 2> log.snappy.time
  echo "snappyHexMesh exit code $?"
  checkMesh </dev/null > log.checkMesh 2>&1
  echo "checkMesh exit code $?"
  postProcess -func writeCellVolumes -time 0 </dev/null > log.pp 2>&1
  foamToVTK -ascii -constant -noInternal -excludePatches "$PATCHES" </dev/null > log.foamToVTK 2>&1
) >> "$OUT" 2>&1
{
  echo "--- time/memory ---"
  grep -E "Elapsed|Maximum resident" "$D/log.snappy.time"
  echo "--- layer table (achieved) ---"
  grep -n -B1 -A3 "overall thickness" "$D/log.snappyHexMesh" | tail -n 8
  echo "--- layer mesh / errors ---"
  grep -E "Layer mesh|FOAM FATAL|Finished meshing|Extruding" "$D/log.snappyHexMesh" | tail -n 6
  echo "--- checkMesh ---"
  grep -E "cells:|Max aspect|Min volume|Mesh non-orthogonality|severely|Max skewness|Error|Failed|Mesh OK" "$D/log.checkMesh"
  echo "--- layer coverage ---"
  $PY "$ROOT/scripts/check_layer_coverage.py" "$D" "$D/VTK/wing/wing_0.vtk" </dev/null
  echo "--- surface fidelity ---"
  $PY "$ROOT/scripts/check_surface_fidelity.py" "$D/VTK/wing/wing_0.vtk" "$D/constant/triSurface/wing.stl" </dev/null
} >> "$OUT" 2>&1
echo "=== layers $TAG : end $(date '+%F %T') ===" >> "$OUT"
echo "LAYERS DONE $TAG" >> "$OUT"
