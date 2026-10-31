#!/bin/bash
# Reconstruct a decomposed steady run at one time and export VTK (internal mesh + all patches).
# usage: cfd_postprocess.sh <case_dir> <time>
# Requires the OpenFOAM 11 environment (reconstructPar, foamToVTK on PATH).
CASE=$1
T=$2
if [ -z "$CASE" ] || [ -z "$T" ]; then echo "usage: $0 <case_dir> <time>"; exit 1; fi
cd "$CASE" || { echo "no such case $CASE"; exit 1; }
if ! command -v reconstructPar >/dev/null 2>&1; then echo "ABORT: OpenFOAM environment not loaded"; exit 1; fi
if [ ! -d "processor0/$T" ]; then echo "ABORT: processor0/$T does not exist"; ls processor0; exit 1; fi
reconstructPar -time "$T" </dev/null > log.reconstruct 2>&1 || { echo "reconstructPar FAILED"; tail -n 20 log.reconstruct; exit 1; }
foamToVTK -time "$T" </dev/null > log.foamToVTK 2>&1 || { echo "foamToVTK FAILED"; tail -n 20 log.foamToVTK; exit 1; }
if [ ! -f "VTK/wing/wing_$T.vtk" ]; then echo "ABORT: VTK/wing/wing_$T.vtk missing"; ls VTK VTK/wing; exit 1; fi
echo "reconstruct + VTK OK"
ls "$T" VTK
ls -la "VTK/wing/wing_$T.vtk"
