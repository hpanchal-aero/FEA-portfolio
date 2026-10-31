#!/bin/bash
# Usage: [MESHSRC=<case dir holding constant/polyMesh>] cfd_run_upwind.sh <alpha_deg> <tag> [nIter=1500] [np=8]
# First-order upwind steady Euler run from the template case simulation/cfd_v1/euler
# Case is created in simulation/cfd_v1/runs/<tag>
# MESHSRC (optional): replace the template mesh by constant/polyMesh of that case (mesh-sensitivity study).
ALPHA=$1
TAG=$2
N=${3:-1500}
NP=${4:-8}
MESHSRC=${MESHSRC:-}
BASE=~/FEA-portfolio/04-crm-wing-cfd-structural/simulation/cfd_v1
TEMPLATE=$BASE/euler
CASE=$BASE/runs/$TAG

if [ -z "$ALPHA" ] || [ -z "$TAG" ]; then
  echo "ABORT: usage: cfd_run_upwind.sh <alpha_deg> <tag> [nIter] [np]"; exit 1
fi
if [ -e "$CASE" ]; then
  echo "ABORT: $CASE already exists"; exit 1
fi
if [ -n "$MESHSRC" ] && [ ! -f "$MESHSRC/constant/polyMesh/boundary" ]; then
  echo "ABORT: no mesh found in $MESHSRC/constant/polyMesh"; exit 1
fi

mkdir -p "$BASE/runs"
cp -r "$TEMPLATE" "$CASE" || { echo "ABORT: copy failed"; exit 1; }
cd "$CASE" || { echo "ABORT: cd failed"; exit 1; }
rm -rf processor* postProcessing* VTK dynamicCode log.*
for d in [0-9]*; do [ "$d" != "0" ] && [ -d "$d" ] && rm -rf "$d"; done

if [ -n "$MESHSRC" ]; then
  rm -rf constant/polyMesh
  cp -r "$MESHSRC/constant/polyMesh" constant/polyMesh || { echo "ABORT: mesh copy failed"; exit 1; }
  echo "mesh taken from $MESHSRC"
fi
echo "mesh in case: $(grep -m1 -o 'nCells:[0-9]*' constant/polyMesh/owner)"

# angle of attack
sed -i "s/^angleOfAttack[[:space:]].*/angleOfAttack $ALPHA;/" 0/U
if [ "$(grep -c "^angleOfAttack $ALPHA;" 0/U)" != "1" ]; then
  echo "ABORT: angleOfAttack not set in 0/U"; exit 1
fi

# first-order upwind schemes
cp system/fvSchemes.upwind system/fvSchemes || { echo "ABORT: no fvSchemes.upwind"; exit 1; }
grep -q "upwind" system/fvSchemes || { echo "ABORT: fvSchemes has no upwind"; exit 1; }

# run control: start from 0, end at N
sed -i "s/^startFrom .*/startFrom       startTime;/" system/controlDict
sed -i "s/^startTime .*/startTime       0;/" system/controlDict
sed -i "s/^endTime .*/endTime         $N;/" system/controlDict
grep -E "^(startFrom|startTime|endTime|writeInterval|purgeWrite)" system/controlDict

decomposePar -force > log.decompose 2>&1 || { echo "ABORT: decomposePar failed (see log.decompose)"; exit 1; }
mpirun --bind-to none -np $NP foamRun -parallel > log.run 2>&1 </dev/null
RC=$?
if [ $RC -ne 0 ]; then
  echo "ABORT: foamRun exit code $RC (see log.run)"; exit 1
fi
echo "DONE: $TAG upwind alpha=$ALPHA, $N iterations"
