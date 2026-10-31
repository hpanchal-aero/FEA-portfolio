#!/bin/bash
# Two-stage Euler run at one angle of attack for the CL search (decision D12).
# usage: cfd_run_alpha.sh <alpha_deg> <tag> [n_first_order=200] [n_tvd=500] [n_procs=8]
# Stage 1: first-order upwind from a uniform start, n_first_order iterations.
# Stage 2: TVD (limitedLinear) restart for n_tvd more iterations.
# The case is copied from simulation/cfd_v1/euler (template) to simulation/cfd_v1/runs/<tag>.
set -u
alpha=$1
tag=$2
n1=${3:-200}
n2=${4:-500}
np=${5:-8}

root=$(cd "$(dirname "$0")/.." && pwd)
tmpl=$root/simulation/cfd_v1/euler
casedir=$root/simulation/cfd_v1/runs/$tag

if [ -e "$casedir" ]; then echo "ABORT: $casedir already exists"; exit 1; fi
for f in system/fvSchemes.upwind system/fvSchemes.tvd constant/polyMesh/boundary 0/U 0/p 0/T \
         system/controlDict system/fvSolution system/fvConstraints system/decomposeParDict; do
  if [ ! -e "$tmpl/$f" ]; then echo "ABORT: template is missing $f"; exit 1; fi
done

mkdir -p "$casedir"
cp -r "$tmpl/0" "$tmpl/constant" "$tmpl/system" "$casedir/"
cd "$casedir" || exit 1

sed -i "s/^angleOfAttack .*/angleOfAttack   $alpha;/" 0/U
cp system/fvSchemes.upwind system/fvSchemes
sed -i 's/^startFrom .*/startFrom       startTime;/' system/controlDict
sed -i "s/^endTime .*/endTime         $n1;/" system/controlDict
sed -i "s/^numberOfSubdomains .*/numberOfSubdomains $np;/" system/decomposeParDict

echo "=== $(date) case $tag alpha=$alpha n1=$n1 n2=$n2 np=$np"
grep -n "^angleOfAttack" 0/U
grep -n "^startFrom\|^endTime" system/controlDict
grep -n "bounded Gauss" system/fvSchemes

decomposePar -force </dev/null > log.decomposePar 2>&1
if [ $? -ne 0 ]; then echo "ABORT: decomposePar failed"; exit 1; fi

mpirun --bind-to none -np "$np" foamRun -parallel </dev/null > log.stage1 2>&1
if [ $? -ne 0 ]; then echo "ABORT: stage 1 failed"; tail -n 20 log.stage1; exit 1; fi
echo "=== $(date) stage 1 finished"; tail -n 3 log.stage1

cp system/fvSchemes.tvd system/fvSchemes
sed -i 's/^startFrom .*/startFrom       latestTime;/' system/controlDict
sed -i "s/^endTime .*/endTime         $((n1 + n2));/" system/controlDict
grep -n "^startFrom\|^endTime" system/controlDict
grep -n "limitedLinear" system/fvSchemes

mpirun --bind-to none -np "$np" foamRun -parallel </dev/null > log.stage2 2>&1
if [ $? -ne 0 ]; then echo "ABORT: stage 2 failed"; tail -n 20 log.stage2; exit 1; fi
echo "=== $(date) stage 2 finished"; tail -n 3 log.stage2
echo DONE
