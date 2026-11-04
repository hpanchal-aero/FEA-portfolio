#!/bin/bash
# Project 05: rerun the early verification cases from the tracked decks in verification/decks/ and compare
# the new result files with the stored ones.
#   bash scripts/run_early_verification.sh
# Cases: V1, V2, V3a (mid-line pin, a recorded FAIL, exit status 1 expected), V3b, V3c, V4, V6 (V6a and V6b).
# Stored copies of the result files are kept in the printed /tmp folder before anything is overwritten.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CCX=/home/harsh/miniconda3/envs/ccx223/bin/ccx
PY=/usr/bin/python3
cd "$ROOT" || exit 2

bash scripts/build_meshes.sh build > /dev/null || { echo "mesh build FAILED"; exit 1; }

BEFORE=$(mktemp -d /tmp/p05_early_before.XXXXXX)
echo "stored copies in $BEFORE"
RES="v1_conduction v2_free_expansion v3a_midline_support v3b_full_edge v3c_membrane v4_gradient v6_buckling v6a_thermal_preload"
for r in $RES; do [ -f "results/$r.txt" ] && cp "results/$r.txt" "$BEFORE/"; done
[ -f simulation/v4_gradient/temperature.inp ] && cp simulation/v4_gradient/temperature.inp "$BEFORE/temperature.inp"
bad=0

solve() {   # NAME
  local n="$1"
  mkdir -p "simulation/$n"
  cp "verification/decks/$n.inp" "simulation/$n/$n.inp"
  ( cd "simulation/$n" && "$CCX" "$n" > ccx.log 2>&1 ) || { echo "ccx FAILED for $n (see simulation/$n/ccx.log)"; exit 1; }
  grep -q "Job finished" "simulation/$n/ccx.log" || { echo "no 'Job finished' in simulation/$n/ccx.log"; exit 1; }
  echo "solved    $n"
}
check() {   # EXPECTED_STATUS LABEL COMMAND...
  local want="$1" label="$2"; shift 2
  "$@" > "$BEFORE/$label.out" 2>&1
  local got=$?
  if [ "$got" -eq "$want" ]; then
    echo "checked   $label (exit status $got, as expected)"
  else
    echo "UNEXPECTED exit status $got (expected $want) for $label"; tail -n 6 "$BEFORE/$label.out"; bad=1
  fi
}

M=mesh/v_n40_z2.inp
solve v1_conduction
check 0 v1  $PY scripts/check_v1.py simulation/v1_conduction/v1_conduction.dat $M 120 20
solve v2_free_expansion
check 0 v2  $PY scripts/check_v2.py simulation/v2_free_expansion/v2_free_expansion.dat $M 120 20
solve v3_restrained
check 1 v3a $PY scripts/check_v3.py simulation/v3_restrained/v3_restrained.dat $M 22 20 v3a_midline_support.txt
solve v3b_full_edge
check 0 v3b $PY scripts/check_v3.py simulation/v3b_full_edge/v3b_full_edge.dat $M 22 20 v3b_full_edge.txt
solve v3c_membrane
check 0 v3c $PY scripts/check_v3.py simulation/v3c_membrane/v3c_membrane.dat $M 22 20 v3c_membrane.txt
mkdir -p simulation/v4_gradient
$PY scripts/make_temp_linear.py $M 20 22 simulation/v4_gradient/temperature.inp > /dev/null || { echo "make_temp_linear FAILED"; exit 1; }
if [ -f "$BEFORE/temperature.inp" ]; then
  if cmp -s "$BEFORE/temperature.inp" simulation/v4_gradient/temperature.inp; then echo "IDENTICAL temperature.inp"; else echo "DIFFERS   temperature.inp"; bad=1; fi
fi
solve v4_gradient
check 0 v4  $PY scripts/check_v4.py simulation/v4_gradient/v4_gradient.dat $M 20 22 20
solve v6_buckling
solve v6_diag_prestress
solve v6b_prescribed
check 0 v6  $PY scripts/check_v6.py

for r in $RES; do
  if [ ! -f "results/$r.txt" ]; then echo "MISSING   $r"; bad=1
  elif [ ! -f "$BEFORE/$r.txt" ]; then echo "NEW       $r (no stored copy)"
  elif cmp -s "$BEFORE/$r.txt" "results/$r.txt"; then echo "IDENTICAL $r"
  else echo "DIFFERS   $r (stored copy: $BEFORE/$r.txt)"; bad=1
  fi
done
if [ "$bad" -eq 0 ]; then echo "EARLY VERIFICATION: REPRODUCED"; else echo "EARLY VERIFICATION: CHECK ABOVE"; exit 1; fi

