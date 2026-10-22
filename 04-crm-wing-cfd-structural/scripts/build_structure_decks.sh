#!/usr/bin/env bash
# Rebuild the uCRM-9 coarse structure-only decks from the NASTRAN BDFs, in the order used for the verified results.
#   1. ucrm_to_ccx.py coarse                 -> simulation/static_v1/coarse/uCRM9_coarse_v1.inp   (shells + aero nodal forces)
#   2. add_masses.py coarse static           -> simulation/static_v3/coarse/uCRM9_coarse_v3.inp   (gravity + point masses)
#      add_masses.py coarse modal --modes 10 -> simulation/modal_v1/coarse/uCRM9_coarse_modal.inp
#   3. variant_truss_mass.py --unit-area     -> static_v3/coarse_tr_f2/uCRM9_tr_f2.inp, modal_v1/coarse_f2/uCRM9_f2.inp
#      (point masses as T3D2 pairs, E = 1e9 Pa, unit area; works around the add_sm_st error of MASS + *EQUATION)
# Usage: bash scripts/build_structure_decks.sh            (build only)
#        RUN=1 bash scripts/build_structure_decks.sh      (build, run CalculiX, post-process)
set -eu
cd "$(dirname "$0")/.."
PY=/usr/bin/python3

"$PY" scripts/ucrm_to_ccx.py coarse </dev/null
"$PY" scripts/add_masses.py coarse static </dev/null
"$PY" scripts/add_masses.py coarse modal --modes 10 </dev/null

mkdir -p simulation/static_v3/coarse_tr_f2 simulation/modal_v1/coarse_f2
"$PY" scripts/variant_truss_mass.py simulation/static_v3/coarse/uCRM9_coarse_v3.inp \
      simulation/static_v3/coarse_tr_f2/uCRM9_tr_f2.inp --unit-area </dev/null
"$PY" scripts/variant_truss_mass.py simulation/modal_v1/coarse/uCRM9_coarse_modal.inp \
      simulation/modal_v1/coarse_f2/uCRM9_f2.inp --unit-area </dev/null

if [ "${RUN:-0}" = "1" ]; then
  FOLLOW=0 bash scripts/run_ccx.sh simulation/static_v3/coarse_tr_f2 uCRM9_tr_f2 </dev/null
  FOLLOW=0 bash scripts/run_ccx.sh simulation/modal_v1/coarse_f2 uCRM9_f2 </dev/null
  "$PY" scripts/post_static.py simulation/static_v3/coarse_tr_f2/uCRM9_tr_f2 2.5239 </dev/null
  "$PY" scripts/check_equilibrium.py simulation/static_v3/coarse_tr_f2/uCRM9_tr_f2 </dev/null
  "$PY" scripts/post_modal.py simulation/modal_v1/coarse_f2/uCRM9_f2.dat 5 5 </dev/null
fi
echo "build_structure_decks: done"
