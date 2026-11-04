#!/bin/bash
# Project 05: rebuild the meshes that no run script builds (and verify the recipe).
#   bash scripts/build_meshes.sh build    build every mesh below that is missing
#   bash scripts/build_meshes.sh verify   rebuild all into scratch tags (prefix zz_), diff them against the
#                                         existing meshes (tag name ignored), delete the scratch files
# Meshes built by run_case.py, run_v5.py, run_v6_conv.py, run_v6_thick.py and slice_study.py are not listed.
# The diagnostic meshes v8x_* are scratch and are not reproduced.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/usr/bin/python3
S="$ROOT/scripts"
MODE="${1:-}"
case "$MODE" in
  build)  P="" ;;
  verify) P="zz_" ;;
  *) echo "usage: bash scripts/build_meshes.sh build|verify"; exit 2 ;;
esac
cd "$ROOT" || exit 2
if [ "$MODE" = verify ]; then trap 'rm -f mesh/zz_*' EXIT; fi

run() { "$@" > /dev/null || { echo "FAILED: $*"; exit 1; }; }
skip() { [ "$MODE" = build ] && [ -f "$1" ]; }

panel() {   # TAG [make_panel_mesh options]
  local tag="$1"; shift
  if skip "mesh/$tag.inp"; then echo "exists   $tag"; return; fi
  run $PY "$S/make_panel_mesh.py" "$P$tag" "$@"
  echo "built    $P$tag"
}
frame() {   # TAG   (membrane restraint on all four faces)
  local tag="$1"
  if skip "mesh/${tag}_frame.inp"; then echo "exists   ${tag}_frame"; return; fi
  run $PY "$S/frame_constraints.py" "$P$tag" --faces xy
  echo "built    $P${tag}_frame"
}
prestrain() {   # TAG   (mechanical equivalent of a 1 K rise, default --dT)
  local tag="$1"
  if skip "mesh/${tag}_framepre.inp"; then echo "exists   ${tag}_framepre"; return; fi
  run $PY "$S/frame_prestrain.py" "$P$tag"
  echo "built    $P${tag}_framepre"
}
imperf() {  # BASE OUT A_OVER_H
  local base="$1" out="$2" amp="$3"
  if skip "mesh/$out.inp"; then echo "exists   $out"; return; fi
  run $PY "$S/make_imperfect_mesh.py" "$P$base" "$P$out" --a_over_h "$amp"
  echo "built    $P$out"
}

# V1 to V4 and V6 plate (C3D20R, 40 x 40 x 2, make_panel_mesh.py defaults)
panel v_n40_z2
frame v_n40_z2
prestrain v_n40_z2
# base meshes of the V7 imperfect meshes (same commands as run_v6_conv.py)
panel v6c_n20_z2 --nx 20 --ny 20 --nz 2
frame v6c_n20_z2
panel v6c_n40_z2 --nx 40 --ny 40 --nz 2
frame v6c_n40_z2
# V7 imperfect meshes
imperf v6c_n20_z2 v7_n20_z2_a0p01 0.01
imperf v6c_n20_z2 v7_n20_z2_a0p1 0.1
imperf v6c_n20_z2 v7_n20_z2_a0p5 0.5
imperf v6c_n40_z2 v7_n40_z2_a0p1 0.1
# V8a thermal meshes (C3D20, full integration)
panel v8a_n20_z2 --nx 20 --ny 20 --nz 2 --etype C3D20
panel v8a_n40_z2 --nx 40 --ny 40 --nz 2 --etype C3D20

if [ "$MODE" = verify ]; then
  bad=0; n=0
  for f in mesh/zz_*; do
    ref="mesh/${f#mesh/zz_}"
    n=$((n + 1))
    if [ ! -f "$ref" ]; then echo "NO REFERENCE $ref"; bad=1; continue; fi
    if diff <(sed 's/zz_//g' "$f") "$ref" > /dev/null; then echo "IDENTICAL $ref"; else echo "DIFFERS   $ref"; bad=1; fi
  done
  echo "compared $n files"
  if [ "$bad" -eq 0 ]; then echo "VERIFY: PASS"; else echo "VERIFY: FAIL"; exit 1; fi
fi
