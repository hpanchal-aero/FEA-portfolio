#!/usr/bin/env bash
# Regenerate the case-dependent figures of Project 04 with explicit case sets (PyVista off-screen, Matplotlib).
# Usage: bash scripts/make_figures.sh
# Needs the CalculiX jobs in simulation/static_v3 and simulation/static_cfd and the CFD VTK files in simulation/cfd_v1/runs.
# The other figures: plot_cp_sections.py <time> <case_dir> <tag>  (tags cfd_cp_sections_*), mesh_study.py / plot_mesh_study.py,
# plot_force_histories.py.
set -eu
cd "$(dirname "$0")/.."
PY=/usr/bin/python3
STC=simulation/static_cfd
RUNS=simulation/cfd_v1/runs
DATASET=simulation/static_v3/coarse_tr_f2/uCRM9_tr_f2
L9=$STC/p0p36_full/uCRM9_cfd_p0p36_full
L10=$STC/L10_p0p36/uCRM9_cfd_L10_p0p36
CORR=$STC/defl_blend_cl0p5_pb/uCRM9_cfd_defl_blend_cl0p5_pb

echo "== displacement contours"
"$PY" scripts/contour_displacement.py figures/contour_displacement_cl0p5.png \
  "Dataset loads=$DATASET" \
  "One-way L9, CFD CL 0.4998=$L9" \
  "Corrected, CL 0.5=$CORR"

echo "== spanwise load, deflection and twist: dataset, one-way L9 and L10"
"$PY" scripts/plot_oneway_results.py figures/oneway_coupling_p0p36.png \
  "dataset loads=$DATASET.inp" \
  "one-way L9=$L9.inp" \
  "one-way L10=$L10.inp"

echo "== spanwise load, deflection and twist: dataset, one-way L9 and corrected CL 0.5"
"$PY" scripts/plot_oneway_results.py figures/coupling_correction_cl0p5.png \
  "dataset loads=$DATASET.inp" \
  "one-way L9=$L9.inp" \
  "corrected, CL 0.5=$CORR.inp"

echo "== surface Cp contours"
"$PY" scripts/contour_cp_surface.py figures/contour_cp_surface_cl0p5.png \
  "One-way L9|CL 0.4998, alpha 0.36 deg=$RUNS/u_p0p36/VTK/wing/wing_1500.vtk" \
  "Corrected|deflected wing, loads interpolated to CL 0.5=$RUNS/u_defl_blend_cl0p5/VTK/wing/wing_blend.vtk"

echo "== Mach slices"
"$PY" scripts/contour_mach_slices.py figures/contour_mach_slices.png \
  "Undeflected, alpha 0.36 deg=$RUNS/u_p0p36/VTK/u_p0p36_1500.vtk" \
  "Deflected, alpha 0.36 deg=$RUNS/u_defl_p0p36/VTK/u_defl_p0p36_1500.vtk" \
  "Deflected, alpha 2.07 deg=$RUNS/u_defl_a2p07/VTK/u_defl_a2p07_1300.vtk"
echo "make_figures done"
