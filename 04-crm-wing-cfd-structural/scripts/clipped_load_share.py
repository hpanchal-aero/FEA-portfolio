#!/usr/bin/env python3
"""Share of the wing load carried by pressure-clipped faces (Euler runs, decision D18).

Usage: clipped_load_share.py <time> <alpha_deg> [case_dir]

Reads VTK/wing/wing_<time>.vtk (foamToVTK output) and recomputes the pressure
force on the wing: F = sum (p - pRef) * Sf, with Sf the face area vector of the
patch as OpenFOAM defines it (pointing out of the fluid domain, into the wing).
Compares Cl and Cd with the solver's forceCoeffs.dat at the same time, then
reports the contribution of faces whose pressure sits at the limitPressure
floor, and the Cl change if their pressure were lower than the floor.
"""
import glob
import os
import sys
import numpy as np
import pyvista as pv

time = sys.argv[1]
alpha = float(sys.argv[2])
case = sys.argv[3] if len(sys.argv) > 3 else "simulation/cfd_v1/euler"

P_INF, Q_INF, AREF = 21662.7, 10955.9, 191.870
P_FLOOR = 2166.27
a = np.radians(alpha)
lift = np.array([-np.sin(a), 0.0, np.cos(a)])
drag = np.array([np.cos(a), 0.0, np.sin(a)])

w = pv.read(f"{case}/VTK/wing/wing_{time}.vtk")
w = w.compute_normals(cell_normals=True, point_normals=False,
                      consistent_normals=False, auto_orient_normals=False,
                      split_vertices=False)
w = w.compute_cell_sizes(length=False, area=True, volume=False)
p = np.asarray(w.cell_data["p"])
A = np.asarray(w.cell_data["Area"])
Sf = np.asarray(w.cell_data["Normals"]) * A[:, None]
force = (p - P_INF)[:, None] * Sf
norm = AREF * Q_INF

print(f"case {case}, time {time}, alpha {alpha} deg")
print(f"wing faces {len(p)}, total area {A.sum():.3f} m2 (STL half-wing area 424.491 m2)")

cl = force.sum(axis=0) @ lift / norm
cd = force.sum(axis=0) @ drag / norm
print(f"from surface pressure: Cl {cl:.5f}  Cd {cd:.5f}")

ref = None
for f in sorted(glob.glob(f"{case}/postProcessing/forceCoeffsCompressible/*/forceCoeffs.dat")):
    for line in open(f):
        s = line.split()
        if line.startswith("#") or len(s) != 6:
            continue
        try:
            if float(s[0]) == float(time):
                ref = (float(s[3]), float(s[2]))
        except ValueError:
            pass
if ref:
    print(f"from forceCoeffs.dat : Cl {ref[0]:.5f}  Cd {ref[1]:.5f}"
          f"   (difference Cl {cl - ref[0]:+.5f}, Cd {cd - ref[1]:+.5f})")
else:
    print("no forceCoeffs.dat row found for this time")

mask = p <= P_FLOOR * 1.0005
print(f"\nclipped faces (p <= 1.0005 x floor): {int(mask.sum())} "
      f"({100 * mask.mean():.2f}% of faces, {100 * A[mask].sum() / A.sum():.2f}% of area)")
dcl_c = force[mask].sum(axis=0) @ lift / norm
dcd_c = force[mask].sum(axis=0) @ drag / norm
print(f"their contribution: Cl {dcl_c:+.5f} ({100 * dcl_c / cl:+.2f}% of total), "
      f"Cd {dcd_c:+.5f} ({100 * dcd_c / cd:+.2f}% of total)")
for label, p_new in (("0.5 x floor", 0.5 * P_FLOOR), ("0 Pa", 0.0)):
    dF = (p_new - P_FLOOR) * Sf[mask].sum(axis=0)
    print(f"if clipped pressure were {label}: Cl change {dF @ lift / norm:+.5f} "
          f"({100 * (dF @ lift / norm) / cl:+.2f}% of total), "
          f"Cd change {dF @ drag / norm:+.5f}")
