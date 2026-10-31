#!/usr/bin/env python3
"""Locate cells clipped by limitPressure / limitTemperature.

Usage: check_floor_cells.py [time] [case_dir]

Reads VTK/euler_<time>.vtk and the wing patch VTK written by foamToVTK.
Reports, for cells at the pressure floor, temperature floor and temperature
ceiling: how many mesh cells, bounding box, spanwise (y) histogram in 2 m bins,
and distance to the nearest wing-surface face centre.
"""
import glob
import sys
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

time = sys.argv[1] if len(sys.argv) > 1 else "600"
case = sys.argv[2] if len(sys.argv) > 2 else "simulation/cfd_v1/euler"

P_FLOOR = 2166.27   # Pa, 0.1 x freestream (limitPressure minFactor)
T_FLOOR = 100.0     # K
T_CEIL = 400.0      # K

mesh = pv.read(sorted(glob.glob(f"{case}/VTK/*_{time}.vtk"))[0])
p = np.asarray(mesh.cell_data["p"])
T = np.asarray(mesh.cell_data["T"])
cid = np.asarray(mesh.cell_data["cellID"])
ctr = mesh.cell_centers().points

wing = pv.read(f"{case}/VTK/wing/wing_{time}.vtk")
tree = cKDTree(wing.cell_centers().points)

print(f"time {time}: {mesh.n_cells} VTK cells, "
      f"{np.unique(cid).size} mesh cells")


def report(name, mask):
    n_mesh = np.unique(cid[mask]).size
    print(f"\n--- {name}: {n_mesh} mesh cells ({int(mask.sum())} VTK cells)")
    if n_mesh == 0:
        return
    c = ctr[mask]
    print("bbox x [%.2f %.2f]  y [%.2f %.2f]  z [%.2f %.2f]" % (
        c[:, 0].min(), c[:, 0].max(), c[:, 1].min(), c[:, 1].max(),
        c[:, 2].min(), c[:, 2].max()))
    d, _ = tree.query(c)
    print("distance to nearest wing face centre [m]: min/25/50/75/max =",
          np.round(np.percentile(d, [0, 25, 50, 75, 100]), 3))
    print("fraction within 0.1 m: %.3f, within 0.5 m: %.3f, beyond 2 m: %.3f" % (
        np.mean(d < 0.1), np.mean(d < 0.5), np.mean(d > 2.0)))
    h, e = np.histogram(c[:, 1], bins=np.arange(0.0, 32.0, 2.0))
    print("y histogram (2 m bins from y = 0):")
    for k in range(h.size):
        if h[k]:
            print(f"   y {e[k]:5.1f}-{e[k+1]:5.1f}: {h[k]}")


report("pressure floor (p <= 1.0005 x floor)", p <= P_FLOOR * 1.0005)
report("temperature floor (T <= 100.05 K)", T <= T_FLOOR + 0.05)
report("temperature ceiling (T >= 399.95 K)", T >= T_CEIL - 0.05)
