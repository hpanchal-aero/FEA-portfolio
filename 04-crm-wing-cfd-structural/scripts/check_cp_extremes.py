#!/usr/bin/python3
"""Diagnose the Cp contour figure: locate Cp extremes, calibrate the upper/lower split,
and print the spanwise trend of upper-surface Cp near the root (y = 0..8 m).

usage: check_cp_extremes.py "label1=wing1.vtk" ["label2=wing2.vtk" ...]
Cp = (p - 21662.7) / 10955.9.  Upper surface = faces whose normal has the same z-sign as the
majority of strongly suctioned faces (Cp < -0.5).
"""
import sys

import numpy as np
import pyvista as pv

PINF, Q = 21662.7, 10955.9

for arg in sys.argv[1:]:
    label, path = arg.split("=", 1)
    m = pv.read(path)
    if not isinstance(m, pv.PolyData):
        m = m.extract_surface()
    cp = (np.asarray(m.cell_data["p"]) - PINF) / Q
    C = np.asarray(m.cell_centers().points)
    nz = np.asarray(m.compute_normals(cell_normals=True, point_normals=False).cell_data["Normals"])[:, 2]
    A = np.asarray(m.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"])
    suc = cp < -0.5
    s = 1.0 if nz[suc].mean() > 0 else -1.0
    up = (s * nz) > 0
    print(f"\n=== {label}")
    print(f"faces {m.n_cells}  upper {up.sum()} ({100 * up.mean():.1f}%)  mean nz of suction faces {nz[suc].mean():+.3f}  total area {A.sum():.1f} m2")
    for thr in (-1.0, -1.1, -1.2):
        k = cp < thr
        if k.any():
            print(f"Cp < {thr:.1f}: {k.sum()} faces, area {A[k].sum():.3f} m2, "
                  f"x {C[k, 0].min():.2f}..{C[k, 0].max():.2f}, y {C[k, 1].min():.2f}..{C[k, 1].max():.2f}, "
                  f"z {C[k, 2].min():.2f}..{C[k, 2].max():.2f}")
        else:
            print(f"Cp < {thr:.1f}: 0 faces")
    order = np.argsort(cp)[:8]
    print("8 lowest-Cp faces (x, y, z, Cp, area):")
    for i in order:
        print(f"  {C[i, 0]:7.2f} {C[i, 1]:6.2f} {C[i, 2]:6.2f}  {cp[i]:+.3f}  {A[i]:.2e}")
    print("upper surface, y-bins (0.5 m): faces, area, area-weighted mean Cp, min Cp")
    for y0 in np.arange(0.0, 8.0, 0.5):
        k = up & (C[:, 1] >= y0) & (C[:, 1] < y0 + 0.5)
        if k.any():
            print(f"  y {y0:4.1f}-{y0 + 0.5:4.1f}: {k.sum():6d}  {A[k].sum():7.3f}  {np.average(cp[k], weights=A[k]):+.4f}  {cp[k].min():+.3f}")
