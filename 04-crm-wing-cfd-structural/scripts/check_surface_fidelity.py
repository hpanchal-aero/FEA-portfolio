#!/usr/bin/env python3
"""Geometric fidelity of a snapped wing patch: area and vertex/face-centre distance to the STL surface.

usage: check_surface_fidelity.py <patch.vtk> <wing.stl>
"""
import sys

import numpy as np
import pyvista as pv

patch = pv.read(sys.argv[1])
stl = pv.read(sys.argv[2])
if patch.n_cells == 0:
    sys.exit("ABORT: patch has no cells")

ta = stl.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"]
tc = stl.cell_centers().points
a_stl_half = float(ta[tc[:, 1] > 0].sum())
a_patch = float(patch.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"].sum())
print(f"patch: {patch.n_cells} faces, {patch.n_points} points, bbox {tuple(round(x, 3) for x in patch.bounds)}")
print(f"area: patch {a_patch:.3f} m2, STL (triangle centres y > 0) {a_stl_half:.3f} m2, ratio {a_patch / a_stl_half:.5f}")


def closest(points):
    try:
        _, cp = stl.find_closest_cell(points, return_closest_point=True)
        return np.asarray(cp)
    except Exception:
        out = np.empty_like(points)
        for i, p in enumerate(points):
            _, out[i] = stl.find_closest_cell(p, return_closest_point=True)
        return out


for label, pts in (("vertices", patch.points), ("face centres", patch.cell_centers().points)):
    d = np.linalg.norm(closest(pts) - pts, axis=1)
    p = np.percentile(d, [50, 90, 99, 99.9])
    imax = int(np.argmax(d))
    print(f"{label}: N={len(d)}  distance to STL [mm]: median {1000*p[0]:.3f}  90% {1000*p[1]:.3f}  99% {1000*p[2]:.3f}  99.9% {1000*p[3]:.3f}  max {1000*d.max():.3f}")
    print(f"   count > 1 mm: {int((d > 0.001).sum())}   > 5 mm: {int((d > 0.005).sum())}   > 10 mm: {int((d > 0.010).sum())}")
    print(f"   max at ({pts[imax,0]:.3f}, {pts[imax,1]:.3f}, {pts[imax,2]:.3f})")
    big = d > 0.005
    if big.any():
        b = pts[big]
        print(f"   points > 5 mm lie in: x {b[:,0].min():.3f}..{b[:,0].max():.3f}  y {b[:,1].min():.3f}..{b[:,1].max():.3f}  z {b[:,2].min():.3f}..{b[:,2].max():.3f}")
