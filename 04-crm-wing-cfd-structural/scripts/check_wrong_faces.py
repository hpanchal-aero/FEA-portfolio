#!/usr/bin/env python3
"""Locate faces of an OpenFOAM face set (converted to VTK) relative to the wing surface and the feature edges.

usage: check_wrong_faces.py <faceset.vtk> <wing.stl> <wing_edgeMesh.obj>
"""
import sys

import numpy as np
import pyvista as pv

fs = pv.read(sys.argv[1])
stl = pv.read(sys.argv[2])
v = []
for ln in open(sys.argv[3]):
    t = ln.split()
    if t and t[0] == "v":
        v.append([float(x) for x in t[1:4]])
v = np.array(v)
if fs.n_cells == 0:
    sys.exit("ABORT: face set VTK has no cells")
sz = fs.compute_cell_sizes(length=False, area=True, volume=False)
cc = fs.cell_centers().points
area = sz.cell_data["Area"]
print(f"{sys.argv[1]}: {fs.n_cells} faces, {fs.n_points} points")
print("face |        x        y        z  |  area [m2]  | dist to wing surface [mm] | dist to feature edge [mm]")
for i, p in enumerate(cc):
    cid, cp = stl.find_closest_cell(p, return_closest_point=True)
    ds = np.linalg.norm(cp - p)
    df = np.linalg.norm(v - p, axis=1).min()
    print(f"{i:4d} | {p[0]:8.3f} {p[1]:8.3f} {p[2]:8.3f} | {area[i]:11.3e} | {1000*ds:25.3f} | {1000*df:25.3f}")
print("bounds of all 4 faces: x %.3f..%.3f  y %.3f..%.3f  z %.3f..%.3f" % (
    fs.bounds[0], fs.bounds[1], fs.bounds[2], fs.bounds[3], fs.bounds[4], fs.bounds[5]))
