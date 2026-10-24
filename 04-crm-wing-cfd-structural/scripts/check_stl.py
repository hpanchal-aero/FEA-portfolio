#!/usr/bin/env python3
"""Independent check of a closed triangulated surface (STL) with PyVista and numpy.

usage: check_stl.py <surface.stl>
Reports points, triangles, bounds, boundary and non-manifold edges, area, signed volume (positive = outward normals),
edge-length statistics, triangle minimum angle (numpy, independent of the PyVista version) and degenerate triangles.
"""
import sys

import numpy as np
import pyvista as pv

m = pv.read(sys.argv[1])
print(f"{sys.argv[1]}: {m.n_points} points, {m.n_cells} cells, bounds {np.round(m.bounds, 3).tolist()}")
nb = m.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False).n_cells
nm = m.extract_feature_edges(boundary_edges=False, non_manifold_edges=True, feature_edges=False, manifold_edges=False).n_cells
print(f"boundary (open) edges {nb}, non-manifold edges {nm}")
f = m.faces.reshape(-1, 4)[:, 1:]
p = m.points
v0, v1, v2 = p[f[:, 0]], p[f[:, 1]], p[f[:, 2]]
tri_area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
vol = np.einsum("ij,ij->i", v0, np.cross(v1, v2)).sum() / 6.0
print(f"area {tri_area.sum():.3f} m2, signed volume {vol:.3f} m3 ({'outward' if vol > 0 else 'INWARD'} normals)")
a, b, c = (np.linalg.norm(v1 - v2, axis=1), np.linalg.norm(v2 - v0, axis=1), np.linalg.norm(v0 - v1, axis=1))
e = np.concatenate([a, b, c])
print("edge length [mm]: min %.2f, 1%% %.2f, median %.2f, 99%% %.2f, max %.2f" % tuple(
    1000 * x for x in (e.min(), *np.percentile(e, [1, 50, 99]), e.max())))


def ang(x, y, z):
    return np.degrees(np.arccos(np.clip((y * y + z * z - x * x) / np.maximum(2 * y * z, 1e-30), -1.0, 1.0)))


amin = np.minimum(np.minimum(ang(a, b, c), ang(b, c, a)), ang(c, a, b))
print("triangle minimum angle [deg]: min %.2f, 1%% %.2f, median %.1f" % (amin.min(), np.percentile(amin, 1), np.median(amin)))
print(f"triangles with minimum angle below 10 deg: {int((amin < 10).sum())} ({100 * (amin < 10).mean():.2f} %); "
      f"below 1 deg: {int((amin < 1).sum())}; zero-area (< 1e-12 m2): {int((tri_area < 1e-12).sum())}")
