#!/usr/bin/env python3
"""Trailing-edge thickness of the wing STL at spanwise stations.

The aft-most x of the surface is a swept line, so at each station y0 it is found by fitting a straight
line x_te(y) to the per-bin maximum x of vertices within +-band of y0 and evaluating it at y0.
Vertical rays are then cast at x = x_te(y0) - offset; thickness = z_max - z_min of the intersections.

usage: te_thickness.py <wing.stl> [offset_m=0.003] [band_m=0.15]
"""
import sys

import numpy as np
import pyvista as pv

stl = pv.read(sys.argv[1])
off = float(sys.argv[2]) if len(sys.argv) > 2 else 0.003
band = float(sys.argv[3]) if len(sys.argv) > 3 else 0.15
P = stl.points
print(f"{sys.argv[1]}: {stl.n_cells} triangles, offset {1000*off:.1f} mm forward of fitted TE, band +-{band} m")
print("   y0 [m]   x_te [m]   n_hits   z_min [m]   z_max [m]   thickness [mm]")
for y0 in (0.5, 3.0, 6.0, 10.0, 15.0, 20.0, 25.0, 28.0, 29.2, 29.35, 29.42):
    m = np.abs(P[:, 1] - y0) < band
    if m.sum() < 50:
        print(f"{y0:9.2f}   too few vertices ({int(m.sum())})")
        continue
    yy, xx = P[m, 1], P[m, 0]
    edges = np.linspace(yy.min(), yy.max(), 21)
    yb, xb = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (yy >= lo) & (yy <= hi)
        if k.any():
            j = np.argmax(xx[k])
            yb.append(yy[k][j])
            xb.append(xx[k][j])
    c = np.polyfit(yb, xb, 1)
    resid = np.abs(np.polyval(c, yb) - np.array(xb)).max()
    xte = np.polyval(c, y0)
    x = xte - off
    pts, _ = stl.ray_trace([x, y0, 1.0], [x, y0, 8.0])
    if len(pts) < 2:
        print(f"{y0:9.2f}  {xte:9.3f}  {len(pts):6d}   (fewer than 2 intersections; fit residual {1000*resid:.1f} mm)")
        continue
    z = np.sort(pts[:, 2])
    print(f"{y0:9.2f}  {xte:9.3f}  {len(pts):6d}  {z[0]:10.4f}  {z[-1]:10.4f}  {1000*(z[-1]-z[0]):12.2f}   (fit residual {1000*resid:.1f} mm)")
