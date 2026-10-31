#!/usr/bin/python3
"""Surface-smoothness check: turning angle of the upper-surface profile at spanwise stations.

usage: check_surface_turning.py "label1=surface1" ["label2=surface2" ...]
Each surface is a VTK/STL wing surface.  At y = 4, 5, 6, 7, 8 m the plane y = const cuts the surface;
the upper envelope z(x) is built on a 0.1 m grid and the turning angle
    turn(x) = atan(slope after x) - atan(slope before x)   (slopes over 0.4 m windows)
is evaluated.  Printed per station: median |turn| and the maximum |turn| over 8%..92% chord, and the
turning at x = 29.9 m (the suction-streak location).  For y = 6 m the signed turning profile from
x = 28.9 to 31.1 m is printed for direct comparison between surfaces.
"""
import sys

import numpy as np
import pyvista as pv

STATIONS = (4.0, 5.0, 6.0, 7.0, 8.0)
DX = 0.1
NW = 4  # window = NW * DX = 0.4 m
X0 = 29.9


def turning(poly, y):
    p = np.asarray(poly.slice(normal=(0.0, 1.0, 0.0), origin=(0.0, y, 0.0)).points)
    xa, xb = p[:, 0].min(), p[:, 0].max()
    edges = np.arange(xa, xb + DX, DX)
    idx = np.digitize(p[:, 0], edges) - 1
    xs, zs = [], []
    for i in range(len(edges) - 1):
        k = np.where(idx == i)[0]
        if k.size:
            j = k[np.argmax(p[k, 2])]
            xs.append(p[j, 0])
            zs.append(p[j, 2])
    xg = np.arange(xa + DX, xb - DX, DX)
    zg = np.interp(xg, xs, zs)
    n = NW
    sb = (zg[n:-n] - zg[:-2 * n]) / (n * DX)
    sa = (zg[2 * n:] - zg[n:-n]) / (n * DX)
    return xg[n:-n], np.degrees(np.arctan(sa) - np.arctan(sb)), xa, xb


surfaces = [a.split("=", 1) for a in sys.argv[1:]]
data = {}
for label, path in surfaces:
    poly = pv.read(path)
    if not isinstance(poly, pv.PolyData):
        poly = poly.extract_surface()
    print(f"\n=== {label}: {path}  ({poly.n_cells} cells)")
    for y in STATIONS:
        x, t, xa, xb = turning(poly, y)
        c = xb - xa
        mask = (x >= xa + 0.08 * c) & (x <= xb - 0.08 * c)
        med = np.median(np.abs(t[mask]))
        im = int(np.argmax(np.abs(t) * mask))
        t0 = float(np.interp(X0, x, t)) if x[0] <= X0 <= x[-1] else float("nan")
        print(f"  y {y:3.0f}  chord {c:6.2f}  median|turn| {med:6.3f} deg   max|turn| {abs(t[im]):6.3f} deg at x {x[im]:6.2f}   turn at x=29.9: {t0:+7.3f} deg")
        if y == 6.0:
            data[label] = (x, t)

xs = np.arange(28.9, 31.11, 0.2)
print("\nsigned turning (deg) at y = 6 m, x =", " ".join(f"{v:5.1f}" for v in xs))
for label, (x, t) in data.items():
    print(f"  {label:20s}", " ".join(f"{np.interp(v, x, t):+5.2f}" for v in xs))
