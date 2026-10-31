#!/usr/bin/python3
"""Frame check: do structural skin nodes lie on the CFD wing STL?

usage: frame_check_surface.py <job.frd> <wing.stl> [xy_radius_m = 0.15] [z_tol_m = 0.04]

For each wing node (x >= 20 m, z >= 2 m, every 3rd node) the STL vertices within xy_radius are found
and the minimum vertical distance min|z_node - z_vertex| is taken.  Nodes with no neighbour are outside the STL footprint.
The test is repeated with the deck shifted by dx in x.
"""
import sys

import numpy as np
from scipy.spatial import cKDTree

frd, stl = sys.argv[1], sys.argv[2]
R = float(sys.argv[3]) if len(sys.argv) > 3 else 0.15
TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 0.04

coords, block = [], None
with open(frd) as f:
    for line in f:
        if line.startswith("    2C"):
            block = "coord"
        elif line.startswith(" -3"):
            block = None
        elif line.startswith(" -1") and block == "coord":
            coords.append([float(line[13 + 12 * i:25 + 12 * i]) for i in range(3)])
X = np.array(coords)
X = X[(X[:, 0] >= 20.0) & (X[:, 2] >= 2.0) & (X[:, 1] >= 0.0)][::3]

rec = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
with open(stl, "rb") as f:
    f.read(80)
    nt = int(np.frombuffer(f.read(4), dtype="<u4")[0])
    t = np.fromfile(f, dtype=rec, count=nt)
V = np.unique(np.round(t["v"].reshape(-1, 3).astype(np.float64), 5), axis=0)
V = V[V[:, 1] >= -0.2]
tree = cKDTree(V[:, :2])
print(f"deck nodes tested {len(X)}   STL vertices (y>=-0.2) {len(V)}   radius {R} m   tol {TOL} m")


def test(dx):
    P = X.copy()
    P[:, 0] += dx
    nb = tree.query_ball_point(P[:, :2], r=R)
    dz = np.full(len(P), np.nan)
    for k, idx in enumerate(nb):
        if idx:
            dz[k] = np.abs(V[idx, 2] - P[k, 2]).min()
    return dz


dz0 = test(0.0)
has = ~np.isnan(dz0)
print(f"\ndx = 0: no STL neighbour {np.sum(~has)} ({100 * np.mean(~has):.1f} %)   "
      f"within tol {100 * np.mean(dz0[has] < TOL):.1f} % of those with neighbour   "
      f"median {np.median(dz0[has]):.4f} m   p90 {np.percentile(dz0[has], 90):.4f} m")
print("\nspanwise (5 m bins) at dx = 0:  y_bin   n   %no-neighbour   %within-tol(of neighbour)")
for lo in range(0, 30, 5):
    m = (X[:, 1] >= lo) & (X[:, 1] < lo + 5)
    h = m & has
    if m.sum():
        print(f"   {lo:2d}-{lo + 5:<2d}  {m.sum():6d}   {100 * np.mean(~has[m]):6.1f}   "
              f"{100 * np.mean(dz0[h] < TOL) if h.sum() else float('nan'):6.1f}")
print("\nshift test:  dx[m]   %within-tol (all nodes)   median dz (nodes with neighbour)")
for dx in np.arange(-1.5, 1.51, 0.25):
    dz = test(dx)
    h = ~np.isnan(dz)
    print(f"   {dx:+5.2f}      {100 * np.mean(np.nan_to_num(dz, nan=9.9) < TOL):6.1f}              "
          f"{np.median(dz[h]) if h.sum() else float('nan'):.4f}")
