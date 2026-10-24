#!/usr/bin/env python3
"""Cell size at the wing feature-edge points of an OpenFOAM mesh.

Reads the ascii volVectorField C (cell centres) and volScalarField V (cell volumes) written by
postProcess, estimates each cell size as V^(1/3) and its refinement level for a 10 m background cell,
then looks up the nearest cell centre of every feature-edge point (KD-tree; approximate, since the
nearest centre is not always the containing cell).

usage: check_te_cells.py <case_dir> <time_dir> <wing_edgeMesh.obj>
"""
import sys

import numpy as np
from scipy.spatial import cKDTree


def read_internal(path, vec):
    lines = open(path).read().split("\n")
    for i, ln in enumerate(lines):
        if ln.startswith("internalField"):
            break
    else:
        sys.exit(f"ABORT: no internalField in {path}")
    if "nonuniform" not in ln:
        sys.exit(f"ABORT: {path} internalField is not 'nonuniform' ascii: {ln}")
    n = int(lines[i + 1])
    if lines[i + 2].strip() != "(" or lines[i + 3 + n].strip() != ")":
        sys.exit(f"ABORT: {path} unexpected list layout (binary or different format?)")
    rows = lines[i + 3:i + 3 + n]
    if vec:
        return np.array([[float(x) for x in r.strip().strip("()").split()] for r in rows])
    return np.array([float(r) for r in rows])


case, t, obj = sys.argv[1], sys.argv[2], sys.argv[3]
C = read_internal(f"{case}/{t}/C", True)
V = read_internal(f"{case}/{t}/Vc", False)
size = np.cbrt(V)
lev = np.rint(np.log2(10.0 / size)).astype(int)
print(f"cells: {len(V)}  total volume {V.sum():.6e} m3")
print("estimated level histogram from V^(1/3) (trimmed cells appear one level too fine):")
for L in range(0, 13):
    k = int((lev == L).sum())
    if k:
        print(f"   level {L:2d}  size {10000/2**L:9.2f} mm  cells {k}")

pts = []
for ln in open(obj):
    tt = ln.split()
    if tt and tt[0] == "v":
        pts.append([float(x) for x in tt[1:4]])
pts = np.array(pts)
pts = pts[pts[:, 1] >= -1e-9]
d, j = cKDTree(C).query(pts)
print(f"feature points with y >= 0: {len(pts)}")
print("nearest-centre distance [mm]: median %.1f  90%% %.1f  max %.1f" % tuple(1000 * np.percentile(d, [50, 90, 100])))
L = lev[j]
print("estimated level of the nearest cell at each feature point:")
for l in sorted(set(L)):
    print(f"   level {l:2d}: {int((L == l).sum()):6d} points ({100.0*(L == l).mean():5.1f} %)")
print("by spanwise range (y):")
for lo, hi in ((0, 3.05), (3.05, 10.9), (10.9, 25.0), (25.0, 29.0), (29.0, 30.0)):
    m = (pts[:, 1] >= lo) & (pts[:, 1] < hi)
    if m.any():
        print(f"   y {lo:5.2f}..{hi:5.2f}: {int(m.sum()):5d} pts  median level {int(np.median(L[m]))}  fraction level >= 10: {100.0*(L[m] >= 10).mean():5.1f} %  median cell size {1000*np.median(size[j][m]):7.2f} mm")
