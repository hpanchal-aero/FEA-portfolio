#!/usr/bin/env python3
"""Locate feature-edge length inside a spanwise window and list the concave-looking
clusters. Edges are grouped by 0.5 m spanwise bins; for each bin the length and the
x-range and z-range of edge midpoints are printed. Also prints the chordwise split
(x below / above the local mid-chord is approximated by a global x median of the bin).

usage: check_features_detail.py <file.obj> [ymin=0] [ymax=3.5]
"""
import sys

import numpy as np

v, e = [], []
for ln in open(sys.argv[1]):
    t = ln.split()
    if not t:
        continue
    if t[0] == "v":
        v.append([float(x) for x in t[1:4]])
    elif t[0] == "l":
        idx = [int(x) - 1 for x in t[1:]]
        e.extend(zip(idx[:-1], idx[1:]))
v = np.array(v)
e = np.array(e)
ymin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
ymax = float(sys.argv[3]) if len(sys.argv) > 3 else 3.5
a, b = v[e[:, 0]], v[e[:, 1]]
ln_ = np.linalg.norm(b - a, axis=1)
mid = 0.5 * (a + b)
d = b - a
print(f"window y {ymin}..{ymax} m (y >= 0 side only)")
print(" y-bin [m]     edges   length[m]   x_mid min..max [m]      z_mid min..max [m]   mean |dy|/len")
edges = np.arange(ymin, ymax + 1e-9, 0.5)
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (mid[:, 1] >= lo) & (mid[:, 1] < hi)
    if not m.any():
        print(f" {lo:5.2f}..{hi:5.2f}  {0:7d}")
        continue
    dy = np.abs(d[m, 1]) / ln_[m]
    print(f" {lo:5.2f}..{hi:5.2f}  {int(m.sum()):7d}  {ln_[m].sum():9.3f}   {mid[m,0].min():7.3f}..{mid[m,0].max():7.3f}   {mid[m,2].min():6.3f}..{mid[m,2].max():6.3f}   {dy.mean():6.3f}")
m = (mid[:, 1] >= ymin) & (mid[:, 1] < ymax)
xs = mid[m, 0]
print("x histogram of edge-midpoints in window (2 m bins):")
h, xe = np.histogram(xs, bins=np.arange(np.floor(xs.min()), np.ceil(xs.max()) + 2, 2.0), weights=ln_[m])
for k, val in enumerate(h):
    if val > 0:
        print(f"   x {xe[k]:6.1f}..{xe[k+1]:6.1f}  length {val:8.3f} m")
