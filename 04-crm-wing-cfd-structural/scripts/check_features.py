#!/usr/bin/env python3
"""Summarise an OpenFOAM feature-edge OBJ file (v and l lines): edges, total length, bounding box, length by spanwise range.

usage: check_features.py <file.obj>
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
if len(e) == 0:
    sys.exit("ABORT: no line elements in the OBJ file")
a, b = v[e[:, 0]], v[e[:, 1]]
ln_ = np.linalg.norm(b - a, axis=1)
mid = 0.5 * (a + b)
print(f"{sys.argv[1]}: {len(v)} points, {len(e)} edges, total length {ln_.sum():.2f} m")
print("bounding box x %.3f..%.3f  y %.3f..%.3f  z %.3f..%.3f" % (v[:, 0].min(), v[:, 0].max(), v[:, 1].min(), v[:, 1].max(), v[:, 2].min(), v[:, 2].max()))
print("edge length [mm]: min %.2f, median %.2f, max %.2f" % tuple(1000 * x for x in (ln_.min(), np.median(ln_), ln_.max())))
for lo, hi in ((-30, -10.9), (-10.9, -3.05), (-3.05, 3.05), (3.05, 10.9), (10.9, 30)):
    m = (mid[:, 1] >= lo) & (mid[:, 1] < hi)
    print(f"  |y| range {lo:6.2f}..{hi:6.2f}: {int(m.sum()):6d} edges, length {ln_[m].sum():8.2f} m")
print(f"half domain (y >= 0) length: {ln_[mid[:, 1] >= 0].sum():.2f} m")
