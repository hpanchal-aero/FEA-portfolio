#!/usr/bin/env python3
"""Test whether feature edges forward of x_cut lie on the y=0 mirror plane.

usage: check_features_symplane.py <file.obj> [x_cut=30.0]
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
        i = [int(x) - 1 for x in t[1:]]
        e.extend(zip(i[:-1], i[1:]))
v = np.array(v)
e = np.array(e)
xcut = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
a, b = v[e[:, 0]], v[e[:, 1]]
L = np.linalg.norm(b - a, axis=1)
m = 0.5 * (a + b)
s = m[:, 0] < xcut
if not s.any():
    sys.exit("ABORT: no edges forward of x_cut")
ay = np.abs(m[s, 1])
print(f"edges with x_mid < {xcut} m: {int(s.sum())}, length {L[s].sum():.3f} m")
print(f"  |y| of midpoints: min {ay.min():.5f}  max {ay.max():.5f} m")
print(f"  edges with |y| < 0.01 m: {int((ay < 0.01).sum())} of {int(s.sum())}")
print(f"  length with |y| < 0.01 m: {L[s][ay < 0.01].sum():.3f} m")
print(f"  x range {m[s,0].min():.3f}..{m[s,0].max():.3f}   z range {m[s,2].min():.3f}..{m[s,2].max():.3f}")
