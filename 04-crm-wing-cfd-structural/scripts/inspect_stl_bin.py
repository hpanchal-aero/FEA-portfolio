#!/usr/bin/python3
"""Inspect a binary STL: triangle count, size consistency, bounding box, unique vertices."""
import sys
import numpy as np

path = sys.argv[1]
with open(path, "rb") as f:
    header = f.read(80)
    n = int(np.frombuffer(f.read(4), dtype="<u4")[0])
    rec = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    tri = np.fromfile(f, dtype=rec, count=n)

import os
expected = 84 + 50 * n
actual = os.path.getsize(path)
print(f"file      : {path}")
print(f"header    : {header[:40]!r}")
print(f"triangles : {n}   size expected {expected}  actual {actual}  match {expected == actual}")
v = tri["v"].reshape(-1, 3).astype(np.float64)
lo, hi = v.min(axis=0), v.max(axis=0)
for k, ax in enumerate("xyz"):
    print(f"  {ax}: {lo[k]:10.4f} .. {hi[k]:10.4f}")
u = np.unique(np.round(v, 5), axis=0)
print(f"unique vertices (5 dp): {len(u)}   Euler V-E+F check needs edges, skipped")
