#!/usr/bin/python3
"""Deform a binary wing STL with per-section kinematics (decision D22, correction step).

usage: deform_wing_stl.py <section_kinematics.csv> <in.stl> <out.stl> [scale = 1.0]

Displacement of a vertex (x, y, z), s = sign(y), y' = |y|:
    ux = ux(y'),  uy = s uy(y'),  uz = a(y') + b(y') (x - xbar(y'))
Section values are linearly interpolated in y' with a zero-motion root point at y' = 0 and constant extrapolation
beyond the last section.  The same function is applied to every vertex, so shared vertices stay shared.
Normals are recomputed from the deformed triangles.
"""
import sys

import numpy as np

kin, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
scale = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0

K = np.genfromtxt(kin, delimiter=",", names=True)
Y = np.concatenate([[0.0], K["y_c"]])
UX = np.concatenate([[0.0], K["ux"]])
UY = np.concatenate([[0.0], K["uy"]])
A = np.concatenate([[0.0], K["uz"]])
B = np.concatenate([[0.0], K["b_slope"]])
XB = np.concatenate([[K["xbar"][0]], K["xbar"]])

rec = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
with open(src, "rb") as f:
    f.read(80)
    nt = int(np.frombuffer(f.read(4), dtype="<u4")[0])
    t = np.fromfile(f, dtype=rec, count=nt)
P = t["v"].reshape(-1, 3).astype(np.float64)

y = np.abs(P[:, 1])
s = np.sign(P[:, 1])
D = np.column_stack([np.interp(y, Y, UX), s * np.interp(y, Y, UY),
                     np.interp(y, Y, A) + np.interp(y, Y, B) * (P[:, 0] - np.interp(y, Y, XB))]) * scale
Q = P + D


def tri_geom(V):
    T = V.reshape(-1, 3, 3)
    c = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    area = 0.5 * np.linalg.norm(c, axis=1)
    nrm = c / (2 * area[:, None])
    vol = np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6.0
    return nrm, area, vol


n0, a0, v0 = tri_geom(P)
n1, a1, v1 = tri_geom(Q)
dot = np.einsum("ij,ij->i", n0, n1)
ratio = a1 / a0
mag = np.linalg.norm(D, axis=1)
print(f"triangles {nt}   scale {scale}")
print(f"displacement magnitude: max {mag.max():.4f} m   mean {mag.mean():.4f} m   max uz {D[:, 2].max():.4f}   min uz {D[:, 2].min():.4f}")
print(f"tip-most vertices (|y| > 29.4): mean uz {D[y > 29.4, 2].mean():.4f}   mean ux {D[y > 29.4, 0].mean():.4f}")
print(f"volume {v0:.4f} -> {v1:.4f} m3   change {100 * (v1 - v0) / v0:+.3f} %")
print(f"area   {a0.sum():.4f} -> {a1.sum():.4f} m2   change {100 * (a1.sum() - a0.sum()) / a0.sum():+.3f} %")
print(f"triangle area ratio: min {ratio.min():.4f}   max {ratio.max():.4f}")
print(f"normal rotation: min cos {dot.min():.5f}   triangles with cos < 0.98: {int((dot < 0.98).sum())}   inverted (cos < 0): {int((dot < 0).sum())}")
print(f"bbox before: x {P[:, 0].min():.3f}..{P[:, 0].max():.3f}  y {P[:, 1].min():.3f}..{P[:, 1].max():.3f}  z {P[:, 2].min():.3f}..{P[:, 2].max():.3f}")
print(f"bbox after : x {Q[:, 0].min():.3f}..{Q[:, 0].max():.3f}  y {Q[:, 1].min():.3f}..{Q[:, 1].max():.3f}  z {Q[:, 2].min():.3f}..{Q[:, 2].max():.3f}")

out = np.zeros(nt, dtype=rec)
out["n"] = n1.astype(np.float32)
out["v"] = Q.reshape(-1, 3, 3).astype(np.float32)
with open(dst, "wb") as f:
    f.write(b"Deformed by scripts/deform_wing_stl.py".ljust(80, b"\0"))
    f.write(np.uint32(nt).tobytes())
    out.tofile(f)
print(f"wrote {dst}")
