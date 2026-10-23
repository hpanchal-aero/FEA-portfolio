#!/usr/bin/env python3
"""Edge matching between the faces of a candidate wing surface (read-only).

usage: wing_edge_matching.py <wbt.igs> [faces = 5,26,11,12,14,16,18,20,22,24,27,30] [scale = 1/25.4]

Each boundary curve is represented by a polyline (4000 points from Gmsh getValue). For every curve A, 200 points on A are
measured against the polyline of every other listed curve B whose bounding box (padded 50 mm) overlaps that of A, with the
exact point-to-segment distance. d(A->B) = max over the 200 points. The nearest curve minimises d(A->B); the reverse
distance d(B->A) is reported for that partner. Polyline chord error is a few micrometres for the curvatures involved.
Distances in mm.
  < 0.1 mm coincident | < 1 mm sewable at 1 mm | < 5 mm small gap | otherwise free edge or gap
A preliminary list shows which curves the OpenCASCADE projection routine (gmsh getClosestPoint) rejects.
"""
import sys

import gmsh
import numpy as np

igs = sys.argv[1]
faces = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [5, 26, 11, 12, 14, 16, 18, 20, 22, 24, 27, 30]
scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0 / 25.4
NA, NB, PAD = 200, 4000, 0.050

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("edges")
gmsh.model.occ.importShapes(igs, highestDimOnly=True)
gmsh.model.occ.synchronize()


def points(c, n):
    lo, hi = gmsh.model.getParametrizationBounds(1, c)
    u = np.linspace(lo[0], hi[0], n)
    return np.array(gmsh.model.getValue(1, c, list(u))).reshape(-1, 3) * scale      # metres


def seg_dist(P, Q):
    """distance of each point of P (m,3) to the polyline Q (n,3)"""
    a = Q[:-1]
    ab = Q[1:] - a
    ab2 = np.maximum((ab ** 2).sum(1), 1e-30)
    ap = P[:, None, :] - a[None, :, :]
    t = np.clip((ap * ab[None]).sum(2) / ab2[None], 0.0, 1.0)
    proj = a[None] + t[..., None] * ab[None]
    return np.linalg.norm(P[:, None, :] - proj, axis=2).min(axis=1)


curves = {}
for f in faces:
    for _, c in gmsh.model.getBoundary([(2, f)], oriented=False, combined=False):
        c = abs(c)
        poly = points(c, NB)
        curves[(f, c)] = dict(A=points(c, NA), poly=poly, lo=poly.min(axis=0), hi=poly.max(axis=0))

print(f"{igs}: faces {faces}, {len(curves)} boundary curves, scale {scale:.6f}")
bad = []
for (f, c) in curves:
    try:
        gmsh.model.getClosestPoint(1, c, list(points(c, 3)[1] / scale))
    except Exception as e:                      # diagnostic only
        bad.append((f, c, str(e).strip()))
print("curves rejected by the OpenCASCADE projection routine:", bad if bad else "none")

print("\nface curve  length[m]     ymin     ymax | nearest face curve  d(A->B)[mm]  d(B->A)[mm]   class")
summary = {}
for (f, c), A in sorted(curves.items()):
    length = float(np.linalg.norm(np.diff(A["poly"], axis=0), axis=1).sum())
    best = (np.inf, None)
    for key, B in curves.items():
        if key[0] == f:
            continue
        if np.any(A["lo"] - PAD > B["hi"]) or np.any(B["lo"] - PAD > A["hi"]):
            continue
        d = float(seg_dist(A["A"], B["poly"]).max())
        if d < best[0]:
            best = (d, key)
    if best[1] is None:
        print(f"{f:4d} {c:5d} {length:10.3f} {A['lo'][1]:8.3f} {A['hi'][1]:8.3f} |       none nearby                                  FREE")
        summary.setdefault(f, []).append("FREE/GAP")
        continue
    rev = float(seg_dist(curves[best[1]]["A"], A["poly"]).max())
    dm, rm = best[0] * 1000.0, rev * 1000.0
    worst = max(dm, rm)
    cls = "coincident" if worst < 0.1 else "sewable<1mm" if worst < 1.0 else "small gap<5mm" if worst < 5.0 else "FREE/GAP"
    summary.setdefault(f, []).append(cls)
    print(f"{f:4d} {c:5d} {length:10.3f} {A['lo'][1]:8.3f} {A['hi'][1]:8.3f} | {best[1][0]:7d} {best[1][1]:5d} {dm:12.4f} {rm:12.4f}   {cls}")

print("\nper-face classes:")
for f in faces:
    cl = summary.get(f, [])
    print(f"  face {f:2d}: " + ", ".join(f"{k} x{cl.count(k)}" for k in sorted(set(cl))))
gmsh.finalize()
