#!/usr/bin/env python3
"""
Where are the bad skin triangles?  Reads the raw Netgen npz, classifies triangles
(face-descriptor idx -> skin piece by area + centroid), and reports min-angle statistics
by class, the faces that own the worst triangles, and the worst individual triangles.

Run in the conda env:
  python3 scripts/diagnose_skin_quality.py mesh/stage3_skin/raw_s2p5_k8.npz
"""
import json
import sys
import numpy as np
from scipy.spatial import cKDTree

RIBS = np.array([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])
d = np.load(sys.argv[1])
P, tris, tri_idx = d["P"], d["tris"], d["tri_idx"]
pieces = json.load(open("mesh/stage3_skin/stage3_skin_faces.json"))
jc = np.array([[p["cx"], p["cy"], p["cz"]] for p in pieces])
ja = np.array([p["area"] for p in pieces])
jsh = np.array([p["shared"] for p in pieces])
tree = cKDTree(jc)

p0, p1, p2 = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
cen = (p0 + p1 + p2) / 3.0
nid = int(tri_idx.max()) + 1
asum = np.bincount(tri_idx, weights=area, minlength=nid)
csum = np.stack([np.bincount(tri_idx, weights=area * cen[:, k], minlength=nid) for k in range(3)], axis=1)
idx_cls = np.zeros(nid, dtype=int)
idx_piece = {}
for i in np.nonzero(asum > 0)[0]:
    dist, j = tree.query(csum[i] / asum[i])
    if dist < 1e-3 and abs(asum[i] - ja[j]) / ja[j] < 1e-5:
        idx_cls[i] = 1 if jsh[j] else 2
        idx_piece[i] = j
tri_cls = idx_cls[tri_idx]

m = tri_cls > 0
sk, sidx, scls = tris[m], tri_idx[m], tri_cls[m]
q0, q1, q2 = P[sk[:, 0]], P[sk[:, 1]], P[sk[:, 2]]
ea, eb, ec = (np.linalg.norm(q1 - q0, axis=1), np.linalg.norm(q2 - q1, axis=1),
              np.linalg.norm(q0 - q2, axis=1))


def ang(a, b, c):
    return np.degrees(np.arccos(np.clip((b * b + c * c - a * a) / (2 * b * c), -1, 1)))


amin = np.minimum(np.minimum(ang(ea, eb, ec), ang(eb, ea, ec)), ang(ec, ea, eb))
ar = 0.5 * np.linalg.norm(np.cross(q1 - q0, q2 - q0), axis=1)
print(f"Skin triangles {len(sk)} (shared {int((scls == 1).sum())}, free {int((scls == 2).sum())})")
for thr in (10.0, 5.0, 1.0):
    b = amin < thr
    print(f"  min angle < {thr:>4g} deg: total {int(b.sum()):5d}  "
          f"(shared {int((b & (scls == 1)).sum())}, free {int((b & (scls == 2)).sum())})")

bad = amin < 5.0
print(f"\nFaces owning triangles with min angle < 5 deg (top 10 by count):")
ids, cnt = np.unique(sidx[bad], return_counts=True)
for k in np.argsort(cnt)[::-1][:10]:
    i = ids[k]
    sel = (tri_idx == i)
    nodes = np.unique(tris[sel])
    ext = P[nodes].max(axis=0) - P[nodes].min(axis=0)
    c = csum[i] / asum[i]
    j = idx_piece[i]
    print(f"  face idx {i:4d}: {cnt[k]:4d} bad of {int(sel.sum()):5d} tris | "
          f"{'shared' if jsh[j] else 'free  '} | piece area {ja[j]:9.2f} mm^2 | "
          f"extent dx {ext[0]:6.2f} dy {ext[1]:6.2f} dz {ext[2]:6.2f} | centroid "
          f"({c[0]:7.2f},{c[1]:6.2f},{c[2]:7.2f})")

print("\nWorst 10 triangles (min angle, class, face idx, edge lengths, centroid):")
for k in np.argsort(amin)[:10]:
    c = (q0[k] + q1[k] + q2[k]) / 3.0
    print(f"  {amin[k]:6.3f} deg  {'shared' if scls[k] == 1 else 'free  '}  idx {sidx[k]:4d}  "
          f"edges {ea[k]:.3f}/{eb[k]:.3f}/{ec[k]:.3f}  area {ar[k]:.4f}  "
          f"centroid ({c[0]:.2f},{c[1]:.2f},{c[2]:.2f})")

zc = (q0[:, 2] + q1[:, 2] + q2[:, 2]) / 3.0
dz = np.abs(zc[:, None] - RIBS[None, :]).min(axis=1)
print("\nDistance (in z) of triangles with min angle < 5 deg to the nearest rib station:")
for lo, hi in ((0, 1), (1, 5), (5, 20), (20, 1e9)):
    sel = bad & (dz >= lo) & (dz < hi)
    print(f"  {lo:>4g} .. {hi:<8g} mm: {int(sel.sum())}")
