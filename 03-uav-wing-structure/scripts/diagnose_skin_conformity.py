#!/usr/bin/env python3
"""
Diagnostic: is the free skin mesh spatially coincident with the solid mesh along the
contact curves, despite sharing no node IDs?  Uses the .npz written by probe_netgen_skin.py.

Reports (a) nearest tet-node distance for every free-skin node, and (b) groups of
nodes with identical coordinates (1e-6 mm), classified by member type.

Run inside the conda env (needs scipy):
  python3 scripts/diagnose_skin_conformity.py mesh/stage3_skin/probe_maxh12.npz
"""
import sys
import numpy as np
from scipy.spatial import cKDTree

d = np.load(sys.argv[1])
P, tris, tets, inb = d["P"], d["tris"], d["tets"], d["in_tet_bnd"]

tet_ids = np.unique(tets)
solid_surf_ids = np.unique(tris[inb])
free_ids = np.unique(tris[~inb])
print(f"Unique nodes: tets {len(tet_ids)}, solid-boundary tris {len(solid_surf_ids)}, "
      f"free tris {len(free_ids)}")
print(f"IDs shared between free-tri nodes and tet nodes: "
      f"{len(np.intersect1d(free_ids, tet_ids))}")

tree = cKDTree(P[tet_ids])
dist, _ = tree.query(P[free_ids])
print("\nNearest tet-node distance for each free-skin node (mm):")
for tol in (1e-9, 1e-6, 1e-3, 1e-1, 1.0):
    print(f"  < {tol:g}: {int((dist < tol).sum())}")
print(f"  min {dist.min():.3e}, median {np.median(dist):.3f}, max {dist.max():.3f}")

used = np.union1d(tet_ids, free_ids)
R = np.round(P[used] / 1e-6).astype(np.int64)
_, inv, cnt = np.unique(R, axis=0, return_inverse=True, return_counts=True)
inv = inv.ravel()
print(f"\nNode IDs in use: {len(used)}; distinct coordinate groups (1e-6 mm): {len(cnt)}")
print(f"Group size histogram: {dict(zip(*np.unique(cnt, return_counts=True)))}")
is_tet = np.isin(used, tet_ids)
is_free = np.isin(used, free_ids)
n_mixed = n_free_only_dup = 0
for g in np.nonzero(cnt > 1)[0]:
    m = inv == g
    if is_tet[m].any() and is_free[m].any():
        n_mixed += 1
    elif is_free[m].all():
        n_free_only_dup += 1
print(f"Groups containing both a tet node and a free-skin node: {n_mixed}")
print(f"Duplicate groups made only of free-skin nodes: {n_free_only_dup}")
