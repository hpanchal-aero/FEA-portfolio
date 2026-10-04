#!/usr/bin/env python3
"""
Weld-and-verify on the Netgen probe mesh (solid + free skin faces).

1. Classify every surface triangle by matching its face-descriptor index to the
   skin-piece list (stage3_skin_faces.json) by area + centroid:
     cls 0 = solid boundary, not skin (rib flat faces, end faces)
     cls 1 = shared skin (skin piece lying on the solid boundary)
     cls 2 = free skin
2. Weld nodes coincident within 1e-6 mm (scipy KD-tree + connected components);
   the surviving ID is the tet node's ID whenever the group contains one.
3. Verify the skin shell mesh (cls 1 + 2) after welding:
     - no degenerate triangles
     - edge-sharing histogram; boundary edges (count 1) must lie only on the
       z = -0.75 and z = 1000.75 end planes (no T-junctions / gaps)
     - how many skin nodes are solid (tet) nodes
Writes <npz>_classified.npz (P, tris, tri_cls, tets, cmap) for the viewer.

Run inside the conda env:
  python3 scripts/weld_verify_skin_probe.py mesh/stage3_skin/probe_maxh12.npz
"""
import json
import sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

Z_LO, Z_HI = -0.75, 1000.75
WELD = 1e-6

path = sys.argv[1]
d = np.load(path)
P, tris, tri_idx, tets, inb = d["P"], d["tris"], d["tri_idx"], d["tets"], d["in_tet_bnd"]
pieces = json.load(open("mesh/stage3_skin/stage3_skin_faces.json"))

# ---------- 1. classify ----------
p0, p1, p2 = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
cen = (p0 + p1 + p2) / 3.0
nid = int(tri_idx.max()) + 1
asum = np.bincount(tri_idx, weights=area, minlength=nid)
csum = np.stack([np.bincount(tri_idx, weights=area * cen[:, k], minlength=nid)
                 for k in range(3)], axis=1)
jc = np.array([[p["cx"], p["cy"], p["cz"]] for p in pieces])
ja = np.array([p["area"] for p in pieces])
jshared = np.array([p["shared"] for p in pieces])
tree = cKDTree(jc)

idx_cls = np.zeros(nid, dtype=int)
matched_j = []
for i in np.nonzero(asum > 0)[0]:
    c = csum[i] / asum[i]
    dist, j = tree.query(c)
    if dist < 1e-3 and abs(asum[i] - ja[j]) / ja[j] < 1e-5:
        idx_cls[i] = 1 if jshared[j] else 2
        matched_j.append(j)
n_desc = int((asum > 0).sum())
print(f"Face descriptors with triangles: {n_desc}")
print(f"Matched to skin pieces: {len(matched_j)} (unique {len(set(matched_j))}); "
      f"shared {int((idx_cls == 1).sum())}, free {int((idx_cls == 2).sum())}, "
      f"non-skin {int(((asum > 0) & (idx_cls == 0)).sum())}")
print(f"JSON pieces: {len(pieces)} (shared {int(jshared.sum())}, free {int((~jshared).sum())})")
tri_cls = idx_cls[tri_idx]
print(f"Triangles: skin-shared {int((tri_cls == 1).sum())}, skin-free {int((tri_cls == 2).sum())}, "
      f"non-skin {int((tri_cls == 0).sum())}")
bad1 = int(((tri_cls == 1) & ~inb).sum())
bad2 = int(((tri_cls == 2) & inb).sum())
print(f"Consistency: shared-skin tris not on tet boundary: {bad1}; "
      f"free-skin tris on tet boundary: {bad2}  (both must be 0)")

# ---------- 2. weld ----------
used = np.unique(np.concatenate([tets.ravel(), tris.ravel()]))
t2 = cKDTree(P[used])
pairs = t2.query_pairs(WELD, output_type="ndarray")
n = len(used)
g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
ncomp, comp = connected_components(g, directed=False)
is_tet = np.isin(used, np.unique(tets))
BIG = 10 ** 9
key = np.where(is_tet, used, used + BIG).astype(np.int64)
mn = np.full(ncomp, 2 ** 62, dtype=np.int64)
np.minimum.at(mn, comp, key)
cmap = np.zeros(int(used.max()) + 1, dtype=np.int64)
cmap[used] = mn[comp] % BIG
sizes = np.bincount(comp)
print(f"\nWeld: {n} node IDs -> {ncomp} distinct nodes (merged groups: {int((sizes > 1).sum())}, "
      f"max group size {sizes.max()})")

# ---------- 3. verify skin shell mesh ----------
skin = tris[tri_cls > 0]
cs = cmap[skin]
degen = int(((cs[:, 0] == cs[:, 1]) | (cs[:, 1] == cs[:, 2]) | (cs[:, 0] == cs[:, 2])).sum())
print(f"\nSkin triangles: {len(cs)}; degenerate after weld: {degen}")
e = np.concatenate([np.sort(cs[:, [0, 1]], axis=1),
                    np.sort(cs[:, [1, 2]], axis=1),
                    np.sort(cs[:, [0, 2]], axis=1)])
ue, ecnt = np.unique(e, axis=0, return_counts=True)
print(f"Skin edge-sharing histogram: {dict(zip(*np.unique(ecnt, return_counts=True)))}")
bnd = ue[ecnt == 1]
zb = P[bnd][:, :, 2]
on_lo = np.all(np.abs(zb - Z_LO) < 1e-3, axis=1)
on_hi = np.all(np.abs(zb - Z_HI) < 1e-3, axis=1)
stray = bnd[~(on_lo | on_hi)]
print(f"Boundary edges: {len(bnd)} (root end {int(on_lo.sum())}, tip end {int(on_hi.sum())}, "
      f"elsewhere {len(stray)})  <- 'elsewhere' must be 0")
for a, b in stray[:5]:
    print("   stray:", np.round(P[a], 4), "->", np.round(P[b], 4))
sk_nodes = np.unique(cs)
n_solid = int(np.isin(sk_nodes, np.unique(tets)).sum())
print(f"Skin nodes: {len(sk_nodes)}; of which solid (tet) nodes: {n_solid}")

out = path.replace(".npz", "_classified.npz")
np.savez(out, P=P, tris=tris, tri_cls=tri_cls, tets=tets, cmap=cmap)
print(f"\nWrote {out}")
