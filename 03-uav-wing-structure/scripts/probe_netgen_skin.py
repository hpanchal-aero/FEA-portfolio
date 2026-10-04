#!/usr/bin/env python3
"""
Probe: does Netgen mesh a solid + free skin faces conformally?
Loads mesh/stage3_skin/stage3_skin_geom.step, meshes at a coarse size, and reports:
  * shape inventory (solids, faces), Face attribute names available
  * counts of 2D (surface) and 3D (tet) elements, face-descriptor indices
  * conformity: surface triangles that are tet faces (solid boundary) vs not (free skin);
    tet-mesh boundary faces vs solid-boundary surface triangles
Writes mesh/stage3_skin/probe_maxh<maxh>.npz (points, tris, tri face index, tets)
for the viewer. Read-only w.r.t. geometry.

Usage (system Python):  /usr/bin/python3 scripts/probe_netgen_skin.py <maxh_mm>
"""
import sys
import time
import numpy as np
from netgen.occ import OCCGeometry

STEP = "mesh/stage3_skin/stage3_skin_geom.step"


def main():
    maxh = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
    geo = OCCGeometry(STEP)
    solids = list(geo.shape.solids)
    faces = list(geo.shape.faces)
    print(f"Shape inventory: {len(solids)} solid(s), {len(faces)} face(s)")
    f0 = faces[0]
    names = [a for a in dir(f0) if not a.startswith("_")]
    print("Face attributes:", names)

    t0 = time.time()
    try:
        mesh = geo.GenerateMesh(maxh=maxh)
    except Exception as e:
        print(f"GenerateMesh FAILED: {type(e).__name__}: {e}")
        sys.exit(1)
    print(f"GenerateMesh OK (maxh={maxh}) in {time.time() - t0:.1f} s")

    pts = mesh.Points()
    P = np.zeros((len(pts) + 1, 3))
    for i, p in enumerate(pts, start=1):
        P[i] = [p.p[0], p.p[1], p.p[2]]

    tris, tri_idx = [], []
    for el in mesh.Elements2D():
        vs = [v.nr for v in el.vertices]
        tris.append(vs)
        tri_idx.append(el.index)
    tets = []
    for el in mesh.Elements3D():
        tets.append([v.nr for v in el.vertices])
    tri_lens = sorted(set(len(t) for t in tris))
    tet_lens = sorted(set(len(t) for t in tets))
    print(f"Points: {len(pts)}; 2D elements: {len(tris)} (vertex counts {tri_lens}); "
          f"3D elements: {len(tets)} (vertex counts {tet_lens})")
    print(f"Distinct 2D face-descriptor indices: {len(set(tri_idx))}")
    if not tets or not tris:
        print("ERROR: missing 2D or 3D elements"); sys.exit(1)

    tris = np.array(tris)
    tets = np.array(tets)
    tri_idx = np.array(tri_idx)

    # tet faces: count occurrences of sorted node triples
    combos = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    tf = np.concatenate([np.sort(tets[:, list(c)], axis=1) for c in combos])
    uniq, cnt = np.unique(tf, axis=0, return_counts=True)
    print(f"Tet face-sharing histogram: {dict(zip(*np.unique(cnt, return_counts=True)))}")
    bfaces = {tuple(r) for r in uniq[cnt == 1]}
    print(f"Tet-mesh boundary faces: {len(bfaces)}")

    ts = np.sort(tris, axis=1)
    in_tet_bnd = np.array([tuple(r) in bfaces for r in ts])
    n_solid_bnd = int(in_tet_bnd.sum())
    n_free = int((~in_tet_bnd).sum())
    print(f"Surface triangles coinciding with tet boundary faces: {n_solid_bnd}")
    print(f"Surface triangles NOT on the tet boundary (free skin candidates): {n_free}")
    print(f"Tet boundary faces with NO surface triangle: {len(bfaces) - n_solid_bnd}")

    # free triangle nodes that touch the tet mesh (should be only along contact edges)
    tet_nodes = set(np.unique(tets).tolist())
    free_tris = tris[~in_tet_bnd]
    if len(free_tris):
        nodes_free = np.unique(free_tris)
        touching = [n for n in nodes_free if n in tet_nodes]
        print(f"Free-skin nodes: {len(nodes_free)}; of which also tet nodes "
              f"(contact-curve nodes): {len(touching)}")

    out = f"mesh/stage3_skin/probe_maxh{maxh:g}.npz".replace(".", "p", 0)
    np.savez(out, P=P, tris=tris, tri_idx=tri_idx, tets=tets, in_tet_bnd=in_tet_bnd)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

