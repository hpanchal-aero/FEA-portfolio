#!/usr/bin/env python3
"""
Raw Netgen mesh of the imprinted skin + solid geometry
(mesh/stage3_skin/stage3_skin_geom.step, from build_stage3_skin_geom.py).

Mesh sizes:
  * Netgen face/solid maxh can only REFINE the global maxh, so the global maxh is the
    (coarser) SKIN size and the solid, its boundary faces and the shared skin faces are
    set to the (finer) SOLID size.
  * free skin faces (matched to stage3_skin_faces.json by centroid + area) keep the skin size.

Writes mesh/stage3_skin/raw_s<solid>_k<skin>.npz (points, tris, tri_idx, tets, in_tet_bnd).
Also prints tet count and edge-length statistics so we can SEE whether the per-face and
per-solid maxh settings took effect (baseline: ~76,140 tets at solid maxh 2.5 for the
spar v2 + ribs solid alone).

Usage (system Python):
  /usr/bin/python3 scripts/mesh_stage3_skin_raw.py <solid_maxh> <skin_maxh>
"""
import json
import os
import sys
import time
import numpy as np
from netgen.occ import OCCGeometry

import variant_tag as _vt
RS_SUFFIX = _vt.configuration(1000.0)[3]
STEP = "mesh/stage3_skin/stage3_skin_geom" + RS_SUFFIX + ".step"
PIECES = "mesh/stage3_skin/stage3_skin_faces" + RS_SUFFIX + ".json"


def xyz(c):
    try:
        return (c.x, c.y, c.z)
    except AttributeError:
        return (c[0], c[1], c[2])


def main():
    if len(sys.argv) not in (3, 4, 5):
        print("Usage: [TE_MODE=edge|vertex|msz] mesh_stage3_skin_raw.py <solid_maxh> <skin_maxh> "
              "[<te_h> [<te_z>]]"); sys.exit(1)
    solid_h, skin_h = float(sys.argv[1]), float(sys.argv[2])
    te_h = float(sys.argv[3]) if len(sys.argv) > 3 else None
    te_z = float(sys.argv[4]) if len(sys.argv) > 4 else 800.0
    te_mode = os.environ.get("TE_MODE", "edge")
    tag = f"raw_s{solid_h:g}_k{skin_h:g}".replace(".", "p") + RS_SUFFIX
    gen_kwargs = {}
    for _it in os.environ.get("NETGEN_KW", "").split():
        _k, _v = _it.split("=", 1)
        gen_kwargs[_k] = None if _v == "None" else ((_v == "True") if _v in ("True", "False") else float(_v))
        tag += "_" + _k + _v.replace(".", "p").lower()
    if gen_kwargs:
        print(f"extra GenerateMesh keywords from NETGEN_KW: {gen_kwargs}")
    xte = 0.7 * (300.0 - 0.15 * te_z)
    if te_h is not None:
        tag += f"_te{te_h:g}_z{te_z:g}".replace(".", "p")
        if te_mode != "edge":
            tag += f"_{te_mode}"
        if te_mode == "msz":
            segs = [(0.0, 0.23), (-3.0, 0.7), (-7.2, -1.22), (-7.2, 1.2), (-12.0, 2.0), (-12.0, -1.0)]
            mszpath = f"mesh/stage3_skin/{tag}.msz"
            with open(mszpath, "w") as fm:
                fm.write("0\n")
                fm.write(f"{len(segs)}\n")
                for dx, y in segs:
                    fm.write(f"{xte + dx:.4f} {y:.4f} {te_z - 2.0:.4f} "
                             f"{xte + dx:.4f} {y:.4f} {te_z + 2.0:.4f} {te_h}\n")
            gen_kwargs["meshsizefilename"] = mszpath
            print(f"local size file {mszpath}: {len(segs)} lines near the TE at z={te_z:g}, h={te_h} mm")

    pieces = json.load(open(PIECES))
    jc = np.array([[p["cx"], p["cy"], p["cz"]] for p in pieces])
    ja = np.array([p["area"] for p in pieces])
    jshared = np.array([p["shared"] for p in pieces])

    geo0 = OCCGeometry(STEP)
    shape = geo0.shape
    faces = list(shape.faces)
    n_free = n_shared = n_non = 0
    for f in faces:
        c = np.array(xyz(f.center))
        a = f.mass
        dist = np.linalg.norm(jc - c, axis=1)
        j = int(dist.argmin())
        if dist[j] < 1e-3 and abs(a - ja[j]) / ja[j] < 1e-5:
            if jshared[j]:
                f.maxh = solid_h; n_shared += 1
            else:
                f.maxh = skin_h; n_free += 1
        else:
            f.maxh = solid_h; n_non += 1
    print(f"Faces: {len(faces)} (free skin {n_free} -> maxh {skin_h}; "
          f"shared skin {n_shared} + non-skin {n_non} -> maxh {solid_h})")
    exp_free = int(np.sum(jshared == 0))
    exp_shared = int(np.sum(jshared != 0))
    exp_non = 72 if RS_SUFFIX == "" else None   # not predicted for other spacings, reported only
    print(f"Expected from the imprint JSON: free {exp_free}, shared {exp_shared}, "
          f"non-skin {exp_non if exp_non is not None else 'not predicted'}")
    if n_free != exp_free or n_shared != exp_shared or (exp_non is not None and n_non != exp_non):
        print("ERROR: face classification does not match the imprint JSON -- STOP."); sys.exit(1)
    for s in shape.solids:
        s.maxh = solid_h

    if te_h is not None and te_mode in ("vertex", "edge"):
        items = shape.vertices if te_mode == "vertex" else shape.edges
        nsel = 0
        for it in items:
            if te_mode == "vertex":
                try:
                    c = xyz(it.p)
                except AttributeError:
                    c = xyz(it.center)
            else:
                c = xyz(it.center)
            if abs(c[2] - te_z) < 1.2 and (xte - 15.0) < c[0] < (xte + 1.0):
                it.maxh = te_h
                nsel += 1
        print(f"{te_mode}.maxh = {te_h} set on {nsel} shape(s) near the TE at z={te_z:g}")
        if nsel == 0:
            print("ERROR: no shapes selected -- STOP."); sys.exit(1)
    geo = OCCGeometry(shape)
    t0 = time.time()
    try:
        mesh = geo.GenerateMesh(maxh=skin_h, **gen_kwargs)
    except Exception as e:
        print(f"GenerateMesh FAILED: {type(e).__name__}: {e}"); sys.exit(1)
    print(f"GenerateMesh OK (global maxh={skin_h}) in {time.time() - t0:.1f} s")

    pts = mesh.Points()
    P = np.zeros((len(pts) + 1, 3))
    for i, p in enumerate(pts, start=1):
        P[i] = [p.p[0], p.p[1], p.p[2]]
    tris, tri_idx = [], []
    for el in mesh.Elements2D():
        tris.append([v.nr for v in el.vertices])
        tri_idx.append(el.index)
    tets = [[v.nr for v in el.vertices] for el in mesh.Elements3D()]
    tris, tri_idx, tets = np.array(tris), np.array(tri_idx), np.array(tets)
    if tris.shape[1] != 3 or tets.shape[1] != 4:
        print("ERROR: unexpected element vertex counts -- STOP."); sys.exit(1)

    combos = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    tf = np.concatenate([np.sort(tets[:, list(c)], axis=1) for c in combos])
    uniq, cnt = np.unique(tf, axis=0, return_counts=True)
    bset = {tuple(r) for r in uniq[cnt == 1]}
    in_bnd = np.array([tuple(r) in bset for r in np.sort(tris, axis=1)])
    print(f"Points {len(pts)}; tets {len(tets)}; surface tris {len(tris)} "
          f"(on tet boundary {int(in_bnd.sum())}, free {int((~in_bnd).sum())}); "
          f"tet face histogram {dict(zip(*np.unique(cnt, return_counts=True)))}; "
          f"tet-boundary faces without a surface tri: {len(bset) - int(in_bnd.sum())}")

    def med_edge(tr):
        p0, p1, p2 = P[tr[:, 0]], P[tr[:, 1]], P[tr[:, 2]]
        e = np.concatenate([np.linalg.norm(p1 - p0, axis=1), np.linalg.norm(p2 - p1, axis=1),
                            np.linalg.norm(p0 - p2, axis=1)])
        return np.median(e), np.percentile(e, 95), e.max()

    for name, m in (("solid-boundary tris", in_bnd), ("free-skin tris", ~in_bnd)):
        med, p95, mx = med_edge(tris[m])
        print(f"  {name}: median edge {med:.2f} mm, 95th pct {p95:.2f}, max {mx:.2f}")

    if te_h is not None:
        c0 = np.array([xte, 0.23, te_z])
        e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [0, 2]]])
        pa, pb = P[e[:, 0]], P[e[:, 1]]
        near = (np.linalg.norm(pa - c0, axis=1) < 6.0) & (np.linalg.norm(pb - c0, axis=1) < 6.0)
        if near.any():
            Ln = np.linalg.norm(pa[near] - pb[near], axis=1)
            print(f"Local refinement check (surface edges within 6 mm of the TE vertex at z={te_z:g}): "
                  f"n={int(near.sum())}, median {np.median(Ln):.3f} mm, min {Ln.min():.3f} mm")
        else:
            print("Local refinement check: no surface edges within 6 mm of the TE vertex")
    out = f"mesh/stage3_skin/{tag}.npz"
    np.savez(out, P=P, tris=tris, tri_idx=tri_idx, tets=tets, in_tet_bnd=in_bnd)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

