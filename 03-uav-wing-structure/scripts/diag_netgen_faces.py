#!/usr/bin/env python3
"""
Diagnostics for Netgen surface meshing of the skin-imprint geometry
(mesh/stage3_skin/stage3_skin_geom[_rs<N>].step; RIB_SPACING env selects the suffix like mesh_stage3_skin_raw.py).
Run with the system python, unbuffered:  /usr/bin/python3 -u scripts/diag_netgen_faces.py <mode> ...

  map <raw.npz>                        check how the triangle face-descriptor index relates to the STEP face order
                                       (area and centroid of each index's triangles vs shape.faces[idx-1] / [idx])
  verbose <solid_h> <skin_h> <level>   GenerateMesh with SetMessageImportance(level) and the same per-face maxh as
                                       mesh_stage3_skin_raw.py (run under `timeout`)
  face <k> [<k2> ...]                  centre, area, class and bounding box of shape.faces[k-1] (1-based, like
                                       the verified numbering of `map`; the bbox is LOOSE, ignore it)
  edges <k> [<k2> ...]                 edge end points and chord lengths of those faces
  single <k> <maxh>                    mesh face k ALONE with Netgen (run under `timeout`)
  topo                                 faces per edge (edges used by one face = possible non-conformity)
"""
import json
import os
import sys
import time
import numpy as np
from netgen.occ import OCCGeometry

_RS = float(os.environ.get("RIB_SPACING", "200"))
RS = "" if _RS == 200.0 else "_rs%g" % _RS
STEP = "mesh/stage3_skin/stage3_skin_geom" + RS + ".step"
PIECES = "mesh/stage3_skin/stage3_skin_faces" + RS + ".json"


def xyz(c):
    try:
        return (c.x, c.y, c.z)
    except AttributeError:
        return (c[0], c[1], c[2])


def load():
    geo = OCCGeometry(STEP)
    shape = geo.shape
    faces = list(shape.faces)
    pieces = json.load(open(PIECES))
    jc = np.array([[p["cx"], p["cy"], p["cz"]] for p in pieces])
    ja = np.array([p["area"] for p in pieces])
    jsh = np.array([p["shared"] for p in pieces])
    return shape, faces, jc, ja, jsh


def classify(face, jc, ja, jsh):
    c = np.array(xyz(face.center))
    a = face.mass
    dist = np.linalg.norm(jc - c, axis=1)
    j = int(dist.argmin())
    if dist[j] < 1e-3 and abs(a - ja[j]) / ja[j] < 1e-5:
        return "shared" if jsh[j] else "free"
    return "non-skin"


def mode_map(path):
    d = np.load(path)
    P, tris, idx = d["P"], d["tris"], d["tri_idx"]
    shape, faces, jc, ja, jsh = load()
    print(f"STEP faces {len(faces)}; descriptor indices {idx.min()}..{idx.max()} ({len(np.unique(idx))} distinct)")
    p0, p1, p2 = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    cen = (p0 + p1 + p2) / 3.0
    hits = {"faces[idx-1]": 0, "faces[idx]": 0}
    total = 0
    bad = []
    for k in np.unique(idx):
        m = idx == k
        A = area[m].sum()
        C = (cen[m] * area[m][:, None]).sum(axis=0) / A
        total += 1
        ok_first = False
        for off, name in ((1, "faces[idx-1]"), (0, "faces[idx]")):
            j = int(k) - off
            if 0 <= j < len(faces):
                f = faces[j]
                if np.linalg.norm(np.array(xyz(f.center)) - C) < 0.05 and abs(f.mass - A) / f.mass < 1e-3:
                    hits[name] += 1
                    if off == 1:
                        ok_first = True
        if not ok_first and len(bad) < 5:
            bad.append((int(k), float(A), C.round(3).tolist()))
    print(f"indices tested {total}: matches {hits}")
    if bad:
        print("first indices NOT matching faces[idx-1] (idx, mesh area, mesh centroid):", bad)


def set_sizes(shape, faces, jc, ja, jsh, solid_h, skin_h):
    n = {"free": 0, "shared": 0, "non-skin": 0}
    for f in faces:
        cl = classify(f, jc, ja, jsh)
        f.maxh = skin_h if cl == "free" else solid_h
        n[cl] += 1
    for s in shape.solids:
        s.maxh = solid_h
    return n


def mode_verbose(solid_h, skin_h, level, kw=None):
    kw = kw or {}
    from netgen.meshing import SetMessageImportance
    shape, faces, jc, ja, jsh = load()
    n = set_sizes(shape, faces, jc, ja, jsh, solid_h, skin_h)
    print(f"faces {len(faces)}: {n}; maxh solid {solid_h}, skin {skin_h}; message importance {level}")
    geo = OCCGeometry(shape)
    SetMessageImportance(level)
    t0 = time.time()
    print(f"extra GenerateMesh keywords: {kw}")
    geo.GenerateMesh(maxh=skin_h, **kw)
    print(f"GenerateMesh returned in {time.time() - t0:.1f} s")


def mode_face(ks):
    shape, faces, jc, ja, jsh = load()
    for k in ks:
        f = faces[k - 1]
        c = xyz(f.center)
        line = f"face {k}: centre ({c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f}), area {f.mass:.6f} mm^2, class {classify(f, jc, ja, jsh)}"
        try:
            bb = f.bounding_box
            line += f", bbox {tuple(round(v, 4) for v in xyz(bb[0]))} .. {tuple(round(v, 4) for v in xyz(bb[1]))}"
        except Exception:
            pass
        print(line)


def mode_edges(ks):
    shape, faces, jc, ja, jsh = load()
    for k in ks:
        f = faces[k - 1]
        es = list(f.edges)
        print(f"face {k}: area {f.mass:.6f}, class {classify(f, jc, ja, jsh)}, {len(es)} edges")
        for i, e in enumerate(es):
            c = xyz(e.center)
            print(f"   edge {i}: centre ({c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f}), mass {getattr(e, 'mass', None)}")


def mode_single(k, maxh):
    shape, faces, jc, ja, jsh = load()
    f = faces[k - 1]
    f.maxh = maxh
    geo = OCCGeometry(f)
    t0 = time.time()
    mesh = geo.GenerateMesh(maxh=maxh)
    n2 = sum(1 for _ in mesh.Elements2D())
    print(f"face {k} alone (maxh {maxh}): GenerateMesh returned in {time.time() - t0:.2f} s, {n2} surface triangles")


def mode_topo():
    from collections import Counter
    shape, faces, jc, ja, jsh = load()
    cnt = Counter()
    for f in faces:
        for e in f.edges:
            c = xyz(e.center)
            key = (round(c[0], 3), round(c[1], 3), round(c[2], 3), round(float(getattr(e, "mass", 0.0)), 3))
            cnt[key] += 1
    hist = Counter(cnt.values())
    print(f"faces {len(faces)}; distinct edge keys {len(cnt)}; shape.edges {len(list(shape.edges))}; "
          f"shape.vertices {len(list(shape.vertices))}")
    print("faces per edge key (count: number of edges):", dict(sorted(hist.items())))
    singles = [k for k, v in cnt.items() if v == 1]
    print(f"edge keys used by ONE face: {len(singles)}; first 10: {singles[:10]}")


def parse_kw(items):
    kw = {}
    for it in items:
        k, v = it.split("=", 1)
        if v == "None":
            kw[k] = None
        elif v in ("True", "False"):
            kw[k] = (v == "True")
        else:
            try:
                kw[k] = int(v)
            except ValueError:
                kw[k] = float(v)
    return kw


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    mode = sys.argv[1]
    if mode == "map" and len(sys.argv) == 3:
        mode_map(sys.argv[2])
    elif mode == "verbose" and len(sys.argv) >= 5:
        mode_verbose(float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), parse_kw(sys.argv[5:]))
    elif mode == "face" and len(sys.argv) >= 3:
        mode_face([int(v) for v in sys.argv[2:]])
    elif mode == "edges" and len(sys.argv) >= 3:
        mode_edges([int(v) for v in sys.argv[2:]])
    elif mode == "single" and len(sys.argv) == 4:
        mode_single(int(sys.argv[2]), float(sys.argv[3]))
    elif mode == "topo" and len(sys.argv) == 2:
        mode_topo()
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
