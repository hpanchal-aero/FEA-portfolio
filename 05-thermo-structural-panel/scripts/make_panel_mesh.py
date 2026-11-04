#!/usr/bin/python3
"""Structured hexahedral mesh of the flat Ti panel (Project 05), written as a CalculiX include file.

Usage:
  /usr/bin/python3 scripts/make_panel_mesh.py TAG [--L 0.5] [--W 0.5] [--h 0.002]
        [--nx 40] [--ny 40] [--nz 2] [--etype C3D20R]

Geometry: box x in [0,L], y in [0,W], z in [0,h] (SI, metres); z is the thickness direction.
Elements: 20-node serendipity hexahedra from Gmsh (type 17), renumbered to the CalculiX/Abaqus
C3D20 node order. The conversion is verified (corner orientation, mid-edge node positions,
total volume) before anything is written.
Output: mesh/TAG.inp (nodes, elements, node sets) and mesh/TAG_info.txt (counts and checks).

Node sets (located by coordinates, tolerance 1e-9 m):
  NALL (all nodes), TOP (z=h), BOTTOM (z=0),
  FRAME_X0, FRAME_X1, FRAME_Y0, FRAME_Y1: edge nodes on the mid-thickness plane (z=h/2).
      The frame support pins only the mid-plane, so the edges can rotate (simple support);
      fixing the whole edge face would also fix the rotation and clamp the edge.
  P_A (0,0,h/2), P_B (L,0,h/2), P_C (0,W,h/2), P_D (L,W,h/2), P_CEN (L/2,W/2,h/2),
  COL_CEN: nodes on the vertical line through the plate centre (sorted by z).
"""
import argparse
import os
import sys

import numpy as np
import gmsh

# CalculiX/Abaqus C3D20 node i corresponds to Gmsh local node CCX_FROM_GMSH[i].
CCX_FROM_GMSH = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15]
# Corner-node pair that each CalculiX mid-edge node (local index 8..19) lies between.
CCX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
TOL = 1e-9


def build_mesh(L, W, h, nx, ny, nz):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("panel")
        gmsh.model.occ.addBox(0.0, 0.0, 0.0, L, W, h)
        gmsh.model.occ.synchronize()
        for dim, tag in gmsh.model.getEntities(1):
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, tag)
            ext = np.array([x1 - x0, y1 - y0, z1 - z0])
            k = int(np.argmax(ext))
            others = np.delete(ext, k)
            if others.max() > 1e-3 * ext[k]:
                raise SystemExit("edge %d is not axis-aligned: extents %s" % (tag, ext))
            gmsh.model.mesh.setTransfiniteCurve(tag, (nx, ny, nz)[k] + 1)
        for dim, tag in gmsh.model.getEntities(2):
            gmsh.model.mesh.setTransfiniteSurface(tag)
            gmsh.model.mesh.setRecombine(2, tag)
        for dim, tag in gmsh.model.getEntities(3):
            gmsh.model.mesh.setTransfiniteVolume(tag)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        etypes, etags, enodes = gmsh.model.mesh.getElements(3)
        if list(etypes) != [17]:
            raise SystemExit("expected only 20-node hexahedra (Gmsh type 17), got %s" % list(etypes))
        etags = np.array(etags[0], dtype=np.int64)
        conn = np.array(enodes[0], dtype=np.int64).reshape(-1, 20)
        ntags, ncoord, _ = gmsh.model.mesh.getNodes()
        ntags = np.array(ntags, dtype=np.int64)
        ncoord = np.array(ncoord, dtype=float).reshape(-1, 3)
        quality = np.array(gmsh.model.mesh.getElementQualities(list(etags), "minSJ"))
        version = gmsh.__version__
    finally:
        gmsh.finalize()
    return ntags, ncoord, etags, conn, quality, version


def write_ids(f, ids):
    for i in range(0, len(ids), 16):
        f.write(", ".join(str(int(v)) for v in ids[i:i + 16]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tag")
    ap.add_argument("--L", type=float, default=0.5)
    ap.add_argument("--W", type=float, default=0.5)
    ap.add_argument("--h", type=float, default=0.002)
    ap.add_argument("--nx", type=int, default=40)
    ap.add_argument("--ny", type=int, default=40)
    ap.add_argument("--nz", type=int, default=2)
    ap.add_argument("--etype", default="C3D20R")
    a = ap.parse_args()
    L, W, h = a.L, a.W, a.h

    ntags, ncoord, etags, conn, quality, gver = build_mesh(L, W, h, a.nx, a.ny, a.nz)

    # contiguous node ids 1..N in Gmsh-tag order, elements 1..E in Gmsh-tag order
    order = np.argsort(ntags)
    ntags_s = ntags[order]
    X = ncoord[order]
    idx = np.searchsorted(ntags_s, conn)
    if not np.all(ntags_s[idx] == conn):
        raise SystemExit("element refers to a node that is not in the node list")
    if np.unique(idx).size != X.shape[0]:
        raise SystemExit("mesh has nodes that no element uses")
    eorder = np.argsort(etags)
    conn_ccx = idx[:, CCX_FROM_GMSH][eorder]
    quality = quality[eorder]
    nE, nN = conn_ccx.shape[0], X.shape[0]

    # --- checks on the converted connectivity ---
    P = X[conn_ccx]  # (E, 20, 3)
    v1, v3, v4 = P[:, 1] - P[:, 0], P[:, 3] - P[:, 0], P[:, 4] - P[:, 0]
    triple = np.einsum("ij,ij->i", np.cross(v1, v3), v4)
    mid_err = 0.0
    for k, (p, q) in enumerate(CCX_EDGES):
        mid_err = max(mid_err, float(np.abs(P[:, 8 + k] - 0.5 * (P[:, p] + P[:, q])).max()))
    ext = P[:, :8].max(axis=1) - P[:, :8].min(axis=1)
    vol = float(ext.prod(axis=1).sum())
    lens = np.stack([np.linalg.norm(P[:, p] - P[:, q], axis=1) for p, q in CCX_EDGES], axis=1)
    aspect = float((lens.max(axis=1) / lens.min(axis=1)).max())
    vol_exact = L * W * h
    vol_err = abs(vol - vol_exact) / vol_exact
    if triple.min() <= 0.0:
        raise SystemExit("CHECK FAILED: non-positive corner orientation, min triple product %g" % triple.min())
    if mid_err > 1e-12:
        raise SystemExit("CHECK FAILED: mid-edge node not at the edge midpoint, max error %g m" % mid_err)
    if vol_err > 1e-12:
        raise SystemExit("CHECK FAILED: volume %g differs from L*W*h = %g" % (vol, vol_exact))

    # --- node sets ---
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    isx0, isx1 = np.abs(x) < TOL, np.abs(x - L) < TOL
    isy0, isy1 = np.abs(y) < TOL, np.abs(y - W) < TOL
    mid = np.abs(z - h / 2) < TOL
    cx, cy = np.abs(x - L / 2) < TOL, np.abs(y - W / 2) < TOL

    def sel(mask):
        return np.nonzero(mask)[0] + 1

    col = sel(cx & cy)
    col = col[np.argsort(X[col - 1, 2])]
    sets = {
        "TOP": sel(np.abs(z - h) < TOL),
        "BOTTOM": sel(np.abs(z) < TOL),
        "FRAME_X0": sel(isx0 & mid),
        "FRAME_X1": sel(isx1 & mid),
        "FRAME_Y0": sel(isy0 & mid),
        "FRAME_Y1": sel(isy1 & mid),
        "FACE_X0": sel(isx0),
        "FACE_X1": sel(isx1),
        "FACE_Y0": sel(isy0),
        "FACE_Y1": sel(isy1),
        "P_A": sel(isx0 & isy0 & mid),
        "P_B": sel(isx1 & isy0 & mid),
        "P_C": sel(isx0 & isy1 & mid),
        "P_D": sel(isx1 & isy1 & mid),
        "P_CEN": sel(cx & cy & mid),
        "COL_CEN": col,
    }
    for name, ids in sets.items():
        if len(ids) == 0:
            raise SystemExit("node set %s is empty" % name)
        if name.startswith("P_") and len(ids) != 1:
            raise SystemExit("node set %s has %d nodes, expected 1" % (name, len(ids)))

    # --- write the include file ---
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "mesh"), exist_ok=True)
    inp = os.path.join(root, "mesh", a.tag + ".inp")
    info = os.path.join(root, "mesh", a.tag + "_info.txt")
    with open(inp, "w") as f:
        f.write("** Project 05 panel mesh %s: L=%g W=%g h=%g m, nx=%d ny=%d nz=%d, Gmsh API %s\n"
                % (a.tag, L, W, h, a.nx, a.ny, a.nz, gver))
        f.write("** generated by scripts/make_panel_mesh.py, SI units (m)\n")
        f.write("*NODE, NSET=NALL\n")
        for i in range(nN):
            f.write("%d, %.12e, %.12e, %.12e\n" % (i + 1, X[i, 0], X[i, 1], X[i, 2]))
        f.write("*ELEMENT, TYPE=%s, ELSET=PANEL\n" % a.etype)
        for e in range(nE):
            nodes = [int(v) + 1 for v in conn_ccx[e]]
            f.write("%d, %s,\n" % (e + 1, ", ".join(str(v) for v in nodes[:15])))
            f.write(", ".join(str(v) for v in nodes[15:]) + "\n")
        for name, ids in sets.items():
            f.write("*NSET, NSET=%s\n" % name)
            write_ids(f, ids)

    lines = [
        "tag %s" % a.tag,
        "Gmsh Python API %s" % gver,
        "geometry L=%g W=%g h=%g m" % (L, W, h),
        "divisions nx=%d ny=%d nz=%d" % (a.nx, a.ny, a.nz),
        "element type written %s (Gmsh type 17, 20 nodes)" % a.etype,
        "nodes %d, elements %d" % (nN, nE),
        "min corner triple product %.6e m3, max %.6e m3" % (triple.min(), triple.max()),
        "max mid-edge node midpoint error %.3e m" % mid_err,
        "element volume sum %.9e m3, exact %.9e m3, relative error %.2e" % (vol, vol_exact, vol_err),
        "max element edge aspect ratio %.3f" % aspect,
        "Gmsh minSJ (scaled Jacobian) min %.6f" % quality.min(),
    ]
    for name, ids in sets.items():
        lines.append("set %-9s %6d nodes" % (name, len(ids)))
    text = "\n".join(lines) + "\n"
    with open(info, "w") as f:
        f.write(text)
    sys.stdout.write(text)
    print("written", inp)
    print("written", info)


if __name__ == "__main__":
    main()
