#!/usr/bin/env python3
"""Locate the faces of a checkMesh faceSet (read-only).
usage: locate_bad_faces.py <case> <setName>
Needs: <case>/constant/polyMesh/boundary, <case>/constant/polyMesh/sets/<setName>,
       <case>/postProcessing/checkMesh/constant/<setName>.vtk (from checkMesh -writeSets vtk),
       a wing STL in <case>/constant/triSurface (distance and chord fraction).
"""
import os
import re
import sys
import glob
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree


def parse_boundary(path):
    txt = open(path).read()
    pat = re.compile(r"(\w+)\s*\{[^{}]*?type\s+(\w+);[^{}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", re.S)
    return [(m.group(1), m.group(2), int(m.group(3)), int(m.group(4))) for m in pat.finditer(txt)]


def parse_set(path):
    txt = open(path).read()
    body = txt.split("}", 1)[1]
    body = body[body.index("("):]
    return [int(x) for x in re.findall(r"\d+", body)]


def main(case, name):
    patches = parse_boundary(os.path.join(case, "constant/polyMesh/boundary"))
    n_internal = min(p[3] for p in patches)
    print("patches:", [(p[0], p[1], p[2], p[3]) for p in patches])
    print("internal faces:", n_internal)
    labels = parse_set(os.path.join(case, "constant/polyMesh/sets", name))
    print("set %s: %d faces" % (name, len(labels)))
    where = []
    for L in labels:
        if L < n_internal:
            where.append("internal")
        else:
            w = [p[0] for p in patches if p[3] <= L < p[3] + p[2]]
            where.append(w[0] if w else "?")
    vtk = os.path.join(case, "postProcessing/checkMesh/constant/%s.vtk" % name)
    m = pv.read(vtk)
    C = np.asarray(m.cell_centers().points)
    print("VTK cells:", m.n_cells, " (pairing with labels assumes the same order)")
    stl = [f for f in glob.glob(os.path.join(case, "constant/triSurface/*")) if f.lower().endswith(".stl")]
    wing = None
    for f in stl:
        if "wing" in os.path.basename(f).lower():
            wing = f
    if wing is None and stl:
        wing = stl[0]
    S = pv.read(wing) if wing else None
    print("STL used:", wing)
    tree = cKDTree(np.asarray(S.points)) if S is not None else None
    for k in range(m.n_cells):
        c = C[k]
        lab = labels[k] if k < len(labels) else -1
        wh = where[k] if k < len(where) else "?"
        line = "face %8d  %-10s centre (%.4f, %.4f, %.4f)" % (lab, wh, c[0], c[1], c[2])
        if tree is not None:
            d, _ = tree.query(c)
            line += "  dist to STL vertex %.4f m" % d
            P = np.asarray(S.points)
            sel = np.abs(P[:, 1] - c[1]) < 0.05
            if sel.sum() > 5:
                x0, x1 = P[sel, 0].min(), P[sel, 0].max()
                z0, z1 = P[sel, 2].min(), P[sel, 2].max()
                line += "  local section x %.3f..%.3f (x/c %.3f) z %.3f..%.3f" % (x0, x1, (c[0] - x0) / (x1 - x0), z0, z1)
        print(line)
    print("VTK bounds:", m.bounds)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
