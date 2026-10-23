#!/usr/bin/env python3
"""Build the wing-alone surface (decision D11, option A) from the uCRM-9 'wbt' IGES geometry.

usage: build_wing_geometry.py <wbt.igs> [sewing_tolerance_m = 0.0015] [output_prefix = geometry/wing_alone]

Model units of the IGES file are kept (metres x 25.4); metric scaling is applied when the surface mesh is exported.
Phase A: import with sewing, map faces to their original tags (bounding box vs results/wbt_geometry_survey.txt), delete all
         faces except the wing faces 5 26 11 12 14 16 18 20 22 24 27, close the inboard trailing-edge slit with a ruled strip,
         write <prefix>_stageA.brep.
Phase B: re-import stage A with sewing (sews the strip), report faces, free edges, neighbours, areas (CAD kernel, 50 mm surface
         mesh, table), write <prefix>.brep.
"""
import os
import re
import sys
from pathlib import Path

import gmsh
import numpy as np

igs = sys.argv[1]
tol_m = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0015
prefix = sys.argv[3] if len(sys.argv) > 3 else "geometry/wing_alone"
table = Path(__file__).resolve().parent.parent / "results" / "wbt_geometry_survey.txt"
S = 1.0 / 25.4
WING = [5, 26, 11, 12, 14, 16, 18, 20, 22, 24, 27]

num = r"(-?\d+\.\d+)"
pat = re.compile(rf"^\s*(\d+)\s+BSpline surface\s+{num}\s+{num}\s+{num}\s+{num}\s+{num}\s+{num}\s+{num}\s+(\d+)\s*$")
ref = {}
for ln in table.read_text().splitlines():
    m = pat.match(ln)
    if m:
        ref[int(m.group(1))] = ([float(m.group(k)) for k in range(2, 8)], float(m.group(8)))
if sorted(ref) != list(range(1, 31)):
    sys.exit(f"ABORT: face table in {table} has {len(ref)} rows, expected 30")


def set_sew():
    for k, v in (("Geometry.OCCSewFaces", 1), ("Geometry.OCCFixDegenerated", 0), ("Geometry.OCCFixSmallEdges", 0),
                 ("Geometry.OCCFixSmallFaces", 0), ("Geometry.OCCMakeSolids", 0), ("Geometry.Tolerance", tol_m / S)):
        gmsh.option.setNumber(k, v)


def map_faces():
    """sewn face tag -> original tag (or None when no original face has a matching bounding box)"""
    out = {}
    for _, t in gmsh.model.getEntities(2):
        b = [v * S for v in gmsh.model.getBoundingBox(2, t)]
        bb = [b[0], b[3], b[1], b[4], b[2], b[5]]
        ranked = sorted((max(abs(x - y) for x, y in zip(bb, ref[o][0])), o) for o in ref)
        out[t] = ranked[0][1] if ranked[0][0] < 0.01 and ranked[1][0] > 3 * ranked[0][0] + 1e-4 else None
    return out


def ends(c):
    lo, hi = gmsh.model.getParametrizationBounds(1, c)
    return np.array(gmsh.model.getValue(1, c, [lo[0], hi[0]])).reshape(2, 3)


# ---------------------------------------------------------------- phase A
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
set_sew()
gmsh.model.add("a")
gmsh.model.occ.importShapes(igs, highestDimOnly=True)
gmsh.model.occ.synchronize()
old = map_faces()
if len(old) != 30 or None in old.values() or len(set(old.values())) != 30:
    sys.exit("ABORT: face mapping not unique in phase A")
gmsh.model.occ.remove([(2, t) for t, o in old.items() if o not in WING], recursive=True)
gmsh.model.occ.synchronize()
left = sorted(old[t] for _, t in gmsh.model.getEntities(2))
if left != sorted(WING):
    sys.exit(f"ABORT: faces left after removal {left}")
print(f"phase A: sewn at {tol_m * 1000:.2f} mm, kept wing faces {left}, curves {len(gmsh.model.getEntities(1))}")

slit = []
for _, c in gmsh.model.getEntities(1):
    up, _ = gmsh.model.getAdjacencies(1, c)
    if len(up) != 1:
        continue
    bb = [v * S for v in gmsh.model.getBoundingBox(1, c)]
    if bb[0] > 36.0 and (bb[4] - bb[1]) > 3.0:
        slit.append(c)
if len(slit) != 2:
    sys.exit(f"ABORT: expected 2 slit curves, found {slit}")
cA, cB = slit
pa, pb = ends(cA), ends(cB)
same = np.linalg.norm(pa[0] - pb[0]) + np.linalg.norm(pa[1] - pb[1])
flip = np.linalg.norm(pa[0] - pb[1]) + np.linalg.norm(pa[1] - pb[0])
print(f"slit curves {cA}, {cB}; end-point mismatch same/flipped {same * S:.4f} / {flip * S:.4f} m -> {'same' if same <= flip else 'reversed'} order")
wA = gmsh.model.occ.addWire([cA])
wB = gmsh.model.occ.addWire([cB if same <= flip else -cB])
strip = gmsh.model.occ.addThruSections([wA, wB], makeSolid=False, makeRuled=True)
gmsh.model.occ.synchronize()
sf = [t for d, t in strip if d == 2]
print(f"strip faces {sf}, area {sum(gmsh.model.occ.getMass(2, t) for t in sf) * S * S:.4f} m2")
os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
gmsh.write(prefix + "_stageA.brep")
gmsh.clear()

# ---------------------------------------------------------------- phase B
set_sew()
gmsh.model.add("b")
gmsh.model.occ.importShapes(prefix + "_stageA.brep", highestDimOnly=True)
gmsh.model.occ.synchronize()
faces = sorted(t for _, t in gmsh.model.getEntities(2))
curves = [t for _, t in gmsh.model.getEntities(1)]
omap = map_faces()
lab = {t: (str(omap[t]) if omap[t] is not None else "STRIP?") for t in faces}
print(f"\nphase B: faces {len(faces)}, curves {len(curves)}")

adj, free, multi = {}, [], []
for c in curves:
    up, _ = gmsh.model.getAdjacencies(1, c)
    adj[c] = sorted(lab[abs(u)] for u in up)
    if len(up) == 1:
        free.append(c)
    elif len(up) > 2:
        multi.append(c)
print(f"free edges {len(free)}, edges shared by 3 or more faces {len(multi)}")
for c in free:
    bb = [v * S for v in gmsh.model.getBoundingBox(1, c)]
    print(f"  free curve {c:4d} face {adj[c][0]:>6s} length {gmsh.model.occ.getMass(1, c) * S:8.4f} m  "
          f"x {bb[0]:.3f}..{bb[3]:.3f}  y {bb[1]:.4f}..{bb[4]:.4f}  z {bb[2]:.3f}..{bb[5]:.3f}")

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05 / S)
gmsh.model.mesh.generate(2)
nt, co, _ = gmsh.model.mesh.getNodes()
nt = np.asarray(nt, dtype=np.int64)
xyz = np.zeros((nt.max() + 1, 3))
xyz[nt] = np.asarray(co).reshape(-1, 3)


def mesh_area(t):
    tot = 0.0
    types, _, ntags = gmsh.model.mesh.getElements(2, t)
    for ty, nn in zip(types, ntags):
        if ty == 2:
            tri = xyz[np.asarray(nn, dtype=np.int64).reshape(-1, 3)]
            tot += 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum()
    return tot * S * S


print("\nface(orig) | CAD-kernel area | mesh area (50 mm) | table area [m2] | triangles | neighbours")
tot_cad = tot_mesh = 0.0
for t in faces:
    a_cad = gmsh.model.occ.getMass(2, t) * S * S
    a_mesh = mesh_area(t)
    ntri = sum(len(e) for ty, e in zip(*gmsh.model.mesh.getElements(2, t)[:2]) if ty == 2)
    tot_cad += a_cad
    tot_mesh += a_mesh
    nb = sorted({g for c in curves if lab[t] in adj[c] for g in adj[c] if g != lab[t]})
    tab = f"{ref[omap[t]][1]:9.3f}" if omap[t] is not None else "      n/a"
    print(f"{lab[t]:>6s}     {a_cad:15.4f} {a_mesh:17.4f} {tab:>16s} {ntri:10d}   {nb}")
print(f"total area: CAD kernel {tot_cad:.3f} m2, mesh {tot_mesh:.3f} m2")
b = [v * S for v in gmsh.model.getBoundingBox(-1, -1)]
print(f"bounding box [m] x {b[0]:.3f}..{b[3]:.3f}  y {b[1]:.4f}..{b[4]:.3f}  z {b[2]:.3f}..{b[5]:.3f}")
gmsh.write(prefix + ".brep")
print(f"wrote {prefix}.brep ({os.path.getsize(prefix + '.brep') / 1e6:.2f} MB, model units = metres x 25.4)")
gmsh.finalize()
