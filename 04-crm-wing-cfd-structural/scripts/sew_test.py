#!/usr/bin/env python3
"""Test sewing of the uCRM-9 'wbt' faces at import (Geometry.OCCSewFaces), read-only.

usage: sew_test.py <wbt.igs> <tolerance_m> [face_table = results/wbt_geometry_survey.txt]

Sewing renumbers the faces, so every sewn face is mapped back to its original tag by bounding box and area, using the face
table of the survey record. All reporting uses original tags. Tolerance in metres, converted to model units (x 25.4).
If the mapping is not unique for all 30 faces the unmapped faces are printed with their nearest candidates and the script stops.
Otherwise it reports counts, curves by number of adjacent faces, and for the wing faces 5 26 11 12 14 16 18 20 22 24 27 (and the
junction section 30) their neighbouring faces and free edges.
"""
import re
import sys
from pathlib import Path

import gmsh

igs = sys.argv[1]
tol_m = float(sys.argv[2])
table = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).resolve().parent.parent / "results" / "wbt_geometry_survey.txt"
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

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
for k, v in (("Geometry.OCCSewFaces", 1), ("Geometry.OCCFixDegenerated", 0), ("Geometry.OCCFixSmallEdges", 0),
             ("Geometry.OCCFixSmallFaces", 0), ("Geometry.OCCMakeSolids", 0), ("Geometry.Tolerance", tol_m / S)):
    gmsh.option.setNumber(k, v)
gmsh.model.add("sew")
gmsh.model.occ.importShapes(igs, highestDimOnly=True)
gmsh.model.occ.synchronize()

faces = sorted(t for _, t in gmsh.model.getEntities(2))
curves = [t for _, t in gmsh.model.getEntities(1)]
print(f"tolerance {tol_m * 1000:.2f} mm ({tol_m / S:.5f} model units): faces {len(faces)}, curves {len(curves)}, "
      f"points {len(gmsh.model.getEntities(0))}")

old, unmapped = {}, []
for t in faces:
    b = [v * S for v in gmsh.model.getBoundingBox(2, t)]
    bb = [b[0], b[3], b[1], b[4], b[2], b[5]]
    a = gmsh.model.occ.getMass(2, t) * S * S
    ranked = sorted((max(abs(x - y) for x, y in zip(bb, ref[o][0])), o) for o in ref)
    ok = ranked                         # the bounding box decides; area is only checked afterwards
    if ok[0][0] < 0.01 and (len(ok) == 1 or ok[1][0] > 3 * ok[0][0] + 1e-4):
        old[t] = ok[0][1]
    else:
        unmapped.append((t, bb, a, ranked[:3], ok[:3]))
claimed = {}
for t, o in old.items():
    claimed.setdefault(o, []).append(t)
dups = {o: ts for o, ts in claimed.items() if len(ts) > 1}
if unmapped or dups:
    print(f"ABORT: mapped {len(old)} of {len(faces)} faces; duplicates {dups}")
    for t, bb, a, ranked, ok in unmapped:
        print(f"  unmapped sewn face {t}: bbox [m] x {bb[0]:.3f}..{bb[1]:.3f} y {bb[2]:.3f}..{bb[3]:.3f} z {bb[4]:.3f}..{bb[5]:.3f}  area {a:.4f} m2")
        print("    nearest by bbox      : " + ", ".join(f"tag {o} (bbox diff {d * 1000:.2f} mm, table area {ref[o][1]:.3f})" for d, o in ranked))
        print("    passing the area test: " + (", ".join(f"tag {o} ({d * 1000:.2f} mm)" for d, o in ok) or "none"))
    missing = sorted(set(ref) - set(claimed))
    print(f"  original tags not claimed by any sewn face: {missing}")
    sys.exit(1)
moved = sum(1 for t, o in old.items() if t != o)
print(f"mapped 30/30 faces to original tags; {moved} faces carry a different tag after sewing")
bad_area = []
for t, o in old.items():
    a = gmsh.model.occ.getMass(2, t) * S * S
    if abs(a - ref[o][1]) > 0.01 * max(ref[o][1], 0.05):
        bad_area.append((o, round(a, 4), ref[o][1]))
print("faces whose sewn area differs from the table by more than 1 %:", bad_area if bad_area else "none")

adj, by_n = {}, {}
for c in curves:
    up, _ = gmsh.model.getAdjacencies(1, c)
    adj[c] = sorted(old[abs(u)] for u in up)
    by_n.setdefault(len(up), []).append(c)
print("curves by number of adjacent faces:", {k: len(v) for k, v in sorted(by_n.items())})
free_all = by_n.get(1, [])
print(f"free edges in total {len(free_all)}; on wing faces or face 30: {sum(1 for c in free_all if adj[c][0] in WING + [30])}")

print("\nface  neighbours (via shared curves)                 free edges")
for f in WING + [30]:
    nb = sorted({g for c in curves if f in adj[c] for g in adj[c] if g != f})
    fr = [c for c in free_all if adj[c][0] == f]
    desc = []
    for c in fr:
        bb = [v * S for v in gmsh.model.getBoundingBox(1, c)]
        desc.append(f"{gmsh.model.occ.getMass(1, c) * S:.3f} m @ x {bb[0]:.2f}..{bb[3]:.2f} y {bb[1]:.2f}..{bb[4]:.2f}")
    print(f"{f:4d}  {str(nb):44s} {len(fr)}  " + "; ".join(desc))
gmsh.finalize()
