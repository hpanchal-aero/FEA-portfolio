#!/usr/bin/env python3
"""List the faces of an IGES surface model: Gmsh type, bounding box [m], area [m2], boundary-curve count.

usage: list_faces.py <file.igs> [scale]      scale = metres per model unit (default 1/25.4, the uCRM-9 'wbt' file; wingbox file: 0.001)
Read-only.
"""
import sys

import gmsh

path = sys.argv[1]
scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0 / 25.4

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("faces")
gmsh.model.occ.importShapes(path, highestDimOnly=True)
gmsh.model.occ.synchronize()

faces = sorted(t for _, t in gmsh.model.getEntities(2))
print(f"{path}: {len(faces)} faces, scale {scale:.6f} m per model unit")
print(" tag type                     xmin     xmax     ymin     ymax     zmin     zmax   area[m2]  curves")
tot = 0.0
for t in faces:
    b = [v * scale for v in gmsh.model.getBoundingBox(2, t)]
    a = gmsh.model.occ.getMass(2, t) * scale * scale
    tot += a
    nc = len(gmsh.model.getBoundary([(2, t)], oriented=False, combined=False))
    print(f"{t:4d} {gmsh.model.getType(2, t):18s} " + " ".join(f"{v:8.3f}" for v in (b[0], b[3], b[1], b[4], b[2], b[5]))
          + f" {a:10.3f} {nc:6d}")
print(f"total area {tot:.2f} m2")
gmsh.finalize()
