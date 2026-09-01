"""
Project 01 - Aerospace Mounting Bracket
Diagnostic: inspect physical groups in the verification plate mesh,
confirming each physical face maps to the correct combined bounding box.
"""

import gmsh
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
mesh_path = os.path.join(script_dir, "..", "mesh", "verification_plate.msh")

gmsh.initialize()
gmsh.open(mesh_path)

print("=== Physical Groups ===")
groups = gmsh.model.getPhysicalGroups()
for dim, tag in groups:
    name = gmsh.model.getPhysicalName(dim, tag)
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    xmin = ymin = zmin = 1e9
    xmax = ymax = zmax = -1e9
    for e in entities:
        bb = gmsh.model.getBoundingBox(dim, e)
        xmin, ymin, zmin = min(xmin, bb[0]), min(ymin, bb[1]), min(zmin, bb[2])
        xmax, ymax, zmax = max(xmax, bb[3]), max(ymax, bb[4]), max(zmax, bb[5])
    print(f"'{name}' (dim={dim}, tag={tag}): entities={list(entities)} "
          f"combined bbox x=[{xmin:.3f},{xmax:.3f}] "
          f"y=[{ymin:.3f},{ymax:.3f}] z=[{zmin:.3f},{zmax:.3f}]")

gmsh.finalize()
