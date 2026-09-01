"""
Project 01 - Aerospace Mounting Bracket
Stage: Mesh -> CalculiX .inp export with correct NSETs.

Purpose: Gmsh's native Abaqus/.inp writer exports Physical Surfaces as
*ELSET (facet elements), not *NSET (nodes). CalculiX boundary conditions
and concentrated loads require *NSET. This script runs geometry + mesh
generation and .inp export in a single Gmsh session so that node
numbering in the exported mesh and the derived node sets are guaranteed
consistent (Delaunay meshing is non-deterministic across separate runs).

Outputs:
  mesh/verification_plate.inp        - mesh (nodes + elements) only
  mesh/verification_plate_nsets.inp  - *NSET definitions for each
                                        named physical surface, to be
                                        *INCLUDE-d in the analysis .inp
"""

import gmsh
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
geo_path = os.path.join(script_dir, "..", "mesh", "01_verification_plate.geo")
mesh_out = os.path.join(script_dir, "..", "mesh", "verification_plate.inp")
nset_out = os.path.join(script_dir, "..", "mesh", "verification_plate_nsets.inp")

gmsh.initialize()
gmsh.open(geo_path)   # executes the .geo script, including "Mesh 3;"

# --- Export mesh (nodes + elements) in Abaqus/CalculiX format ---
gmsh.write(mesh_out)

# --- Derive node sets from the SAME in-memory mesh, guaranteed consistent ---
named_surfaces = ["root", "tip", "bottom", "top", "side_y0", "side_y40"]

groups = gmsh.model.getPhysicalGroups(dim=2)
name_to_tag = {gmsh.model.getPhysicalName(d, t): t for d, t in groups}

with open(nset_out, "w") as f:
    for name in named_surfaces:
        if name not in name_to_tag:
            raise RuntimeError(f"Physical surface '{name}' not found in model")
        tag = name_to_tag[name]
        node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
        node_tags = sorted(int(n) for n in node_tags)

        print(f"NSET '{name}': {len(node_tags)} nodes")

        f.write(f"*NSET,NSET={name}\n")
        for i in range(0, len(node_tags), 8):
            chunk = node_tags[i:i+8]
            f.write(",".join(str(n) for n in chunk) + ",\n")

print(f"\nWrote: {mesh_out}")
print(f"Wrote: {nset_out}")

gmsh.finalize()
