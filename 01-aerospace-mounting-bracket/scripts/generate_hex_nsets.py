"""
Project 01 - Aerospace Mounting Bracket
Generate *NSET definitions for the structured hex mesh, in the SAME
Gmsh session as mesh generation to guarantee node numbering
consistency with mesh_hex.inp (same rationale as
generate_calculix_mesh.py for the tet mesh).
"""

import gmsh
import os

GEO_PATH = os.path.abspath("../mesh/hex_test/case_hex.geo")
NSET_OUT = os.path.abspath("../mesh/hex_test/nsets_hex.inp")

gmsh.initialize()
gmsh.open(GEO_PATH)

named_surfaces = ["root", "tip", "bottom", "top", "side_y0", "side_y40"]
groups = gmsh.model.getPhysicalGroups(dim=2)
name_to_tag = {gmsh.model.getPhysicalName(d, t): t for d, t in groups}

with open(NSET_OUT, "w") as f:
    for name in named_surfaces:
        tag = name_to_tag[name]
        node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
        node_tags = sorted(int(n) for n in node_tags)
        print(f"NSET '{name}': {len(node_tags)} nodes")
        f.write(f"*NSET,NSET={name}\n")
        for i in range(0, len(node_tags), 8):
            chunk = node_tags[i:i+8]
            f.write(",".join(str(n) for n in chunk) + ",\n")

print(f"\nWrote: {NSET_OUT}")
gmsh.finalize()
