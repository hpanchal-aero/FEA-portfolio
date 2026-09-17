"""
Test for coordinate-rounding collisions: check whether the 6-decimal
rounding used when writing node coordinates to the .inp causes any
two DISTINCT node tags to receive IDENTICAL rounded coordinates,
which would create a mathematically degenerate element invisible to
Gmsh's own full-precision quality check but fatal to CalculiX.
"""

import gmsh
import numpy as np
from collections import defaultdict

MSH_FILE = "mesh/full_bracket_study/stage5_level1_v3.msh"

gmsh.initialize()
gmsh.open(MSH_FILE)

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_tags = np.array(node_tags, dtype=int)
coords = np.array(node_coords).reshape(-1, 3)

print(f"Total nodes: {len(node_tags)}")

# Round to 6 decimals, exactly as build_full_bracket_inp.py does
rounded = np.round(coords, 6)

seen = defaultdict(list)
for tag, c in zip(node_tags, rounded):
    key = (c[0], c[1], c[2])
    seen[key].append(int(tag))

collisions = {k: v for k, v in seen.items() if len(v) > 1}
print(f"\nCoordinate collisions after 6-decimal rounding: {len(collisions)}")

if collisions:
    print("\nSample collisions (rounded coord -> colliding node tags):")
    for i, (k, v) in enumerate(collisions.items()):
        if i >= 20:
            print(f"  ... and {len(collisions)-20} more")
            break
        print(f"  {k}: nodes {v}")

    # Check true (unrounded) separation for a few collision groups
    print("\nTrue (full-precision) separation for first 5 collision groups:")
    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}
    for i, (k, v) in enumerate(collisions.items()):
        if i >= 5:
            break
        c0 = coords[tag_to_idx[v[0]]]
        c1 = coords[tag_to_idx[v[1]]]
        dist = np.linalg.norm(c0 - c1)
        print(f"  nodes {v[0]} vs {v[1]}: true distance = {dist:.10f} mm "
              f"(full coords: {c0} vs {c1})")
else:
    print("No collisions found -- this hypothesis is not the cause.")

gmsh.finalize()
