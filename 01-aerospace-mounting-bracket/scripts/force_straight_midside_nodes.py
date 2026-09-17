"""
Post-process fix: force EVERY midside node in the mesh to sit exactly
at the straight-line midpoint of its two parent corner nodes, by
direct coordinate overwrite -- rather than relying on Gmsh's
SecondOrderLinear option, whose behavior in this mesh did not match
its documented effect (some midside nodes remained off-center for
reasons not yet understood, including nodes classified as volume-
interior with no curved-geometry justification).

This guarantees, by direct construction, that every element is
straight-sided, which mathematically guarantees a Jacobian identical
to the (already-verified-positive) linear tetrahedron's Jacobian --
eliminating the CalculiX nonpositive-Jacobian failures regardless of
their underlying cause.
"""

import gmsh
import numpy as np

MSH_IN = "mesh/full_bracket_study/stage5_level1_v3.msh"
MSH_OUT = "mesh/full_bracket_study/stage5_level1_v3_straightened.msh"

edge_pairs = [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]

gmsh.initialize()
gmsh.open(MSH_IN)

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_tags = np.array(node_tags, dtype=int)
coords = np.array(node_coords).reshape(-1, 3)
coord_map = {int(t): c for t, c in zip(node_tags, coords)}

etypes, etags, enodes = gmsh.model.mesh.getElements(3)
tet_idx = [i for i, et in enumerate(etypes) if et == 11][0]
etags_arr = np.array(etags[tet_idx])
enodes_arr = np.array(enodes[tet_idx]).reshape(len(etags_arr), 10)

print(f"Processing {len(etags_arr)} elements...")

corrections = {}
max_dev_found = 0.0
n_corrected = 0

for nn in enodes_arr:
    for i, (a, b) in enumerate(edge_pairs):
        corner_a = coord_map[int(nn[a])]
        corner_b = coord_map[int(nn[b])]
        mid_tag = int(nn[4 + i])
        straight_mid = (corner_a + corner_b) / 2.0
        current = coord_map[mid_tag]
        dev = np.linalg.norm(current - straight_mid)
        if dev > max_dev_found:
            max_dev_found = dev
        if dev > 1e-9:
            corrections[mid_tag] = straight_mid
            n_corrected += 1

print(f"Max deviation found before correction: {max_dev_found:.6f} mm")
print(f"Midside nodes requiring correction: {n_corrected} / "
      f"{len(node_tags) - len(etags_arr)*4//10} (approx total midside count)")

# Apply corrections via setNode (overwrite coordinates directly)
for tag, new_coord in corrections.items():
    gmsh.model.mesh.setNode(tag, new_coord.tolist(), [])

print(f"Applied {len(corrections)} corrections.")

gmsh.write(MSH_OUT)
print(f"Saved corrected mesh: {MSH_OUT}")

gmsh.finalize()
