"""
Diagnose the 462 persistent nonpositive-Jacobian elements by checking
DIRECT NODE MEMBERSHIP against the actual curved surface meshes (toe
fillet=surf 3, R3 fillet=surf 12, hole bore=surf 15), rather than
distance to an approximate center point -- the previous check's blind
spot. A coarse element with one node on a tight-radius curved surface
can be far (by centroid) from any "center" while still triggering the
curvilinear midside-node projection failure.
"""

import gmsh
import numpy as np

MSH_FILE = "mesh/full_bracket_study/stage5_level1_v3.msh"
BAD_START, BAD_END = 170259, 170720

gmsh.initialize()
gmsh.open(MSH_FILE)

# --- Collect node IDs belonging to each curved surface's 2D mesh ---
surface_node_sets = {}
for tag, name in [(3, "toe fillet"), (12, "R3 fillet"), (15, "hole bore")]:
    etypes, etags, enodes = gmsh.model.mesh.getElements(2, tag)
    node_set = set()
    for et, ntags in zip(etypes, enodes):
        node_set.update(int(n) for n in ntags)
    surface_node_sets[name] = node_set
    print(f"{name} (surf {tag}): {len(node_set)} nodes in its surface mesh")

# --- Get the bad elements' node connectivity ---
etypes3, etags3, enodes3 = gmsh.model.mesh.getElements(3)
tet_idx = [i for i, et in enumerate(etypes3) if et == 11][0]
etags = np.array(etags3[tet_idx])
enodes = np.array(enodes3[tet_idx]).reshape(len(etags), 10)

bad_mask = (etags >= BAD_START) & (etags <= BAD_END)
bad_etags = etags[bad_mask]
bad_enodes = enodes[bad_mask]
print(f"\nBad elements found in mesh: {len(bad_etags)}")

touch_counts = {name: 0 for name in surface_node_sets}
touch_none = 0
for tag, nn in zip(bad_etags, bad_enodes):
    elem_nodes = set(int(n) for n in nn)
    touched_any = False
    for name, node_set in surface_node_sets.items():
        if elem_nodes & node_set:
            touch_counts[name] += 1
            touched_any = True
    if not touched_any:
        touch_none += 1

print("\nBad elements touching each curved surface (by shared node membership):")
for name, count in touch_counts.items():
    print(f"  {name}: {count} / {len(bad_etags)} ({100*count/len(bad_etags):.1f}%)")
print(f"  touching NONE of the three curved surfaces: {touch_none} / {len(bad_etags)} "
      f"({100*touch_none/len(bad_etags):.1f}%)")

gmsh.finalize()
