"""
Direct test: for a sample of the persistently-failing elements
(170259-170720), check whether their midside nodes sit exactly at
the straight-line edge midpoint of their parent corner nodes.

If YES (deviation ~0): SecondOrderLinear=1 took effect correctly,
and per the geometric identity (straight midside nodes -> quadratic
map == linear map), these elements are mathematically guaranteed
non-degenerate -- meaning the Jacobian failure must have some OTHER
cause we haven't found (e.g. a parsing/export bug, not a geometric one).

If NO (nonzero deviation): SecondOrderLinear did NOT actually apply,
despite being set -- meaning curving is still happening somehow, and
we need to find out why the option didn't take effect.
"""

import gmsh
import numpy as np

MSH_FILE = "mesh/full_bracket_study/stage5_level1_v3.msh"
SAMPLE_TAGS = [170300, 170400, 170500, 170600, 170700]

gmsh.initialize()
gmsh.open(MSH_FILE)

etypes, etags, enodes = gmsh.model.mesh.getElements(3)
tet_idx = [i for i, et in enumerate(etypes) if et == 11][0]
etags_arr = np.array(etags[tet_idx])
enodes_arr = np.array(enodes[tet_idx]).reshape(len(etags_arr), 10)

# C3D10/gmsh mid-edge node order (positions 4-9) corresponds to edges:
# (0,1),(1,2),(2,0),(0,3),(1,3),(2,3)  [gmsh 10-node tet convention]
edge_pairs = [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]

for target_tag in SAMPLE_TAGS:
    match = np.where(etags_arr == target_tag)[0]
    if len(match) == 0:
        print(f"Element {target_tag}: NOT FOUND in this mesh (tag may not exist)")
        continue
    idx = match[0]
    nn = enodes_arr[idx]
    coords = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in nn])

    print(f"\nElement {target_tag}:")
    max_dev = 0.0
    for i, (a, b) in enumerate(edge_pairs):
        corner_mid = (coords[a] + coords[b]) / 2.0
        actual_mid = coords[4 + i]
        dev = np.linalg.norm(actual_mid - corner_mid)
        max_dev = max(max_dev, dev)
        print(f"  edge({a},{b}) -> midside node {int(nn[4+i])}: "
              f"deviation from straight midpoint = {dev:.6f} mm")
    print(f"  MAX deviation this element: {max_dev:.6f} mm "
          f"[{'STRAIGHT (as expected)' if max_dev < 1e-6 else 'CURVED -- SecondOrderLinear NOT applied!'}]")

gmsh.finalize()
