"""
Stage 5 -- diagnose the 18 poor-quality (minSICN < 0.1) elements found
in the Level 1 v2 mesh, including the near-degenerate minSICN=0.0064
element. Locate their centroids to determine whether they cluster at
a specific feature or at an overlap boundary between refinement fields.
"""

import gmsh
import numpy as np

MSH_FILE = "mesh/full_bracket_study/stage5_level1.msh"

gmsh.initialize()
gmsh.open(MSH_FILE)

elem_types, elem_tags, node_tags_flat = gmsh.model.mesh.getElements(3)

# Find the tet element type block (type 11 = 10-node tet, since order=2)
for et, etags, ntags in zip(elem_types, elem_tags, node_tags_flat):
    quals = np.array(gmsh.model.mesh.getElementQualities(etags, "minSICN"))
    bad_mask = quals < 0.1
    n_bad = bad_mask.sum()
    print(f"Element type {et}: {len(etags)} elements, {n_bad} with minSICN < 0.1")

    if n_bad == 0:
        continue

    bad_indices = np.where(bad_mask)[0]
    n_nodes_per_elem = len(ntags) // len(etags)
    ntags_reshaped = np.array(ntags).reshape(-1, n_nodes_per_elem)

    print(f"\n{'ElemTag':>10} {'minSICN':>10} {'Centroid (x,y,z)':>30}")
    for idx in bad_indices:
        tag = etags[idx]
        q = quals[idx]
        node_ids = ntags_reshaped[idx][:4]  # corner nodes only, sufficient for centroid
        coords = []
        for nid in node_ids:
            coord, _, _, _ = gmsh.model.mesh.getNode(int(nid))
            coords.append(coord)
        coords = np.array(coords)
        centroid = coords.mean(axis=0)
        print(f"{tag:>10} {q:>10.6f}   ({centroid[0]:7.3f}, {centroid[1]:7.3f}, {centroid[2]:7.3f})")

    print("\nProximity check -- distance from each bad element's centroid to each")
    print("feature's nominal location (helps identify which field boundary is implicated):")
    R3_CENTER_APPROX = np.array([56.5, 20.0, 3.5])   # rough R3 fillet region center
    TOE_APPROX = np.array([38.0, 20.0, 2.0])          # toe fillet location
    HOLE_APPROX = np.array([20.0, 20.0, 0.0])         # hole center

    for idx in bad_indices:
        tag = etags[idx]
        node_ids = ntags_reshaped[idx][:4]
        coords = []
        for nid in node_ids:
            coord, _, _, _ = gmsh.model.mesh.getNode(int(nid))
            coords.append(coord)
        centroid = np.array(coords).mean(axis=0)
        d_r3 = np.linalg.norm(centroid - R3_CENTER_APPROX)
        d_toe = np.linalg.norm(centroid - TOE_APPROX)
        d_hole = np.linalg.norm(centroid - HOLE_APPROX)
        nearest = min([("R3 fillet", d_r3), ("toe fillet", d_toe), ("hole", d_hole)],
                       key=lambda p: p[1])
        print(f"  elem {tag}: nearest feature = {nearest[0]} (dist ~{nearest[1]:.2f}mm), "
              f"[R3={d_r3:.2f}, toe={d_toe:.2f}, hole={d_hole:.2f}]")

gmsh.finalize()
