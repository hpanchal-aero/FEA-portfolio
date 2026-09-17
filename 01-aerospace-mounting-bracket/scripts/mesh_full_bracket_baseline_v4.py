"""
Stage 5 -- Meshing, Level 1, v4.

Root cause of persistent order-2 curving bug (v3): setting
Mesh.ElementOrder=2 directly before generate(3) did not reliably
respect Mesh.SecondOrderLinear for all nodes (462 elements remained
curved for reasons not fully diagnosed). A subsequent attempt to
manually force-straighten midside nodes via direct setNode() calls
CATASTROPHICALLY corrupted the mesh (151,058 degenerate elements) --
that approach is abandoned entirely, not reused.

v4 uses Gmsh's documented order-elevation path instead:
  1. Generate mesh at ElementOrder=1 (linear tets only)
  2. Set Mesh.SecondOrderLinear=1
  3. Explicitly call gmsh.model.mesh.setOrder(2)
This is Gmsh's tested, separate code path for order elevation, distinct
from generating directly at order 2.

Same field/sizing setup as v3 (toe end-cap relief balls) -- only the
order-2 generation MECHANISM changes.
"""

import gmsh
import sys
import os
import numpy as np

INPUT_BREP = "mesh/full_bracket_study/full_bracket_geometry_final.brep"
MSH_OUT = "mesh/full_bracket_study/stage5_level1_v4.msh"
INP_RAW_OUT = "mesh/full_bracket_study/stage5_level1_v4_raw.inp"

HOLE_BORE_SURF = 15
R3_FILLET_SURF = 12
TOE_FILLET_SURF = 3

SIZE_FAR = 6.0
HOLE_SIZE_MIN = 0.5
HOLE_DIST_MIN = 2.0
HOLE_DIST_MAX = 5.0
R3_SIZE_MIN = 0.5
R3_DIST_MIN = 2.0
R3_DIST_MAX = 6.0
TOE_SIZE_MIN = 0.25
TOE_DIST_MIN = 1.5
TOE_DIST_MAX = 4.0

RELIEF_X = 38.0
RELIEF_Z = 1.0
RELIEF_RADIUS = 3.0
RELIEF_THICKNESS = 2.0
RELIEF_SIZE = 0.6

def main():
    if not os.path.exists(INPUT_BREP):
        print(f"ERROR: input file not found: {INPUT_BREP}")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("stage5_level1_v4")
    gmsh.model.occ.importShapes(INPUT_BREP)
    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=3)
    if len(volumes) != 1:
        print(f"ERROR: expected 1 volume, found {len(volumes)}")
        sys.exit(1)
    print("Geometry import OK: 1 volume confirmed.")

    for tag, expected_area, name in [
        (HOLE_BORE_SURF, 62.8319, "hole bore"),
        (R3_FILLET_SURF, 188.4956, "R3 fillet"),
        (TOE_FILLET_SURF, 31.4159, "R1 toe fillet"),
    ]:
        area = gmsh.model.occ.getMass(2, tag)
        ok = abs(area - expected_area) < 0.1
        print(f"  surf {tag} ({name}): area={area:.4f} [{'OK' if ok else 'MISMATCH -- STOP'}]")
        if not ok:
            gmsh.finalize()
            sys.exit(1)

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [HOLE_BORE_SURF])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "SurfacesList", [R3_FILLET_SURF])
    gmsh.model.mesh.field.setNumber(2, "Sampling", 100)
    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "SurfacesList", [TOE_FILLET_SURF])
    gmsh.model.mesh.field.setNumber(3, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 1)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", HOLE_SIZE_MIN)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", SIZE_FAR)
    gmsh.model.mesh.field.setNumber(4, "DistMin", HOLE_DIST_MIN)
    gmsh.model.mesh.field.setNumber(4, "DistMax", HOLE_DIST_MAX)
    gmsh.model.mesh.field.setNumber(4, "Sigmoid", 1)

    gmsh.model.mesh.field.add("Threshold", 5)
    gmsh.model.mesh.field.setNumber(5, "InField", 2)
    gmsh.model.mesh.field.setNumber(5, "SizeMin", R3_SIZE_MIN)
    gmsh.model.mesh.field.setNumber(5, "SizeMax", SIZE_FAR)
    gmsh.model.mesh.field.setNumber(5, "DistMin", R3_DIST_MIN)
    gmsh.model.mesh.field.setNumber(5, "DistMax", R3_DIST_MAX)
    gmsh.model.mesh.field.setNumber(5, "Sigmoid", 1)

    gmsh.model.mesh.field.add("Threshold", 6)
    gmsh.model.mesh.field.setNumber(6, "InField", 3)
    gmsh.model.mesh.field.setNumber(6, "SizeMin", TOE_SIZE_MIN)
    gmsh.model.mesh.field.setNumber(6, "SizeMax", SIZE_FAR)
    gmsh.model.mesh.field.setNumber(6, "DistMin", TOE_DIST_MIN)
    gmsh.model.mesh.field.setNumber(6, "DistMax", TOE_DIST_MAX)
    gmsh.model.mesh.field.setNumber(6, "Sigmoid", 1)

    gmsh.model.mesh.field.add("Ball", 8)
    gmsh.model.mesh.field.setNumber(8, "XCenter", RELIEF_X)
    gmsh.model.mesh.field.setNumber(8, "YCenter", 0.0)
    gmsh.model.mesh.field.setNumber(8, "ZCenter", RELIEF_Z)
    gmsh.model.mesh.field.setNumber(8, "Radius", RELIEF_RADIUS)
    gmsh.model.mesh.field.setNumber(8, "Thickness", RELIEF_THICKNESS)
    gmsh.model.mesh.field.setNumber(8, "VIn", RELIEF_SIZE)
    gmsh.model.mesh.field.setNumber(8, "VOut", 0.0)

    gmsh.model.mesh.field.add("Ball", 9)
    gmsh.model.mesh.field.setNumber(9, "XCenter", RELIEF_X)
    gmsh.model.mesh.field.setNumber(9, "YCenter", 40.0)
    gmsh.model.mesh.field.setNumber(9, "ZCenter", RELIEF_Z)
    gmsh.model.mesh.field.setNumber(9, "Radius", RELIEF_RADIUS)
    gmsh.model.mesh.field.setNumber(9, "Thickness", RELIEF_THICKNESS)
    gmsh.model.mesh.field.setNumber(9, "VIn", RELIEF_SIZE)
    gmsh.model.mesh.field.setNumber(9, "VOut", 0.0)

    gmsh.model.mesh.field.add("Max", 10)
    gmsh.model.mesh.field.setNumbers(10, "FieldsList", [6, 8, 9])

    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [4, 5, 10])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # --- Step 1: generate at ORDER 1 (linear tets only) ---
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    print("\nGenerating LINEAR (order 1) 3D mesh...")
    gmsh.model.mesh.generate(3)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_lin_elements = sum(len(t) for t in elem_tags)
    print(f"Linear mesh: {len(node_tags)} nodes, {n_lin_elements} elements")

    # --- Step 2: set SecondOrderLinear, THEN explicitly elevate order ---
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    print("\nElevating to order 2 via explicit setOrder(2) call "
          "(documented order-elevation path, not direct order-2 generation)...")
    gmsh.model.mesh.setOrder(2)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    n_nodes = len(node_tags)
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(t) for t in elem_tags)
    print(f"\nOrder-2 mesh: {n_nodes} nodes, {n_elements} 3D elements")

    est_equations = 3 * n_nodes
    ceiling = 1_180_000
    pct = 100 * est_equations / ceiling
    print(f"Rough equation estimate: {est_equations} ({pct:.1f}% of ceiling)")

    tet_idx = [i for i, et in enumerate(elem_types) if et == 11][0]
    quals = np.array(gmsh.model.mesh.getElementQualities(elem_tags[tet_idx], "minSICN"))
    print(f"\nElement quality (minSICN): min={quals.min():.4f}, "
          f"mean={quals.mean():.4f}, max={quals.max():.4f}")
    n_bad = int((quals < 0.1).sum())
    print(f"Elements with minSICN < 0.1: {n_bad} ({100*n_bad/len(quals):.3f}%)")

    # --- Direct verification: check straightness on the SAME sample tags
    # that were problematic in v3, to directly confirm this mechanism worked ---
    print("\nVerifying straightness on previously-problematic element tags...")
    etags_arr = np.array(elem_tags[tet_idx])
    enodes_arr = np.array(gmsh.model.mesh.getElements(3)[2][tet_idx]).reshape(len(etags_arr), 10)
    edge_pairs = [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
    for target_tag in [170300, 170400, 170500, 170600, 170700]:
        match = np.where(etags_arr == target_tag)[0]
        if len(match) == 0:
            print(f"  elem {target_tag}: not found (tag numbering likely shifted -- expected)")
            continue
        idx = match[0]
        nn = enodes_arr[idx]
        coords = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in nn])
        max_dev = 0.0
        for i, (a, b) in enumerate(edge_pairs):
            straight_mid = (coords[a] + coords[b]) / 2.0
            dev = np.linalg.norm(coords[4+i] - straight_mid)
            max_dev = max(max_dev, dev)
        print(f"  elem {target_tag}: max midside deviation = {max_dev:.6f} mm")

    os.makedirs(os.path.dirname(MSH_OUT), exist_ok=True)
    gmsh.write(MSH_OUT)
    print(f"\nSaved: {MSH_OUT}")
    gmsh.write(INP_RAW_OUT)
    print(f"Saved: {INP_RAW_OUT}")

    gmsh.finalize()

if __name__ == "__main__":
    main()
