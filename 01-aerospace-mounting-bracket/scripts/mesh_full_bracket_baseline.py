"""
Stage 5 -- Meshing, Level 1 (baseline/coarsest), REVISED SIZING.

v1 overshot the 1.18M-equation ceiling by 278.9% (1,096,817 nodes).
Root cause: DistMax falloff zones for the three features (8-12mm) were
too large relative to the feature separations, causing overlapping
fine-mesh regions to cover most of the bracket volume, not just a thin
skin around each feature. Quality was excellent (minSICN min=0.20,
mean=0.82) -- this is a sizing-magnitude fix only, not a field-pattern
change.

Changes vs v1:
  - SIZE_FAR: 3.0 -> 6.0mm (coarser background, most of the volume)
  - All DistMax values roughly halved (tighter falloff zones)
  - Near-feature SizeMin values held close to original (toe fillet
    resolution is not sacrificed, since it has no convergence precedent)
"""

import gmsh
import sys
import os

INPUT_BREP = "mesh/full_bracket_study/full_bracket_geometry_final.brep"
MSH_OUT = "mesh/full_bracket_study/stage5_level1.msh"
INP_RAW_OUT = "mesh/full_bracket_study/stage5_level1_raw.inp"

HOLE_BORE_SURF = 15
R3_FILLET_SURF = 12
TOE_FILLET_SURF = 3

# --- Level 1 (baseline), REVISED sizing ---
SIZE_FAR = 6.0  # was 3.0

HOLE_SIZE_MIN = 0.5   # unchanged
HOLE_DIST_MIN = 2.0   # was 3.0
HOLE_DIST_MAX = 5.0   # was 10.0

R3_SIZE_MIN = 0.5     # was 0.4 (slightly relaxed)
R3_DIST_MIN = 2.0     # was 3.0
R3_DIST_MAX = 6.0     # was 12.0

TOE_SIZE_MIN = 0.25   # was 0.2 (kept tight -- no convergence precedent)
TOE_DIST_MIN = 1.5    # was 2.0
TOE_DIST_MAX = 4.0    # was 8.0

def main():
    if not os.path.exists(INPUT_BREP):
        print(f"ERROR: input file not found: {INPUT_BREP}")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("stage5_level1_v2")
    gmsh.model.occ.importShapes(INPUT_BREP)
    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=3)
    if len(volumes) != 1:
        print(f"ERROR: expected 1 volume, found {len(volumes)}")
        sys.exit(1)
    print(f"Geometry import OK: 1 volume confirmed.")

    for tag, expected_area, name in [
        (HOLE_BORE_SURF, 62.8319, "hole bore"),
        (R3_FILLET_SURF, 188.4956, "R3 fillet"),
        (TOE_FILLET_SURF, 31.4159, "R1 toe fillet"),
    ]:
        area = gmsh.model.occ.getMass(2, tag)
        ok = abs(area - expected_area) < 0.1
        print(f"  surf {tag} ({name}): area={area:.4f}, expected={expected_area:.4f} "
              f"[{'OK' if ok else 'MISMATCH -- STOP'}]")
        if not ok:
            print("ERROR: surface tag/area mismatch. Aborting.")
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

    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [4, 5, 6])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    print("\nGenerating 3D mesh (revised sizing)...")
    gmsh.model.mesh.generate(3)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    n_nodes = len(node_tags)
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(t) for t in elem_tags)
    print(f"\nMesh generated: {n_nodes} nodes, {n_elements} 3D elements")

    est_equations = 3 * n_nodes
    ceiling = 1_180_000
    pct_of_ceiling = 100 * est_equations / ceiling
    print(f"Rough equation estimate (3 x nodes): {est_equations}")
    print(f"Percent of ceiling: {pct_of_ceiling:.1f}%")
    if est_equations > ceiling:
        print("*** WARNING: still exceeds ceiling. Do NOT solve. Coarsen further. ***")
    elif pct_of_ceiling > 70:
        print("CAUTION: within 70% of ceiling -- limited headroom for refinement levels.")
    else:
        print("OK: comfortably under ceiling for a Level 1 baseline mesh.")

    import numpy as np
    quals = np.array(gmsh.model.mesh.getElementQualities(elem_tags[0], "minSICN"))
    print(f"\nElement quality (minSICN): min={quals.min():.4f}, "
          f"mean={quals.mean():.4f}, max={quals.max():.4f}")
    n_poor = int((quals < 0.1).sum())
    print(f"Elements with minSICN < 0.1: {n_poor} ({100*n_poor/len(quals):.2f}%)")

    os.makedirs(os.path.dirname(MSH_OUT), exist_ok=True)
    gmsh.write(MSH_OUT)
    print(f"\nSaved: {MSH_OUT}")
    gmsh.write(INP_RAW_OUT)
    print(f"Saved: {INP_RAW_OUT}")
    print("\nRAW .inp -- needs cleanup pass before solving, same as before.")

    gmsh.finalize()

if __name__ == "__main__":
    main()
