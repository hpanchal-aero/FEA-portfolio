"""
Stage 5 -- Meshing, Level 1, LINEAR (C3D4) TEMPORARY BASELINE.

Uses the same verified geometry, field setup, and sizing as v3/v4,
but stops at ElementOrder=1 (no order elevation) -- deliberately
avoiding the still-unexplained order-2 midside-node defect (462
elements with deterministic, curving-setting-independent deviations,
suspected periodic-entity artifact from the y=0/y=40 extrusion,
not yet root-caused after 9 diagnostic rounds).

STATUS: TEMPORARY. C3D4 (linear tet) elements are known to exhibit
shear locking / excess bending stiffness -- NOT appropriate for
Stage 5's final reported results on a bracket under bending load.
This mesh exists solely to unblock baseline-solve verification
(equilibrium check, sanity-level deflection/stress) while the
order-2 periodicity issue is investigated separately. Must be
superseded by a working C3D10 mesh before Stage 5 conclusions
are finalized -- flag explicitly in README, do not let this
silently become the reported result.
"""

import gmsh
import sys
import os
import numpy as np

INPUT_BREP = "mesh/full_bracket_study/full_bracket_geometry_final.brep"
MSH_OUT = "mesh/full_bracket_study/stage5_level1_LINEAR_TEMP.msh"
INP_RAW_OUT = "mesh/full_bracket_study/stage5_level1_LINEAR_TEMP_raw.inp"

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
    gmsh.model.add("stage5_level1_linear_temp")
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
    gmsh.option.setNumber("Mesh.ElementOrder", 1)  # LINEAR -- deliberate, temporary

    print("\nGenerating LINEAR (C3D4) 3D mesh -- TEMPORARY baseline only...")
    gmsh.model.mesh.generate(3)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    n_nodes = len(node_tags)
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(t) for t in elem_tags)
    print(f"\nLinear mesh: {n_nodes} nodes, {n_elements} elements")

    est_equations = 3 * n_nodes
    ceiling = 1_180_000
    print(f"Rough equation estimate: {est_equations} ({100*est_equations/ceiling:.1f}% of ceiling)")

    tet_idx = [i for i, et in enumerate(elem_types) if et == 4][0]  # type 4 = 4-node tet
    quals = np.array(gmsh.model.mesh.getElementQualities(elem_tags[tet_idx], "minSICN"))
    print(f"Element quality (minSICN): min={quals.min():.4f}, mean={quals.mean():.4f}, "
          f"max={quals.max():.4f}")
    n_bad = int((quals < 0.1).sum())
    print(f"Elements with minSICN < 0.1: {n_bad} ({100*n_bad/len(quals):.3f}%)")

    os.makedirs(os.path.dirname(MSH_OUT), exist_ok=True)
    gmsh.write(MSH_OUT)
    print(f"\nSaved: {MSH_OUT}")
    gmsh.write(INP_RAW_OUT)
    print(f"Saved: {INP_RAW_OUT}")
    print("\n*** REMINDER: this is a TEMPORARY C3D4 baseline. Not the final")
    print("*** Stage 5 mesh formulation. Order-2 periodicity issue remains open.")

    gmsh.finalize()

if __name__ == "__main__":
    main()
