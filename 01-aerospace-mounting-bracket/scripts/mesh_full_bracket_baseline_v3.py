"""
Stage 5 -- Meshing, Level 1, v3: targeted local relief at the toe-fillet
curve's two end-cap vertices (y=0 and y=40), where 14 of 18 poor-quality
elements clustered (including the two worst, minSICN 0.0064/0.0074).

Mechanism (from diagnosis): the toe-fillet's curved profile meets the
flat y=0/y=40 end-cap faces at a tight trihedral corner. The aggressive
TOE_SIZE_MIN=0.25mm target there, combined with the surrounding size
gradient, produces near-degenerate tets that Delaunay can't avoid.

Fix: two Ball fields, centered on the toe-fillet region at y=0 and
y=40, that LOCALLY RELAX the size floor (via Max-combination with the
existing toe threshold field) only within a small radius of those two
vertices. Everywhere else, sizing is unchanged from v2.

The 2 remaining bad elements near the R3 fillet (moderate severity,
0.049/0.060, mid-span not at a corner) are NOT addressed by this fix --
flagged to re-check after this run, not assumed resolved.
"""

import gmsh
import sys
import os
import numpy as np

INPUT_BREP = "mesh/full_bracket_study/full_bracket_geometry_final.brep"
MSH_OUT = "mesh/full_bracket_study/stage5_level1_v3.msh"
INP_RAW_OUT = "mesh/full_bracket_study/stage5_level1_v3_raw.inp"

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

# --- New: local relief at toe-fillet end-cap vertices ---
# Toe fillet surface bbox (confirmed earlier): x=[37.5858,38.2929],
# z=[2.0000,2.2929]. Relief ball centered at the toe-fillet region,
# at each end-cap plane.
RELIEF_X = 38.0
RELIEF_Z = 1.0   # spans a bit below/above the fillet to cover the
                  # broader bad-element footprint (x=35-39, z=0-2)
RELIEF_RADIUS = 3.0
RELIEF_THICKNESS = 2.0
RELIEF_SIZE = 0.6   # relaxed floor at these two points only

def main():
    if not os.path.exists(INPUT_BREP):
        print(f"ERROR: input file not found: {INPUT_BREP}")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("stage5_level1_v3")
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

    # --- Distance + Threshold fields, same as v2 ---
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

    # --- New: relief Ball fields at y=0 and y=40 toe end-caps ---
    gmsh.model.mesh.field.add("Ball", 8)
    gmsh.model.mesh.field.setNumber(8, "XCenter", RELIEF_X)
    gmsh.model.mesh.field.setNumber(8, "YCenter", 0.0)
    gmsh.model.mesh.field.setNumber(8, "ZCenter", RELIEF_Z)
    gmsh.model.mesh.field.setNumber(8, "Radius", RELIEF_RADIUS)
    gmsh.model.mesh.field.setNumber(8, "Thickness", RELIEF_THICKNESS)
    gmsh.model.mesh.field.setNumber(8, "VIn", RELIEF_SIZE)
    gmsh.model.mesh.field.setNumber(8, "VOut", 0.0)  # inert (small, not large) outside radius+thickness

    gmsh.model.mesh.field.add("Ball", 9)
    gmsh.model.mesh.field.setNumber(9, "XCenter", RELIEF_X)
    gmsh.model.mesh.field.setNumber(9, "YCenter", 40.0)
    gmsh.model.mesh.field.setNumber(9, "ZCenter", RELIEF_Z)
    gmsh.model.mesh.field.setNumber(9, "Radius", RELIEF_RADIUS)
    gmsh.model.mesh.field.setNumber(9, "Thickness", RELIEF_THICKNESS)
    gmsh.model.mesh.field.setNumber(9, "VIn", RELIEF_SIZE)
    gmsh.model.mesh.field.setNumber(9, "VOut", 0.0)

    # Max-combine: picks the LARGER (coarser) of field6 and the two relief
    # balls -- i.e. relief only ever coarsens locally, never overrides
    # field6's finer sizing anywhere the balls aren't active.
    gmsh.model.mesh.field.add("Max", 10)
    gmsh.model.mesh.field.setNumbers(10, "FieldsList", [6, 8, 9])

    # --- Combine via Min, using field 10 in place of field 6 ---
    gmsh.model.mesh.field.add("Min", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [4, 5, 10])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)  # straight-sided midside nodes -- eliminates curving-induced invalid elements

    print("\nGenerating 3D mesh (v3, with toe end-cap relief)...")
    gmsh.model.mesh.generate(3)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    n_nodes = len(node_tags)
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(t) for t in elem_tags)
    print(f"\nMesh generated: {n_nodes} nodes, {n_elements} 3D elements")

    est_equations = 3 * n_nodes
    ceiling = 1_180_000
    pct = 100 * est_equations / ceiling
    print(f"Rough equation estimate: {est_equations} ({pct:.1f}% of ceiling)")
    if est_equations > ceiling:
        print("*** WARNING: exceeds ceiling. Do NOT solve. ***")
    else:
        print("OK: under ceiling.")

    quals = np.array(gmsh.model.mesh.getElementQualities(elem_tags[0], "minSICN"))
    print(f"\nElement quality (minSICN): min={quals.min():.4f}, "
          f"mean={quals.mean():.4f}, max={quals.max():.4f}")
    n_bad = int((quals < 0.1).sum())
    print(f"Elements with minSICN < 0.1: {n_bad} ({100*n_bad/len(quals):.3f}%)")

    if n_bad > 0:
        print("\nStill-remaining bad elements -- locations:")
        etags = elem_tags[0]
        n_nodes_per_elem = len(node_tags_flat := gmsh.model.mesh.getElements(3)[2][0]) // len(etags)
        ntags_reshaped = np.array(node_tags_flat).reshape(-1, n_nodes_per_elem)
        bad_idx = np.where(quals < 0.1)[0]
        for idx in bad_idx:
            tag = etags[idx]
            q = quals[idx]
            node_ids = ntags_reshaped[idx][:4]
            coords = [gmsh.model.mesh.getNode(int(nid))[0] for nid in node_ids]
            centroid = np.array(coords).mean(axis=0)
            print(f"  elem {tag}: minSICN={q:.6f}  centroid=({centroid[0]:.3f},"
                  f"{centroid[1]:.3f},{centroid[2]:.3f})")

    os.makedirs(os.path.dirname(MSH_OUT), exist_ok=True)
    gmsh.write(MSH_OUT)
    print(f"\nSaved: {MSH_OUT}")
    gmsh.write(INP_RAW_OUT)
    print(f"Saved: {INP_RAW_OUT}")

    gmsh.finalize()

if __name__ == "__main__":
    main()
