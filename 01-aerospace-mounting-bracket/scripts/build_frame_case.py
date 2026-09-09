"""
Stage 4a (right-angle frame, sharp corner, no gusset) - geometry
construction and baseline mesh generation.

Geometry:
  Flange A: centerline x in [0,60], z=0, thickness in z: [-2,2]
            (fixed at x=0)
  Flange B: centerline z in [0,60], x=60, thickness in x: [58,62]
            (free tip at z=60)
  Both: 40mm constant width (y), same material/section as Stages 1-3.

Built via two OCC boxes + BooleanUnion. Sharp reentrant (interior)
corner at (x=58, z=2) - deliberately no fillet, no gusset, in this
sub-stage. A stress singularity is expected at this reentrant corner
(a known, unavoidable feature of a zero-radius interior corner, not
a modeling defect) - quantitative comparison points are all chosen
away from it.

Outer corner comes out as a small unmitered step (boxes don't meet
at a 45 degree miter) rather than a single sharp point - a minor
construction detail, flagged in the session record, not expected to
affect the quantitative verification (all comparison points are on
the flanges, away from both corners).

Load: F=235.44N in the -x direction at Flange B's tip (z=60), kept
parallel to Flange A's axis so the whole problem stays planar (pure
bending in Flange B, combined bending+axial in Flange A) - matches
frame_analytical_reference.py's locked closed-form targets.
"""

import gmsh
import sys

TOL = 1e-6

# ---- Locked geometry parameters ----
L = 60.0           # mm, length of each flange (centerline)
WIDTH = 40.0       # mm, constant width (y)
THICK = 4.0        # mm, constant thickness
HALF_T = THICK / 2.0  # 2.0

# Flange A: x in [0,L], z in [-HALF_T, HALF_T]
# Flange B: x in [L-HALF_T, L+HALF_T], z in [0,L]

# ---- Baseline mesh sizing ----
FAR_SIZE = 3.0
CORNER_SIZE = 0.5     # first-pass local refinement near reentrant corner
TRANSITION_WIDTH = 4.0

OUT_MESH = "01-aerospace-mounting-bracket/mesh/frame_study/mesh_frame_raw.inp"

gmsh.initialize()
gmsh.model.add("frame_case")
occ = gmsh.model.occ

flange_a = occ.addBox(0.0, 0.0, -HALF_T, L, WIDTH, THICK)
flange_b = occ.addBox(L - HALF_T, 0.0, 0.0, THICK, WIDTH, L)
occ.synchronize()

print("=== PRE-FUSE SANITY CHECK ===")
print(f"Flange A bbox: {gmsh.model.getBoundingBox(3, flange_a)}  (expect x:0-60, y:0-40, z:-2-2)")
print(f"Flange B bbox: {gmsh.model.getBoundingBox(3, flange_b)}  (expect x:58-62, y:0-40, z:0-60)")

out, out_map = occ.fuse([(3, flange_a)], [(3, flange_b)])
occ.synchronize()
fused_vol = out[0][1]

print("\n=== POST-FUSE SANITY CHECK ===")
bb = gmsh.model.getBoundingBox(3, fused_vol)
print(f"Fused volume bbox: {bb}  (expect x:0-62, y:0-40, z:-2-60)")

# ---- Identify the reentrant corner edge (for local refinement) ----
# Expected at x=58 (L-HALF_T), z=2 (HALF_T), spanning full y width -
# a straight edge, not a curved surface this time.
all_curves = gmsh.model.getEntities(1)
corner_edge = None
for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    x_at_corner = abs(xmin - (L - HALF_T)) < TOL and abs(xmax - (L - HALF_T)) < TOL
    z_at_corner = abs(zmin - HALF_T) < TOL and abs(zmax - HALF_T) < TOL
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if x_at_corner and z_at_corner and spans_width:
        corner_edge = tag
        break

print(f"\n=== REENTRANT CORNER EDGE IDENTIFICATION ===")
print(f"Corner edge (x={L-HALF_T}, z={HALF_T}, spanning y): {corner_edge}")

if corner_edge is None:
    print("\nFATAL: could not identify the reentrant corner edge. Stopping.")
    gmsh.finalize()
    sys.exit(1)

# ---- Identify root face (x=0) and tip face (z=L, on flange B) ----
surfaces = gmsh.model.getEntities(2)
root_face_tag = None
tip_face_tag = None
for dim, tag in surfaces:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    if abs(xmin) < TOL and abs(xmax) < TOL:
        root_face_tag = tag
    if abs(zmin - L) < TOL and abs(zmax - L) < TOL:
        tip_face_tag = tag

print(f"\n=== ROOT / TIP FACE IDENTIFICATION ===")
print(f"Root face (x=0): {root_face_tag}")
print(f"Tip face (z={L}): {tip_face_tag}")

if root_face_tag is None or tip_face_tag is None:
    print("\nFATAL: could not identify root and/or tip face. Stopping.")
    gmsh.finalize()
    sys.exit(1)

gmsh.model.addPhysicalGroup(3, [fused_vol], name="frame")
gmsh.model.addPhysicalGroup(2, [root_face_tag], name="root")
gmsh.model.addPhysicalGroup(2, [tip_face_tag], name="tip")

# ---- Mesh size field: Distance from the reentrant corner edge,
# sigmoid Threshold (same fix already validated in Stage 3) ----
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "CurvesList", [corner_edge])
gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", CORNER_SIZE)
gmsh.model.mesh.field.setNumber(2, "SizeMax", FAR_SIZE)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
gmsh.model.mesh.field.setNumber(2, "DistMax", TRANSITION_WIDTH)
gmsh.model.mesh.field.setNumber(2, "Sigmoid", 1)
gmsh.model.mesh.field.setAsBackgroundMesh(2)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

print("\n=== GENERATING MESH ===")
gmsh.model.mesh.generate(3)
gmsh.write(OUT_MESH)
print(f"\nMesh written to: {OUT_MESH}")

n_nodes = len(gmsh.model.mesh.getNodes()[0])
elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
n_elements = sum(len(t) for t in elem_tags)
print(f"Nodes: {n_nodes}")
print(f"Volume elements: {n_elements}")

# ---- NSET export, same session ----
NSET_OUT = "01-aerospace-mounting-bracket/mesh/frame_study/nsets_frame.inp"
named_surfaces = ["root", "tip"]
groups = gmsh.model.getPhysicalGroups(dim=2)
name_to_tag = {gmsh.model.getPhysicalName(d, t): t for d, t in groups}
with open(NSET_OUT, "w") as f:
    for name in named_surfaces:
        tag = name_to_tag[name]
        node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
        node_tags = sorted(int(nd) for nd in node_tags)
        print(f"NSET '{name}': {len(node_tags)} nodes")
        f.write(f"*NSET,NSET={name}\n")
        for i in range(0, len(node_tags), 8):
            chunk = node_tags[i:i + 8]
            f.write(",".join(str(nd) for nd in chunk) + ",\n")
print(f"Wrote: {NSET_OUT}")

gmsh.finalize()
