"""
Stage 4b (right-angle frame WITH gusset) - geometry construction and
baseline mesh generation.

Adds a triangular-prism gusset to the exact Stage 4a frame geometry,
bridging the concave (inside-of-bend) faces:
  - Leg 1: along Flange A's concave face (z=2), from the corner
    (x=58) back to x=38 (20mm leg) - flush/coplanar with Flange A's
    top surface.
  - Leg 2: along Flange B's concave face (x=58), from the corner
    (z=2) up to z=22 (20mm leg) - flush/coplanar with Flange B's
    inner surface.
  - Hypotenuse: connects (38,y,2) to (58,y,22) directly - the only
    NEW surface feature this geometry introduces, since both legs
    are coplanar with existing flange faces and disappear into them
    after the union.
  - Full 40mm width (y), matching both flanges.

Expected consequence (to be confirmed against the actual meshed
geometry, not assumed): the original Stage 4a reentrant corner at
(58, y, 2) becomes fully enclosed inside the solid (the gusset's own
right-angle vertex sits there) and should no longer exist as a
surface feature - i.e., the known singularity from 4a should
disappear. The hypotenuse edge and its two endpoints are the new
region of interest for local refinement and stress inspection.

Same material, load (F=235.44N at Flange B tip, -x direction),
coordinate system, and root BC as Stage 4a.
"""

import gmsh
import sys

TOL = 1e-6

L = 60.0
WIDTH = 40.0
THICK = 4.0
HALF_T = THICK / 2.0

GUSSET_LEG = 20.0  # mm

FAR_SIZE = 3.0
CORNER_SIZE = 0.5
TRANSITION_WIDTH = 4.0

OUT_MESH = "01-aerospace-mounting-bracket/mesh/frame_gusset_study/mesh_frame_gusset_raw.inp"

gmsh.initialize()
gmsh.model.add("frame_gusset_case")
occ = gmsh.model.occ

# ---- Flanges (identical to Stage 4a) ----
flange_a = occ.addBox(0.0, 0.0, -HALF_T, L, WIDTH, THICK)
flange_b = occ.addBox(L - HALF_T, 0.0, 0.0, THICK, WIDTH, L)
occ.synchronize()

print("=== PRE-FUSE SANITY CHECK (flanges) ===")
print(f"Flange A bbox: {gmsh.model.getBoundingBox(3, flange_a)}  (expect x:0-60, y:0-40, z:-2-2)")
print(f"Flange B bbox: {gmsh.model.getBoundingBox(3, flange_b)}  (expect x:58-62, y:0-40, z:0-60)")

# ---- Gusset triangular prism ----
corner_x = L - HALF_T   # 58.0
corner_z = HALF_T       # 2.0

p1 = occ.addPoint(corner_x, 0.0, corner_z)                          # (58, 0, 2)
p2 = occ.addPoint(corner_x - GUSSET_LEG, 0.0, corner_z)             # (38, 0, 2)
p3 = occ.addPoint(corner_x, 0.0, corner_z + GUSSET_LEG)             # (58, 0, 22)

l1 = occ.addLine(p1, p2)
l2 = occ.addLine(p2, p3)
l3 = occ.addLine(p3, p1)
wire = occ.addCurveLoop([l1, l2, l3])
gusset_face = occ.addPlaneSurface([wire])
occ.synchronize()

gusset_extrude = occ.extrude([(2, gusset_face)], 0.0, WIDTH, 0.0)
occ.synchronize()

gusset_vol = None
for dim, tag in gusset_extrude:
    if dim == 3:
        gusset_vol = tag
        break

if gusset_vol is None:
    print("\nFATAL: gusset extrusion did not produce a volume. Stopping.")
    gmsh.finalize()
    sys.exit(1)

print("\n=== PRE-FUSE SANITY CHECK (gusset) ===")
print(f"Gusset bbox: {gmsh.model.getBoundingBox(3, gusset_vol)}  (expect x:38-58, y:0-40, z:2-22)")

# ---- Three-body fuse ----
out, out_map = occ.fuse([(3, flange_a)], [(3, flange_b), (3, gusset_vol)])
occ.synchronize()
fused_vol = out[0][1]

print("\n=== POST-FUSE SANITY CHECK ===")
bb = gmsh.model.getBoundingBox(3, fused_vol)
print(f"Fused volume bbox: {bb}  (expect x:0-62, y:0-40, z:-2-60)")

# ---- Identify the hypotenuse FACE (not an edge - a straight edge
# cannot simultaneously span the full width AND be diagonal in x-z;
# only the flat sloped FACE has both properties, same correction
# Stage 3 needed when it moved from edges to surfaces for the fillet) ----
all_surfaces_search = gmsh.model.getEntities(2)
hyp_face = None
for dim, tag in all_surfaces_search:
    stype = gmsh.model.getType(dim, tag)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    x_diagonal = abs(xmin - (corner_x - GUSSET_LEG)) < TOL and abs(xmax - corner_x) < TOL
    z_diagonal = abs(zmin - corner_z) < TOL and abs(zmax - (corner_z + GUSSET_LEG)) < TOL
    if stype == "Plane" and spans_width and x_diagonal and z_diagonal:
        hyp_face = tag
        break

print(f"\n=== HYPOTENUSE FACE IDENTIFICATION ===")
print(f"Hypotenuse face (x:38-58, z:2-22, spanning y, Plane type): {hyp_face}")

if hyp_face is None:
    print("\nFATAL: could not identify the hypotenuse face. Stopping.")
    gmsh.finalize()
    sys.exit(1)

all_curves = gmsh.model.getEntities(1)
# ---- Confirm the OLD reentrant-corner edge no longer exists as a
# distinct surface feature (expected consequence, checked not assumed) ----
old_corner_found = False
for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    x_at_corner = abs(xmin - corner_x) < TOL and abs(xmax - corner_x) < TOL
    z_at_corner = abs(zmin - corner_z) < TOL and abs(zmax - corner_z) < TOL
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if x_at_corner and z_at_corner and spans_width and tag != hyp_edge:
        old_corner_found = True
        print(f"NOTE: an edge still exists at the old corner location (tag {tag}) - "
              f"the singularity may NOT have been fully eliminated. Investigate before assuming otherwise.")

if not old_corner_found:
    print("Confirmed: no distinct edge remains at the old reentrant-corner location "
          "(58, z=2) - consistent with the predicted enclosure by the gusset.")

# ---- Root / tip faces ----
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

gmsh.model.addPhysicalGroup(3, [fused_vol], name="frame_gusset")
gmsh.model.addPhysicalGroup(2, [root_face_tag], name="root")
gmsh.model.addPhysicalGroup(2, [tip_face_tag], name="tip")

# ---- Mesh size field: Distance from the hypotenuse edge, sigmoid Threshold ----
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [hyp_face])
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
NSET_OUT = "01-aerospace-mounting-bracket/mesh/frame_gusset_study/nsets_frame_gusset.inp"
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
