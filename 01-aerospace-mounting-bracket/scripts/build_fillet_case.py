"""
Stage 3 (fillet) — geometry construction and baseline mesh generation.

Geometry: 60mm cantilever, 40mm constant width, symmetric opposite
shoulder fillet transition at x=30mm:
  - Root section (x=0-30mm):  depth D=6mm, z in [-3, +3]  (centroid z=0)
  - Tip section  (x=30-60mm): depth d=4mm,  z in [-2, +2]  (centroid z=0)
  - Fillet radius r=0.95mm on BOTH the top notch (x=30, z=+2) and
    bottom notch (x=30, z=-2) edges.

  r/d = 0.95/4 = 0.2375 (locked target was r=1.0mm, r/d=0.25 exactly;
  see history below for why r was reduced).

Mesh sizing: local refinement via a Distance+Threshold field pair
anchored on the two fillet (Cylinder-type) surfaces, identified
robustly by OCC surface type + z-sign — NOT a full-cross-section Box
field (an earlier attempt at that produced a 1.44M-element mesh).

CONSTRUCTION HISTORY (why r=0.95mm, not the originally locked
r=1.0mm) - kept here since it explains a real, non-obvious geometric
constraint that will recur if this pattern (thin riser + fillet) is
reused in a later stage:
  1. r=1.0mm exactly equals the riser height h=(D-d)/2=1.0mm - the
     fillet tangent point coincides with the adjacent sharp-corner
     vertex, a BRep degeneracy. OCC's fillet solver fails outright
     ("Could not compute fillet").
  2. r=0.99mm succeeds numerically but leaves a residual sliver face
     (h-r=0.01mm tall) smaller than the finest intended mesh size
     (0.15mm) - forces pathological element aspect ratios (up to
     AR=62.77) in its immediate neighborhood on the flat root-side
     face adjacent to the fillet (not on the fillet arc itself).
  3. gmsh.model.occ.healShapes() was tried to remove the sliver post-
     fillet. It did NOT repair the geometry in place - it duplicated
     the surface set (12 -> 22 surfaces, with the original un-healed
     topology left orphaned alongside a new one) and shifted the
     bounding box by ~0.015mm (comparable to the heal tolerance used).
     This silently corrupted downstream surface identification (it
     picked up stale, orphaned duplicate surfaces). Reverted entirely.
  4. r=0.95mm: sliver grows to h-r=0.05mm - large enough to mesh
     without pathological aspect ratios, no healing needed. This is
     the active configuration below.

This script is the authoritative, reproducible geometry+mesh source
for Stage 3. Does NOT modify any Stage 1 or Stage 2 files.
"""

import gmsh
import sys

TOL = 1e-6

# ---- Locked geometry parameters (Stage 3 specification, signed off) ----
LENGTH = 60.0
WIDTH  = 40.0
D_ROOT = 6.0
D_TIP  = 4.0
X_TRANS = 30.0
FILLET_R = 0.95

HALF_ROOT = D_ROOT / 2.0
HALF_TIP  = D_TIP / 2.0

# ---- Healing tolerance ----
HEAL_TOLERANCE = 0.05  # mm - above the 0.01mm sliver, below fillet R and mesh sizes

# ---- Baseline mesh sizing ----
FAR_SIZE    = 3.0
FILLET_SIZE = 0.15
TRANSITION_WIDTH = 1.5

OUT_MESH = "01-aerospace-mounting-bracket/mesh/fillet_study/mesh_fillet_raw.inp"

gmsh.initialize()
gmsh.model.add("fillet_case")

occ = gmsh.model.occ

# ---- Step 1: build the two boxes ----
root_box = occ.addBox(0.0, 0.0, -HALF_ROOT, X_TRANS, WIDTH, D_ROOT)
tip_box  = occ.addBox(X_TRANS, 0.0, -HALF_TIP, LENGTH - X_TRANS, WIDTH, D_TIP)
occ.synchronize()

print("=== PRE-FUSE SANITY CHECK ===")
print(f"Root box bbox: {gmsh.model.getBoundingBox(3, root_box)}  (expect x:0-30, y:0-40, z:-3-3)")
print(f"Tip  box bbox: {gmsh.model.getBoundingBox(3, tip_box)}  (expect x:30-60, y:0-40, z:-2-2)")

# ---- Step 2: boolean union ----
out, out_map = occ.fuse([(3, root_box)], [(3, tip_box)])
occ.synchronize()
fused_vol = out[0][1]

print("\n=== POST-FUSE SANITY CHECK ===")
print(f"Fused volume bbox: {gmsh.model.getBoundingBox(3, fused_vol)}  (expect x:0-60, y:0-40, z:-3-3)")

# ---- Step 3: identify the two reentrant (notch) edges to fillet ----
all_curves = gmsh.model.getEntities(1)
top_edges, bottom_edges, excluded_edges = [], [], []

for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    x_at_trans = abs(xmin - X_TRANS) < TOL and abs(xmax - X_TRANS) < TOL
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if not (x_at_trans and spans_width):
        continue
    z_mid = 0.5 * (zmin + zmax)
    if (zmax - zmin) >= TOL:
        continue
    if abs(z_mid - HALF_TIP) < TOL:
        top_edges.append(tag)
    elif abs(z_mid + HALF_TIP) < TOL:
        bottom_edges.append(tag)
    elif abs(z_mid - HALF_ROOT) < TOL or abs(z_mid + HALF_ROOT) < TOL:
        excluded_edges.append(tag)

print("\n=== FILLET EDGE IDENTIFICATION ===")
print(f"Top notch edges (to be filleted):    {top_edges}")
print(f"Bottom notch edges (to be filleted): {bottom_edges}")
print(f"Excluded outer-corner edges (NOT filleted): {excluded_edges}")

if not top_edges or not bottom_edges:
    print("\nFATAL: could not identify both notch edges. Stopping before fillet.")
    gmsh.finalize()
    sys.exit(1)

fillet_edges = top_edges + bottom_edges
radii = [FILLET_R] * len(fillet_edges)

# ---- Step 4: apply the fillet ----
filleted = occ.fillet([fused_vol], fillet_edges, radii, True)
occ.synchronize()
final_vol = filleted[0][1]

print("\n=== POST-FILLET SANITY CHECK ===")
print(f"Filleted volume bbox: {gmsh.model.getBoundingBox(3, final_vol)}  (expect x:0-60, y:0-40, z:-3-3, unchanged)")
n_surfaces_final = len(gmsh.model.getEntities(2))
print(f"Surface count: {n_surfaces_final} (expect 12; includes a 0.05mm sliver face at r=0.95mm - see docstring history)")

# ---- Step 5 (heal attempt REMOVED - see session record: it duplicated
# the geometry rather than fixing it, and shifted the bounding box by
# ~0.015mm. Reverted; FILLET_R=0.95mm used instead as the fallback.) ----

# ---- Step 6: identify the two fillet (Cylinder) surfaces + root face ----
surfaces = gmsh.model.getEntities(2)
top_fillet_surf = None
bottom_fillet_surf = None
root_face_tag = None

print("\n=== FULL SURFACE LIST (type + bbox, for verification) ===")
for dim, tag in sorted(surfaces, key=lambda s: s[1]):
    stype = gmsh.model.getType(dim, tag)
    bbox = gmsh.model.getBoundingBox(dim, tag)
    print(f"Surface {tag}: type={stype}  bbox={bbox}")

for dim, tag in surfaces:
    stype = gmsh.model.getType(dim, tag)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    if stype == "Cylinder":
        z_mid = 0.5 * (zmin + zmax)
        if z_mid > 0:
            top_fillet_surf = tag
        elif z_mid < 0:
            bottom_fillet_surf = tag
    if abs(xmin) < TOL and abs(xmax) < TOL:
        root_face_tag = tag

print("\n=== FILLET SURFACE IDENTIFICATION (for mesh sizing) ===")
print(f"Top fillet surface:    {top_fillet_surf}")
print(f"Bottom fillet surface: {bottom_fillet_surf}")
print(f"Root face (x=0):       {root_face_tag}")

if top_fillet_surf is None or bottom_fillet_surf is None:
    print("\nFATAL: could not identify both fillet surfaces. Stopping before meshing.")
    gmsh.finalize()
    sys.exit(1)

if root_face_tag is None:
    print("\nFATAL: could not identify root face (x=0). Stopping.")
    gmsh.finalize()
    sys.exit(1)

# ---- Step 7: physical groups ----
gmsh.model.addPhysicalGroup(3, [final_vol], name="plate")
gmsh.model.addPhysicalGroup(2, [root_face_tag], name="root")

# ---- Step 8: mesh size field ----
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [top_fillet_surf, bottom_fillet_surf])
gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", FILLET_SIZE)
gmsh.model.mesh.field.setNumber(2, "SizeMax", FAR_SIZE)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
gmsh.model.mesh.field.setNumber(2, "DistMax", TRANSITION_WIDTH)
gmsh.model.mesh.field.setNumber(2, "Sigmoid", 1)
gmsh.model.mesh.field.setAsBackgroundMesh(2)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# ---- Step 9: mesh settings ----
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

# ---- Step 10: NSET export, same session (guarantees numbering
# consistency - no reopen needed, unlike the Stage 2 .geo-reopen
# approach, since this geometry has no persisted .geo equivalent) ----
gmsh.model.addPhysicalGroup(2, [top_fillet_surf], name="fillet_top")
gmsh.model.addPhysicalGroup(2, [bottom_fillet_surf], name="fillet_bottom")

NSET_OUT = "01-aerospace-mounting-bracket/mesh/fillet_study/nsets_fillet.inp"
named_surfaces = ["root", "fillet_top", "fillet_bottom"]
groups = gmsh.model.getPhysicalGroups(dim=2)
name_to_tag = {gmsh.model.getPhysicalName(d, t): t for d, t in groups}

with open(NSET_OUT, "w") as f:
    for name in named_surfaces:
        tag = name_to_tag[name]
        node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
        node_tags = sorted(int(n) for n in node_tags)
        print(f"NSET '{name}': {len(node_tags)} nodes")
        f.write(f"*NSET,NSET={name}\n")
        for i in range(0, len(node_tags), 8):
            chunk = node_tags[i:i + 8]
            f.write(",".join(str(n) for n in chunk) + ",\n")

print(f"Wrote: {NSET_OUT}")

gmsh.finalize()
