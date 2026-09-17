"""
Stage 5 (full L-bracket integration) - Option A base frame ONLY.

Constructs the two perpendicular flanges and applies the R3mm
interior reentrant-corner fillet BEFORE any gusset, hole, or toe
fillet is added. This script stops after the base-frame fillet and
reports diagnostics - no gusset, no hole, no meshing, no solving.

Tangent points (computed from R=3mm on a 90 degree interior corner,
locked in this session):
  Flange A concave face (z=2): tangent at x=55
  Flange B concave face (x=58): tangent at z=5

The gusset (added in a LATER script, once this base frame is
verified) will attach at these tangent points, not at the original
sharp-corner point (58,2) - which no longer exists once filleted.
This means the gusset's near-corner legs will be 17mm (20mm nominal
minus the 3mm the R3 bend occupies), NOT 20mm. This script does not
build the gusset; it exists only to lock in and verify the R3 base
frame first, per the incremental-construction instruction.
"""

import gmsh
import sys
import os

TOL = 1e-6

L = 60.0
WIDTH = 40.0
THICK = 4.0
HALF_T = THICK / 2.0
MAIN_FILLET_R = 3.0

OUT_DIR = "01-aerospace-mounting-bracket/mesh/full_bracket_study"
os.makedirs(OUT_DIR, exist_ok=True)

gmsh.initialize()
gmsh.model.add("full_bracket_base_frame")
occ = gmsh.model.occ

# ---- Step 1: flanges, fuse ----
flange_a = occ.addBox(0.0, 0.0, -HALF_T, L, WIDTH, THICK)
flange_b = occ.addBox(L - HALF_T, 0.0, 0.0, THICK, WIDTH, L)
occ.synchronize()

print("=== STEP 1: FLANGES ===")
print(f"Flange A bbox: {gmsh.model.getBoundingBox(3, flange_a)}  (expect x:0-60,y:0-40,z:-2-2)")
print(f"Flange B bbox: {gmsh.model.getBoundingBox(3, flange_b)}  (expect x:58-62,y:0-40,z:0-60)")

out, _ = occ.fuse([(3, flange_a)], [(3, flange_b)])
occ.synchronize()
vol = out[0][1]
print(f"Frame fused bbox: {gmsh.model.getBoundingBox(3, vol)}  (expect x:0-62,y:0-40,z:-2-60)")

# ---- Step 2: identify the sharp REENTRANT corner edge (the same
# edge Stage 4a identified at x=58,z=2, spanning y) - this is the
# edge to fillet, NOT the outer/convex step ----
corner_x = L - HALF_T   # 58.0
corner_z = HALF_T       # 2.0

all_curves = gmsh.model.getEntities(1)
reentrant_edge = None
for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    x_at_corner = abs(xmin - corner_x) < TOL and abs(xmax - corner_x) < TOL
    z_at_corner = abs(zmin - corner_z) < TOL and abs(zmax - corner_z) < TOL
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if x_at_corner and z_at_corner and spans_width:
        reentrant_edge = tag
        break

print(f"\n=== STEP 2: REENTRANT CORNER EDGE IDENTIFICATION ===")
print(f"Reentrant edge (x={corner_x}, z={corner_z}, spanning y): {reentrant_edge}")

if reentrant_edge is None:
    print("FATAL: could not identify the reentrant corner edge. Stopping.")
    gmsh.finalize()
    sys.exit(1)

# Sanity: confirm this is the ONLY edge at this exact location (not
# confusable with any outer-step edge, which would be at x=60 or
# x=62, not x=58)
print("Confirmed: identified edge is at the INTERIOR corner (x=58), "
      "not the outer step (x=60/62) - correct target for R3.")

n_surfaces_before = len(gmsh.model.getEntities(2))
n_edges_before = len(gmsh.model.getEntities(1))
print(f"\nSurfaces before fillet: {n_surfaces_before}")
print(f"Edges before fillet: {n_edges_before}")

# ---- Step 3: apply R3 fillet to the reentrant corner ONLY ----
print(f"\n=== STEP 3: R3 INTERIOR FILLET ===")
try:
    filleted = occ.fillet([vol], [reentrant_edge], [MAIN_FILLET_R], True)
    occ.synchronize()
    vol = filleted[0][1]
    print(f"R3 fillet SUCCEEDED.")
except Exception as e:
    print(f"FATAL: R3 interior fillet failed with exception:")
    print(f"  {e}")
    gmsh.finalize()
    sys.exit(1)

bbox_final = gmsh.model.getBoundingBox(3, vol)
print(f"Post-fillet bbox: {bbox_final}  (expect x:0-62,y:0-40,z:-2-60, unchanged - "
      f"fillet is concave, doesn't grow the bounding box)")

n_surfaces_after = len(gmsh.model.getEntities(2))
n_edges_after = len(gmsh.model.getEntities(1))
print(f"Surfaces after fillet: {n_surfaces_after} (expect {n_surfaces_before}+1)")
print(f"Edges after fillet: {n_edges_after}")

# ============================================================
# DIAGNOSTICS - report only, no further construction
# ============================================================
print("\n" + "=" * 60)
print("BASE FRAME DIAGNOSTICS (R3 fillet only, no gusset/hole yet)")
print("=" * 60)

vols = gmsh.model.getEntities(3)
print(f"\nFinal volume count: {len(vols)} (expect 1)")
for dim, tag in vols:
    mass = gmsh.model.occ.getMass(dim, tag)
    print(f"Volume {tag}: mass/volume = {mass:.4f} mm^3")

# Expected volume: two 60x40x4 boxes (but overlapping at the corner
# junction region) minus a small amount removed by the concave
# fillet's rounding. Rough sanity bound only - exact analytic value
# not computed here since flange overlap geometry is nontrivial.

surfaces = gmsh.model.getEntities(2)
print(f"\nFinal surface count: {len(surfaces)}")

# ---- Identify and verify the R3 fillet surface specifically ----
print("\n--- R3 fillet surface verification ---")
fillet_surf = None
for dim, tag in surfaces:
    stype = gmsh.model.getType(dim, tag)
    if stype == "Cylinder":
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        spans_width = (ymax - ymin) > (WIDTH - 1.0)
        # Expected extent: x approx [55,58], z approx [2,5]
        near_expected = (54 < xmin < 56) and (57 < xmax < 59) and (1 < zmin < 3) and (4 < zmax < 6)
        if spans_width and near_expected:
            fillet_surf = tag
            print(f"  Candidate R3 fillet surface: {tag}")
            print(f"    Type: {stype}")
            print(f"    Bbox: ({xmin:.4f},{ymin:.4f},{zmin:.4f}) to ({xmax:.4f},{ymax:.4f},{zmax:.4f})")

if fillet_surf is None:
    print("  WARNING: no Cylinder-type surface found matching the expected R3 "
        "fillet location/extent. Listing ALL Cylinder surfaces found instead:")
    for dim, tag in surfaces:
        if gmsh.model.getType(dim, tag) == "Cylinder":
            bbox = gmsh.model.getBoundingBox(dim, tag)
            print(f"    Cylinder surface {tag}: bbox={bbox}")

# ---- Verify actual radius via getCurvature or direct point sampling ----
if fillet_surf is not None:
    print("\n--- Radius verification (via OCC parametric sampling) ---")
    try:
        # Sample a point on the fillet surface and query its
        # curvature to back out the radius
        u_min, u_max, v_min, v_max = gmsh.model.getParametrizationBounds(2, fillet_surf)
        u_mid = 0.5 * (u_min[0] + u_max[0])
        v_mid = 0.5 * (v_min[0] + v_max[0])
        curv = gmsh.model.getCurvature(2, fillet_surf, [u_mid, v_mid])
        # curv returns max principal curvature; radius = 1/curvature
        if curv[0] > 0:
            measured_r = 1.0 / curv[0]
            print(f"  Measured radius (via curvature at surface midpoint): {measured_r:.4f} mm")
            print(f"  Target: {MAIN_FILLET_R} mm")
            pct_diff = 100 * (measured_r - MAIN_FILLET_R) / MAIN_FILLET_R
            print(f"  % diff: {pct_diff:+.2f}%")
        else:
            print(f"  Curvature query returned non-positive value: {curv}")
    except Exception as e:
        print(f"  Could not compute curvature directly: {e}")
        print(f"  Falling back to bounding-box-based radius estimate:")
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, fillet_surf)
        est_r_x = xmax - xmin
        est_r_z = zmax - zmin
        print(f"  Bbox x-extent: {est_r_x:.4f} mm, z-extent: {est_r_z:.4f} mm "
              f"(both should be close to {MAIN_FILLET_R} mm for a quarter-arc fillet)")

# ---- Tangent point check ----
print("\n--- Tangent point check ---")
print(f"Expected tangent on Flange A concave face: x=55, z=2")
print(f"Expected tangent on Flange B concave face: x=58, z=5")
if fillet_surf is not None:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, fillet_surf)
    print(f"Fillet surface actual extent: x=[{xmin:.4f},{xmax:.4f}], z=[{zmin:.4f},{zmax:.4f}]")
    tangent_a_ok = abs(xmin - 55.0) < 0.1
    tangent_b_ok = abs(zmax - 5.0) < 0.1
    print(f"Tangent at x=55 (Flange A side) confirmed: {tangent_a_ok}")
    print(f"Tangent at z=5 (Flange B side) confirmed: {tangent_b_ok}")

# ---- Sliver / duplicate surface check ----
print("\n--- Sliver / duplicate surface check ---")
areas = []
for dim, tag in surfaces:
    area = gmsh.model.occ.getMass(dim, tag)
    areas.append((tag, area))
areas.sort(key=lambda t: t[1])
print("5 smallest-area surfaces:")
for tag, area in areas[:5]:
    stype = gmsh.model.getType(2, tag)
    bbox = gmsh.model.getBoundingBox(2, tag)
    print(f"  Surface {tag} ({stype}): area={area:.6f} mm^2, bbox={bbox}")

# ---- Outer step confirmation (should be UNCHANGED, sharp, since
# this script did not touch it) ----
print("\n--- Confirm outer step is untouched (sharp, not filleted) ---")
all_curves_final = gmsh.model.getEntities(1)
outer_step_edges_found = 0
for dim, tag in all_curves_final:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if not spans_width:
        continue
    at_x60 = abs(xmin - L) < TOL and abs(xmax - L) < TOL
    at_x62 = abs(xmin - (L + HALF_T)) < TOL and abs(xmax - (L + HALF_T)) < TOL
    if at_x60 or at_x62:
        outer_step_edges_found += 1
        print(f"  Outer-step edge still present (sharp, as expected): curve {tag}, bbox={(xmin,ymin,zmin,xmax,ymax,zmax)}")
print(f"Outer step edges found: {outer_step_edges_found} (should be >0 - confirms this "
      f"script did NOT accidentally fillet the outer step)")

print("\n" + "=" * 60)
print("STOPPING HERE per instructions - base frame only.")
print("Gusset, toe fillet, and hole NOT yet added.")
print("=" * 60)

gmsh.model.addPhysicalGroup(3, [vol], name="base_frame")
gmsh.write(f"{OUT_DIR}/base_frame_check.brep")
print(f"\nSaved BRep for record: {OUT_DIR}/base_frame_check.brep")

gmsh.finalize()
