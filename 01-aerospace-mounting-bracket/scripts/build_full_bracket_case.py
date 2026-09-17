"""
Stage 5 (full L-bracket integration) - geometry construction and
diagnostics ONLY. Does NOT mesh, solve, or optimize.

Combines, onto the Stage 4a/4b right-angle frame:
  - R3mm fillet at the OUTER (convex) flange-junction step (the
    "unmitered step" flagged in Stage 4a - the gusset occupies the
    INNER/concave corner, so this is the only remaining interpretation
    of "flange-junction fillet" that doesn't conflict with it).
  - 20x20x4mm concave-side triangular gusset (identical to Stage 4b).
  - R1mm fillet at the gusset toe (the Stage 4b singularity location)
    - first construction attempt; if this fails, STOP and report
      exactly why, do not silently change the radius.
  - Ø5mm through-hole in Flange A at x=20, y=20 (matches Stage 2's
    diameter and placement logic).

KNOWN RISK, flagged before running: the outer step is only 2mm tall
(both edges are 2mm long, set by wall thickness). R3 > step height is
a similar category of risk to Stage 3's r=h degeneracy - reported
honestly below, not glossed over.

Construction order: flanges -> fuse -> gusset -> fuse -> toe fillet ->
outer-step fillet -> hole cut. Each stage is checked before
proceeding to the next.
"""

import gmsh
import sys
import os

TOL = 1e-6

L = 60.0
WIDTH = 40.0
THICK = 4.0
HALF_T = THICK / 2.0
GUSSET_LEG = 20.0
TOE_FILLET_R = 1.0
MAIN_FILLET_R = 3.0
HOLE_D = 5.0
HOLE_X = 20.0
HOLE_Y = 20.0

OUT_DIR = "01-aerospace-mounting-bracket/mesh/full_bracket_study"
os.makedirs(OUT_DIR, exist_ok=True)

gmsh.initialize()
gmsh.model.add("full_bracket_case")
occ = gmsh.model.occ

# ---- Step 1: flanges, fuse ----
flange_a = occ.addBox(0.0, 0.0, -HALF_T, L, WIDTH, THICK)
flange_b = occ.addBox(L - HALF_T, 0.0, 0.0, THICK, WIDTH, L)
occ.synchronize()

print("=== STEP 1: FLANGES ===")
print(f"Flange A bbox: {gmsh.model.getBoundingBox(3, flange_a)}")
print(f"Flange B bbox: {gmsh.model.getBoundingBox(3, flange_b)}")

out, _ = occ.fuse([(3, flange_a)], [(3, flange_b)])
occ.synchronize()
vol = out[0][1]
print(f"Frame fused bbox: {gmsh.model.getBoundingBox(3, vol)}  (expect x:0-62,y:0-40,z:-2-60)")

# ---- Step 2: gusset, fuse ----
corner_x = L - HALF_T   # 58.0
corner_z = HALF_T       # 2.0

p1 = occ.addPoint(corner_x, 0.0, corner_z)
p2 = occ.addPoint(corner_x - GUSSET_LEG, 0.0, corner_z)
p3 = occ.addPoint(corner_x, 0.0, corner_z + GUSSET_LEG)
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

print(f"\n=== STEP 2: GUSSET ===")
print(f"Gusset bbox: {gmsh.model.getBoundingBox(3, gusset_vol)}  (expect x:38-58,y:0-40,z:2-22)")

out, _ = occ.fuse([(3, vol)], [(3, gusset_vol)])
occ.synchronize()
vol = out[0][1]
print(f"Frame+gusset fused bbox: {gmsh.model.getBoundingBox(3, vol)}")

n_surfaces_pre_fillets = len(gmsh.model.getEntities(2))
n_edges_pre_fillets = len(gmsh.model.getEntities(1))
print(f"Surfaces before any fillet: {n_surfaces_pre_fillets}")
print(f"Edges before any fillet: {n_edges_pre_fillets}")

# ---- Step 3: gusset toe fillet, r=1mm ----
# Toe edge: x = corner_x - GUSSET_LEG = 38, z = corner_z = 2, spanning y
all_curves = gmsh.model.getEntities(1)
toe_edge = None
for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    x_at_toe = abs(xmin - (corner_x - GUSSET_LEG)) < TOL and abs(xmax - (corner_x - GUSSET_LEG)) < TOL
    z_at_toe = abs(zmin - corner_z) < TOL and abs(zmax - corner_z) < TOL
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if x_at_toe and z_at_toe and spans_width:
        toe_edge = tag
        break

print(f"\n=== STEP 3: GUSSET TOE FILLET (R={TOE_FILLET_R}mm) ===")
print(f"Toe edge identified: {toe_edge}")

if toe_edge is None:
    print("FATAL: could not identify the gusset toe edge before filleting. Stopping.")
    gmsh.finalize()
    sys.exit(1)

try:
    filleted = occ.fillet([vol], [toe_edge], [TOE_FILLET_R], True)
    occ.synchronize()
    vol = filleted[0][1]
    print(f"Toe fillet SUCCEEDED. Post-toe-fillet bbox: {gmsh.model.getBoundingBox(3, vol)}")
except Exception as e:
    print(f"FATAL: gusset toe fillet (R={TOE_FILLET_R}mm) failed with exception:")
    print(f"  {e}")
    print("Stopping. Radius NOT silently changed - awaiting explicit decision.")
    gmsh.finalize()
    sys.exit(1)

n_surfaces_after_toe = len(gmsh.model.getEntities(2))
print(f"Surfaces after toe fillet: {n_surfaces_after_toe}")

# ---- Step 4: outer (convex) step fillet, r=3mm ----
# Two edges expected: vertical riser at x=60 (z:-2 to 0), and
# horizontal ledge at z=0 (x:60 to 62), both spanning y.
print("\n--- FULL CURVE DUMP (post-toe-fillet) for diagnosis ---")
all_curves_dump = gmsh.model.getEntities(1)
for dim, tag in sorted(all_curves_dump, key=lambda c: c[1]):
    bbox = gmsh.model.getBoundingBox(dim, tag)
    ctype = gmsh.model.getType(dim, tag)
    print(f"  Curve {tag} ({ctype}): bbox={bbox}")

# CORRECTED per diagnostic curve dump: the outer step is bounded by
# THREE edges, not two - a small riser FACE (x=60, z:-2 to 0) and a
# small ledge FACE (z=0, x:60 to 62) each contribute a boundary edge,
# plus a third edge where the ledge meets Flange B's outer face.
all_curves = gmsh.model.getEntities(1)
edge_bottom_to_riser = None    # x=60, z=-2, spans y
edge_riser_to_ledge = None     # x=60, z=0, spans y
edge_ledge_to_outer = None     # x=62, z=0, spans y
for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if not spans_width:
        continue
    if abs(xmin - L) < TOL and abs(xmax - L) < TOL and abs(zmin - (-HALF_T)) < TOL and abs(zmax - (-HALF_T)) < TOL:
        edge_bottom_to_riser = tag
    if abs(xmin - L) < TOL and abs(xmax - L) < TOL and abs(zmin - 0.0) < TOL and abs(zmax - 0.0) < TOL:
        edge_riser_to_ledge = tag
    if abs(xmin - (L + HALF_T)) < TOL and abs(xmax - (L + HALF_T)) < TOL and abs(zmin - 0.0) < TOL and abs(zmax - 0.0) < TOL:
        edge_ledge_to_outer = tag

print(f"\n=== STEP 4: OUTER STEP FILLET (R={MAIN_FILLET_R}mm) ===")
print(f"Edge bottom-face -> riser-face (x=60,z=-2): {edge_bottom_to_riser}")
print(f"Edge riser-face -> ledge-face (x=60,z=0): {edge_riser_to_ledge}")
print(f"Edge ledge-face -> outer-face (x=62,z=0): {edge_ledge_to_outer}")

fillet_edges = [e for e in [edge_bottom_to_riser, edge_riser_to_ledge, edge_ledge_to_outer] if e is not None]
print(f"Attempting fillet on {len(fillet_edges)} edge(s) at R={MAIN_FILLET_R}mm.")
if not fillet_edges:
    print("FATAL: could not identify either outer-step edge before filleting. Stopping.")
    gmsh.finalize()
    sys.exit(1)
if len(fillet_edges) < 2:
    print("WARNING: only one of the two expected outer-step edges was found. "
          "Proceeding with what was found, but this may indicate a geometry "
          "or identification problem - reported for review, not hidden.")

try:
    filleted = occ.fillet([vol], fillet_edges, [MAIN_FILLET_R] * len(fillet_edges), True)
    occ.synchronize()
    vol = filleted[0][1]
    print(f"Outer-step fillet SUCCEEDED. Post-main-fillet bbox: {gmsh.model.getBoundingBox(3, vol)}")
except Exception as e:
    print(f"FATAL: outer-step fillet (R={MAIN_FILLET_R}mm) failed with exception:")
    print(f"  {e}")
    print("Stopping. Radius NOT silently changed - awaiting explicit decision.")
    gmsh.finalize()
    sys.exit(1)

n_surfaces_after_main = len(gmsh.model.getEntities(2))
print(f"Surfaces after outer-step fillet: {n_surfaces_after_main}")

# ---- Step 5: hole cut, d=5mm at x=20,y=20, through Flange A (z) ----
hole_r = HOLE_D / 2.0
hole_cyl = occ.addCylinder(HOLE_X, HOLE_Y, -HALF_T - 1.0, 0.0, 0.0, THICK + 2.0, hole_r)
occ.synchronize()
print(f"\n=== STEP 5: HOLE CUT (D={HOLE_D}mm at x={HOLE_X}, y={HOLE_Y}) ===")
print(f"Hole cutting cylinder bbox: {gmsh.model.getBoundingBox(3, hole_cyl)}")

try:
    cut_out, _ = occ.cut([(3, vol)], [(3, hole_cyl)])
    occ.synchronize()
    vol = cut_out[0][1]
    print(f"Hole cut SUCCEEDED. Final bbox: {gmsh.model.getBoundingBox(3, vol)}")
except Exception as e:
    print(f"FATAL: hole cut failed with exception:")
    print(f"  {e}")
    gmsh.finalize()
    sys.exit(1)

# ============================================================
# GEOMETRY DIAGNOSTICS (per explicit instruction - report, do
# not proceed to meshing)
# ============================================================
print("\n" + "=" * 60)
print("GEOMETRY DIAGNOSTICS")
print("=" * 60)

# BRep validity
try:
    is_valid = gmsh.model.occ.getEntities(3)
    print(f"\nFinal volume count: {len(is_valid)} (expect 1 - single connected solid)")
except Exception as e:
    print(f"Could not query volume entities: {e}")

vols = gmsh.model.getEntities(3)
for dim, tag in vols:
    mass = gmsh.model.occ.getMass(dim, tag)
    print(f"Volume {tag}: mass/volume = {mass:.4f} mm^3")

final_surfaces = gmsh.model.getEntities(2)
final_curves = gmsh.model.getEntities(1)
print(f"\nFinal surface count: {len(final_surfaces)}")
print(f"Final curve count: {len(final_curves)}")

# Check for sliver/tiny-area surfaces
print("\n--- Surface area check (looking for slivers) ---")
areas = []
for dim, tag in final_surfaces:
    area = gmsh.model.occ.getMass(dim, tag)
    areas.append((tag, area))
areas.sort(key=lambda t: t[1])
print("5 smallest-area surfaces (checking for slivers):")
for tag, area in areas[:5]:
    stype = gmsh.model.getType(2, tag)
    bbox = gmsh.model.getBoundingBox(2, tag)
    print(f"  Surface {tag} ({stype}): area={area:.6f} mm^2, bbox={bbox}")

# Confirm fillet surfaces exist at expected locations/radii
print("\n--- Fillet surface confirmation ---")
found_toe_fillet = False
found_main_fillet = False
found_hole = False
for dim, tag in final_surfaces:
    stype = gmsh.model.getType(dim, tag)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    if stype == "Cylinder":
        # Toe fillet: near x=38, z=2, small extent ~TOE_FILLET_R
        if abs(xmin - (corner_x - GUSSET_LEG - TOE_FILLET_R)) < 0.5 or \
           (37 < xmin < 40 and 1 < zmin < 4):
            print(f"  Possible TOE fillet surface: {tag}, bbox={(xmin,ymin,zmin,xmax,ymax,zmax)}")
            found_toe_fillet = True
        # Main fillet: near x=60, z=0
        if 57 < xmin < 63 and -3 < zmin < 3:
            print(f"  Possible MAIN (outer) fillet surface: {tag}, bbox={(xmin,ymin,zmin,xmax,ymax,zmax)}")
            found_main_fillet = True
        # Hole: near x=20,y=20, radius 2.5
        if abs(xmin - (HOLE_X - hole_r)) < 0.5 and abs(xmax - (HOLE_X + hole_r)) < 0.5:
            print(f"  Possible HOLE surface: {tag}, bbox={(xmin,ymin,zmin,xmax,ymax,zmax)}")
            found_hole = True

print(f"\nToe fillet surface found: {found_toe_fillet}")
print(f"Main (outer) fillet surface found: {found_main_fillet}")
print(f"Hole surface found: {found_hole}")

# Dimension checks
overall_bbox = gmsh.model.getBoundingBox(3, vol)
print(f"\n--- Overall dimension check ---")
print(f"Final bbox: {overall_bbox}")
print(f"Expected: x approx [0,62], y=[0,40], z approx [-2,60] (allowing for R3 fillet slightly rounding the x=62/z=0 corner)")

print("\n" + "=" * 60)
print("STOPPING HERE per instructions - no meshing, solving, or")
print("optimization until geometry is explicitly signed off.")
print("=" * 60)

gmsh.finalize()
