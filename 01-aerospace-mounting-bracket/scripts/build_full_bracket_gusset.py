"""
Stage 5 — Full L-bracket integration.
Increment 2 (CORRECTED): add the gusset (quadrilateral cross-section,
Option A geometry) to the verified R3-filleted base frame.

CORRECTION vs prior attempt: the gusset cross-section spans the full
flange width (y=0 to y=40mm), NOT 4mm. GUSSET_THICKNESS=4mm was a
leftover descriptive label (matching flange plate thickness) and was
incorrectly used as the y-extrusion length in the previous version.
The extrusion length is now FRAME_WIDTH=40mm, independent of the
4mm plate-thickness value.

Does NOT yet include: toe fillet, hole. Those are separate increments.

Gusset cross-section vertices (x-z plane, extruded across full y-width):
    Toe  = (38, 2)
    At   = (55, 2)   -- tangent point, Flange A concave face
    Bt   = (58, 5)   -- tangent point, Flange B concave face
    Far2 = (58, 22)  -- unchanged from Stage 4b

Edge order for the closed loop: Toe -> At -> Bt -> Far2 -> Toe
"""

import gmsh
import sys
import os

BASE_FRAME_BREP = "mesh/full_bracket_study/base_frame_check.brep"
OUTPUT_BREP = "mesh/full_bracket_study/frame_plus_gusset_check.brep"

FRAME_WIDTH = 40.0  # mm, in y -- full flange width, CORRECT extrusion length
Y0 = 0.0
Y1 = FRAME_WIDTH

# Locked cross-section vertices (x, z)
TOE  = (38.0, 2.0)
AT   = (55.0, 2.0)
BT   = (58.0, 5.0)
FAR2 = (58.0, 22.0)

def main():
    if not os.path.exists(BASE_FRAME_BREP):
        print(f"ERROR: base frame file not found: {BASE_FRAME_BREP}")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("frame_plus_gusset")

    # --- Import verified base frame ---
    gmsh.model.occ.importShapes(BASE_FRAME_BREP)
    gmsh.model.occ.synchronize()

    base_volumes = gmsh.model.getEntities(dim=3)
    print(f"Base frame import: {len(base_volumes)} volume(s)")
    if len(base_volumes) != 1:
        print("ERROR: expected exactly 1 volume from base frame import.")
        sys.exit(1)
    base_vol_tag = base_volumes[0][1]
    base_mass = gmsh.model.occ.getMass(3, base_vol_tag)
    print(f"Base frame mass (pre-fuse): {base_mass:.4f} mm^3")

    # --- Build gusset cross-section as a 2D polygon at y=0 ---
    p1 = gmsh.model.occ.addPoint(TOE[0], Y0, TOE[1])
    p2 = gmsh.model.occ.addPoint(AT[0],  Y0, AT[1])
    p3 = gmsh.model.occ.addPoint(BT[0],  Y0, BT[1])
    p4 = gmsh.model.occ.addPoint(FAR2[0],Y0, FAR2[1])

    l1 = gmsh.model.occ.addLine(p1, p2)  # Toe -> At
    l2 = gmsh.model.occ.addLine(p2, p3)  # At  -> Bt
    l3 = gmsh.model.occ.addLine(p3, p4)  # Bt  -> Far2
    l4 = gmsh.model.occ.addLine(p4, p1)  # Far2 -> Toe

    loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()

    # --- Extrude across the FULL frame width (y=0 to y=40) ---
    extruded = gmsh.model.occ.extrude([(2, surf)], 0, FRAME_WIDTH, 0)
    gmsh.model.occ.synchronize()

    gusset_vol_tags = [e[1] for e in extruded if e[0] == 3]
    print(f"Gusset extrusion produced {len(gusset_vol_tags)} volume(s): {gusset_vol_tags}")
    if len(gusset_vol_tags) != 1:
        print("ERROR: expected exactly 1 gusset volume from extrusion.")
        sys.exit(1)
    gusset_vol_tag = gusset_vol_tags[0]
    gusset_mass = gmsh.model.occ.getMass(3, gusset_vol_tag)
    print(f"Gusset volume mass (pre-fuse, standalone): {gusset_mass:.4f} mm^3")
    expected_gusset_mass = 195.5 * FRAME_WIDTH
    print(f"Expected gusset mass (195.5mm^2 x {FRAME_WIDTH}mm): {expected_gusset_mass:.4f} mm^3")

    # --- Fuse gusset onto base frame ---
    fused, _ = gmsh.model.occ.fuse(
        [(3, base_vol_tag)], [(3, gusset_vol_tag)]
    )
    gmsh.model.occ.synchronize()

    print(f"Fuse result: {len(fused)} volume(s)")
    if len(fused) != 1:
        print("ERROR: fuse did not produce a single solid volume.")
        sys.exit(1)

    fused_vol_tag = fused[0][1]

    # --- Diagnostics ---
    mass = gmsh.model.occ.getMass(3, fused_vol_tag)
    mass_increase = mass - base_mass
    print(f"\nFused solid mass: {mass:.4f} mm^3")
    print(f"Mass increase over base frame: {mass_increase:.4f} mm^3")
    print(f"(Compare to expected gusset mass {expected_gusset_mass:.4f} mm^3 -- "
          f"should match closely if no overlap/gap with fillet region)")

    bbox = gmsh.model.occ.getBoundingBox(3, fused_vol_tag)
    print(f"\nBounding box: xmin={bbox[0]:.4f} ymin={bbox[1]:.4f} zmin={bbox[2]:.4f} "
          f"xmax={bbox[3]:.4f} ymax={bbox[4]:.4f} zmax={bbox[5]:.4f}")

    y_span_ok = abs(bbox[1] - 0.0) < 1e-6 and abs(bbox[4] - 40.0) < 1e-6
    print(f"Gusset spans full y-width (0 to 40mm)? "
          f"{'YES' if y_span_ok else 'NO -- CHECK'} "
          f"(ymin={bbox[1]:.4f}, ymax={bbox[4]:.4f})")

    # --- Check original corner enclosure across full width ---
    # The old Stage 4a reentrant corner ran along the edge x=58,z=2, y=0..40.
    # Sample points along this line and check they are now INSIDE the solid
    # (not on a boundary/edge), consistent with "fully enclosed."
    print("\nChecking original reentrant-corner line (x=58,z=2) is enclosed across y=0..40:")
    # Use a simple bounding-box/topology proxy: check whether the fused solid's
    # surfaces still contain the sharp x=58,z=2 edge anywhere. We do this by
    # inspecting curve entities for any curve lying exactly on that line.
    curves = gmsh.model.getEntities(dim=1)
    corner_line_found = False
    for (dim, tag) in curves:
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        # A curve exactly on x=58,z=2 would have bb xmin=xmax=58, zmin=zmax=2
        if abs(bb[0] - 58.0) < 1e-6 and abs(bb[3] - 58.0) < 1e-6 and \
           abs(bb[2] - 2.0) < 1e-6 and abs(bb[5] - 2.0) < 1e-6:
            corner_line_found = True
            print(f"  WARNING: curve {tag} still lies exactly on the old corner line "
                  f"(y-span {bb[1]:.4f} to {bb[4]:.4f})")
    if not corner_line_found:
        print("  No curve entity found lying exactly on the old x=58,z=2 corner line.")
        print("  Consistent with full enclosure across the width.")

    surfaces = gmsh.model.getEntities(dim=2)
    print(f"\nTotal surfaces on fused solid: {len(surfaces)}")
    print("Surface areas (sorted, smallest first):")
    areas = []
    for (dim, tag) in surfaces:
        area = gmsh.model.occ.getMass(2, tag)
        areas.append((area, tag))
    areas.sort()
    for area, tag in areas:
        print(f"  surface {tag}: area = {area:.4f} mm^2")

    slivers = [t for a, t in areas if a < 1.0]
    if slivers:
        print(f"\nWARNING: possible sliver surfaces (<1.0 mm^2): {slivers}")
    else:
        print("\nNo sliver surfaces detected (<1.0 mm^2 threshold).")

    # Duplicate surface check: any two surfaces with identical area AND bbox
    # (crude but useful proxy for accidental coincident-surface duplication)
    print("\nChecking for possible duplicate/coincident surfaces...")
    bboxes = {}
    dup_found = False
    for (dim, tag) in surfaces:
        bb = tuple(round(v, 4) for v in gmsh.model.occ.getBoundingBox(dim, tag))
        if bb in bboxes:
            print(f"  WARNING: surface {tag} has identical bbox to surface {bboxes[bb]}")
            dup_found = True
        else:
            bboxes[bb] = tag
    if not dup_found:
        print("  No duplicate/coincident surfaces detected.")

    # Manifold check: gmsh occ doesn't have a direct isManifold call, but
    # getMass succeeding on the volume + fuse returning exactly 1 volume
    # + all surfaces having finite positive area is the practical check here.
    non_manifold_flag = any(a <= 0 or a != a for a, t in areas)  # a!=a catches NaN
    print(f"\nNon-manifold/degenerate surface indicators (area<=0 or NaN): "
          f"{'FOUND -- PROBLEM' if non_manifold_flag else 'none found'}")

    # Save
    os.makedirs(os.path.dirname(OUTPUT_BREP), exist_ok=True)
    gmsh.write(OUTPUT_BREP)
    print(f"\nSaved fused geometry to {OUTPUT_BREP}")

    gmsh.finalize()

if __name__ == "__main__":
    main()
