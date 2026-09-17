"""
Stage 5 — Full L-bracket integration.
Increment 3: apply the R1mm gusset-toe fillet to the verified
frame+gusset solid (frame_plus_gusset_check.brep).

Toe edge location: x=38, z=2, running the full y=0..40mm width --
the sharp edge where the gusset's Toe->At face meets its Far2->Toe
(hypotenuse) face. This is geometrically the same toe configuration
that succeeded at R1mm in Stage 4b (the toe end of the gusset was
NOT changed by the Option A rework -- only the near-corner end was).

Does NOT yet include: the Ø5mm hole. That is the next increment.
"""

import gmsh
import sys
import os

INPUT_BREP = "mesh/full_bracket_study/frame_plus_gusset_check.brep"
OUTPUT_BREP = "mesh/full_bracket_study/frame_gusset_toefillet_check.brep"

TOE_X = 38.0
TOE_Z = 2.0
FILLET_RADIUS = 1.0
TOL = 1e-4

def main():
    if not os.path.exists(INPUT_BREP):
        print(f"ERROR: input file not found: {INPUT_BREP}")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("frame_gusset_toefillet")

    gmsh.model.occ.importShapes(INPUT_BREP)
    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=3)
    print(f"Import: {len(volumes)} volume(s)")
    if len(volumes) != 1:
        print("ERROR: expected exactly 1 volume on import.")
        sys.exit(1)
    vol_tag = volumes[0][1]
    pre_mass = gmsh.model.occ.getMass(3, vol_tag)
    print(f"Pre-fillet mass: {pre_mass:.4f} mm^3")

    # --- Locate the toe edge by bounding box ---
    # Expect: a straight line at x=38, z=2, spanning y=0 to y=40.
    print(f"\nSearching for toe edge at x={TOE_X}, z={TOE_Z}, y=0..40 ...")
    curves = gmsh.model.getEntities(dim=1)
    candidates = []
    for (dim, tag) in curves:
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bb
        if (abs(xmin - TOE_X) < TOL and abs(xmax - TOE_X) < TOL and
            abs(zmin - TOE_Z) < TOL and abs(zmax - TOE_Z) < TOL and
            abs(ymin - 0.0) < TOL and abs(ymax - 40.0) < TOL):
            candidates.append(tag)
            print(f"  candidate curve {tag}: bbox x=[{xmin:.4f},{xmax:.4f}] "
                  f"y=[{ymin:.4f},{ymax:.4f}] z=[{zmin:.4f},{zmax:.4f}]")

    if len(candidates) == 0:
        print("ERROR: no curve found matching the expected toe edge location.")
        print("Full curve dump for manual inspection:")
        for (dim, tag) in curves:
            bb = gmsh.model.occ.getBoundingBox(dim, tag)
            print(f"  curve {tag}: bbox = {[round(v,4) for v in bb]}")
        gmsh.finalize()
        sys.exit(1)
    if len(candidates) > 1:
        print(f"ERROR: {len(candidates)} curves matched -- ambiguous, need manual inspection.")
        gmsh.finalize()
        sys.exit(1)

    toe_edge_tag = candidates[0]
    print(f"\nToe edge identified: curve {toe_edge_tag}")

    # --- Apply fillet ---
    print(f"\nApplying R{FILLET_RADIUS}mm fillet to curve {toe_edge_tag} ...")
    try:
        filleted = gmsh.model.occ.fillet(
            [vol_tag], [toe_edge_tag], [FILLET_RADIUS], removeVolume=True
        )
        gmsh.model.occ.synchronize()
    except Exception as e:
        print(f"ERROR: fillet operation raised an exception: {e}")
        gmsh.finalize()
        sys.exit(1)

    fillet_volumes = [e for e in filleted if e[0] == 3]
    print(f"Fillet result: {len(fillet_volumes)} volume(s)")
    if len(fillet_volumes) != 1:
        print("ERROR: fillet did not produce a single solid volume.")
        gmsh.finalize()
        sys.exit(1)

    result_vol_tag = fillet_volumes[0][1]

    # --- Diagnostics ---
    post_mass = gmsh.model.occ.getMass(3, result_vol_tag)
    mass_removed = pre_mass - post_mass
    print(f"\nPost-fillet mass: {post_mass:.4f} mm^3")
    print(f"Mass removed by fillet: {mass_removed:.4f} mm^3")
    # Sanity expectation: quarter-circle fillet removes a prism of cross-section
    # (r^2 - pi*r^2/4) along the edge length. For r=1mm, edge length 40mm:
    expected_removed = (FILLET_RADIUS**2 - (3.14159265 * FILLET_RADIUS**2 / 4)) * 40.0
    print(f"Rough expected mass removed (corner-cut approx, r=1mm x 40mm length): "
          f"{expected_removed:.4f} mm^3")

    bbox = gmsh.model.occ.getBoundingBox(3, result_vol_tag)
    print(f"\nBounding box: xmin={bbox[0]:.4f} ymin={bbox[1]:.4f} zmin={bbox[2]:.4f} "
          f"xmax={bbox[3]:.4f} ymax={bbox[4]:.4f} zmax={bbox[5]:.4f}")

    surfaces = gmsh.model.getEntities(dim=2)
    print(f"\nTotal surfaces: {len(surfaces)}")
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

    # Identify the new fillet surface: should be a cylindrical surface,
    # bbox extent in x and z both ~1.0mm (quarter-circle, r=1mm), y-span 40mm.
    print("\nSearching for the new toe-fillet cylindrical surface...")
    found_fillet_surf = False
    for (dim, tag) in surfaces:
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        xext = bb[3] - bb[0]
        zext = bb[5] - bb[2]
        yext = bb[4] - bb[1]
        if abs(xext - FILLET_RADIUS) < 0.01 and abs(zext - FILLET_RADIUS) < 0.01 \
           and abs(yext - 40.0) < 0.01:
            print(f"  surface {tag}: x-extent={xext:.4f} z-extent={zext:.4f} "
                  f"y-extent={yext:.4f} -- matches expected R{FILLET_RADIUS}mm quarter-circle fillet")
            found_fillet_surf = True
    if not found_fillet_surf:
        print("  WARNING: no surface found matching expected fillet bbox signature.")
        print("  This does not necessarily mean the fillet failed -- verify manually.")

    os.makedirs(os.path.dirname(OUTPUT_BREP), exist_ok=True)
    gmsh.write(OUTPUT_BREP)
    print(f"\nSaved to {OUTPUT_BREP}")

    gmsh.finalize()

if __name__ == "__main__":
    main()
