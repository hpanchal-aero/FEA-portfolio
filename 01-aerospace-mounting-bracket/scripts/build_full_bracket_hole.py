"""
Stage 5 -- Full L-bracket integration.
Increment 4: add the Ø5mm through-hole to Flange A at x=20, y=20,
through the verified frame+gusset+toe-fillet solid.

Hole axis: through Flange A's thickness direction (z), same logic as
Stage 2's isolated hole study, new in-plane placement (x=20, y=20)
per the Stage 5 spec. Flange A spans z=[-2,2], so the hole cutter
must fully traverse that range with margin on both sides to guarantee
a clean boolean cut.

Does NOT add anything further -- this is the last geometry increment
before meshing begins.
"""

import gmsh
import sys
import os

INPUT_BREP = "mesh/full_bracket_study/frame_gusset_toefillet_check.brep"
OUTPUT_BREP = "mesh/full_bracket_study/full_bracket_geometry_final.brep"

HOLE_X = 20.0
HOLE_Y = 20.0
HOLE_DIA = 5.0
HOLE_R = HOLE_DIA / 2.0

# Flange A thickness is z=[-2,2]; cutter must extend past both faces
CUTTER_Z0 = -5.0
CUTTER_LEN = 10.0  # spans z=-5 to z=5, comfortably past z=[-2,2]

def main():
    if not os.path.exists(INPUT_BREP):
        print(f"ERROR: input file not found: {INPUT_BREP}")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("full_bracket_with_hole")

    gmsh.model.occ.importShapes(INPUT_BREP)
    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=3)
    print(f"Import: {len(volumes)} volume(s)")
    if len(volumes) != 1:
        print("ERROR: expected exactly 1 volume on import.")
        sys.exit(1)
    vol_tag = volumes[0][1]
    pre_mass = gmsh.model.occ.getMass(3, vol_tag)
    print(f"Pre-hole mass: {pre_mass:.4f} mm^3")

    # --- Sanity: confirm hole location is actually on Flange A material ---
    # Flange A footprint: x in [0,60] roughly, y in [0,40], z in [-2,2].
    # x=20, y=20 should be well within plain flange material, away from
    # the gusset (gusset occupies x=[38,58] region) and away from edges.
    print(f"\nHole center (x={HOLE_X}, y={HOLE_Y}) sanity check:")
    print(f"  Flange A nominal footprint: x=[0,60], y=[0,40], z=[-2,2]")
    print(f"  Gusset footprint: x=[38,58] approx -- hole at x=20 is clear of gusset")
    print(f"  Edge distances: x-edge(0)={HOLE_X}mm, y-edges={HOLE_Y}mm and {40-HOLE_Y}mm")
    print(f"  All edge distances > 10mm (Ø5mm hole, r=2.5mm) -- no edge-proximity concern")

    # --- Build cylindrical cutter ---
    cutter_tag = gmsh.model.occ.addCylinder(
        HOLE_X, HOLE_Y, CUTTER_Z0,   # base center
        0, 0, CUTTER_LEN,             # axis direction + length
        HOLE_R
    )
    gmsh.model.occ.synchronize()
    cutter_mass = gmsh.model.occ.getMass(3, cutter_tag)
    print(f"\nCutter cylinder mass (pre-cut, standalone): {cutter_mass:.4f} mm^3")
    expected_cutter_mass = 3.14159265 * HOLE_R**2 * CUTTER_LEN
    print(f"Expected (pi*r^2*L): {expected_cutter_mass:.4f} mm^3")

    # --- Boolean cut ---
    cut_result, _ = gmsh.model.occ.cut([(3, vol_tag)], [(3, cutter_tag)])
    gmsh.model.occ.synchronize()

    cut_volumes = [e for e in cut_result if e[0] == 3]
    print(f"\nCut result: {len(cut_volumes)} volume(s)")
    if len(cut_volumes) != 1:
        print("ERROR: cut did not produce a single solid volume.")
        print("This usually means the cutter did not fully traverse the material,")
        print("or intersected an unexpected feature (e.g. gusset edge).")
        gmsh.finalize()
        sys.exit(1)

    final_vol_tag = cut_volumes[0][1]

    # --- Diagnostics ---
    post_mass = gmsh.model.occ.getMass(3, final_vol_tag)
    mass_removed = pre_mass - post_mass
    expected_removed = 3.14159265 * HOLE_R**2 * 4.0  # actual material thickness = 4mm (z=-2 to 2)
    print(f"\nPost-hole mass: {post_mass:.4f} mm^3")
    print(f"Mass removed: {mass_removed:.4f} mm^3")
    print(f"Expected removed (pi*r^2*4mm flange thickness): {expected_removed:.4f} mm^3")

    bbox = gmsh.model.occ.getBoundingBox(3, final_vol_tag)
    print(f"\nBounding box: xmin={bbox[0]:.4f} ymin={bbox[1]:.4f} zmin={bbox[2]:.4f} "
          f"xmax={bbox[3]:.4f} ymax={bbox[4]:.4f} zmax={bbox[5]:.4f}")
    # Should be unchanged from pre-hole bbox -- a through-hole in the interior
    # of a face should not alter the overall bounding box.
    print("(Should match pre-hole bbox exactly -- interior hole doesn't change overall extent)")

    surfaces = gmsh.model.getEntities(dim=2)
    print(f"\nTotal surfaces: {len(surfaces)}")
    areas = []
    for (dim, tag) in surfaces:
        area = gmsh.model.occ.getMass(2, tag)
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        areas.append((area, tag, bb))
    areas.sort()
    for area, tag, bb in areas:
        print(f"  surf {tag}: area={area:.4f}  bbox=[{bb[0]:.4f},{bb[1]:.4f},{bb[2]:.4f}] "
              f"-> [{bb[3]:.4f},{bb[4]:.4f},{bb[5]:.4f}]")

    slivers = [t for a, t, bb in areas if a < 1.0]
    if slivers:
        print(f"\nWARNING: possible sliver surfaces (<1.0 mm^2): {slivers}")
    else:
        print("\nNo sliver surfaces detected (<1.0 mm^2 threshold).")

    # --- Confirm the hole bore surface exists with correct geometry ---
    # Expected: a cylindrical surface, radius 2.5mm, spanning z=-2 to z=2
    # (the actual flange thickness), centered at x=20,y=20.
    print("\nSearching for hole bore surface...")
    bore_found = False
    expected_bore_area = 2 * 3.14159265 * HOLE_R * 4.0  # circumference x thickness
    for area, tag, bb in areas:
        zext = bb[5] - bb[2]
        xctr = (bb[0] + bb[3]) / 2
        yctr = (bb[1] + bb[4]) / 2
        if abs(zext - 4.0) < 0.05 and abs(xctr - HOLE_X) < 0.1 and abs(yctr - HOLE_Y) < 0.1:
            print(f"  surf {tag}: area={area:.4f} (expected ~{expected_bore_area:.4f}), "
                  f"z-span={bb[2]:.4f} to {bb[5]:.4f}, centered at ({xctr:.4f},{yctr:.4f})")
            print(f"  -- matches expected hole bore signature")
            bore_found = True
    if not bore_found:
        print("  WARNING: no surface matched expected bore signature -- verify manually.")

    os.makedirs(os.path.dirname(OUTPUT_BREP), exist_ok=True)
    gmsh.write(OUTPUT_BREP)
    print(f"\nSaved to {OUTPUT_BREP}")
    print("\nThis is the FINAL Stage 5 geometry (fillet + gusset + toe fillet + hole).")
    print("Next step after this is verified: meshing with local refinement fields.")

    gmsh.finalize()

if __name__ == "__main__":
    main()
