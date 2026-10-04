"""
Stage 3 buckling pipeline check: curved skin panel geometry, root bay
upper surface, spar (x=0, 30% chord) to trailing edge, z=0 to 200mm
(root rib to first interior rib). NOT the full wing skin -- this is
the first curved-geometry test of the shell/*BUCKLE pipeline proven on
the flat plate.

Builds an OPEN surface (loft between two upper-surface curves), not a
closed volume -- this is a shell panel, matching the flat-plate
verification's use of shell elements.

Geometry construction and validation ONLY. Does not mesh.
"""
import gmsh
import math
import os
import sys

TOL = 1e-6
ROOT_CHORD, TIP_CHORD = 300.0, 150.0
SPAN = 1000.0
SPAR_CHORD_FRAC = 0.30
Z_STATIONS = [0.0, 200.0]  # root bay: root rib to first interior rib
N_PTS = 30  # cosine-spaced points, spar to trailing edge only

OUT_DIR = "mesh/stage3_panel_verify"
os.makedirs(OUT_DIR, exist_ok=True)


def chord_at(z):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * (z / SPAN)


def naca0012_half_thickness(x_over_c, t=0.12):
    return 5 * t * (0.2969 * math.sqrt(x_over_c) - 0.1260 * x_over_c
                     - 0.3516 * x_over_c**2 + 0.2843 * x_over_c**3
                     - 0.1015 * x_over_c**4)


def upper_surface_curve(chord):
    """Points from the spar (x=0, i.e. x/c=0.30) to the trailing edge
    (x=0.70*chord, x/c=1.0), on the upper surface."""
    pts = []
    for i in range(N_PTS):
        frac = i / (N_PTS - 1)
        x_over_c = SPAR_CHORD_FRAC + frac * (1.0 - SPAR_CHORD_FRAC)
        x = x_over_c * chord - SPAR_CHORD_FRAC * chord  # shift so spar is at x=0
        y = chord * naca0012_half_thickness(x_over_c)
        pts.append((x, y))
    return pts


def main():
    print("=" * 70)
    print("STAGE 3 PANEL VERIFY: root-bay upper-surface panel, z=0-200mm")
    print("=" * 70)

    gmsh.initialize()
    gmsh.model.add("panel_verify")
    occ = gmsh.model.occ

    wire_tags = []
    for z in Z_STATIONS:
        c = chord_at(z)
        pts = upper_surface_curve(c)
        pt_tags = [occ.addPoint(x, y, z) for x, y in pts]
        spline = occ.addSpline(pt_tags)
        wire = occ.addWire([spline])
        wire_tags.append(wire)
        print(f"z={z:.0f}mm: chord={c:.2f}mm, {len(pts)} points, "
              f"x range [{pts[0][0]:.2f}, {pts[-1][0]:.2f}], "
              f"y range [{pts[0][1]:.4f}, {pts[-1][1]:.4f}]")

    surf_tags = occ.addThruSections(wire_tags, makeSolid=False, makeRuled=False)
    occ.synchronize()
    surfaces = [t for d, t in surf_tags if d == 2]
    print(f"\naddThruSections returned {len(surf_tags)} entities, "
          f"{len(surfaces)} of dimension 2 (surface)")
    if len(surfaces) != 1:
        print(f"ERROR: expected 1 surface, got {len(surfaces)} -- STOP.")
        gmsh.finalize(); sys.exit(1)
    surf = surfaces[0]

    area = occ.getMass(2, surf)
    bbox = gmsh.model.getBoundingBox(2, surf)
    print(f"Surface area: {area:.4f} mm^2")
    print(f"Bbox: x=[{bbox[0]:.3f},{bbox[3]:.3f}] y=[{bbox[1]:.3f},{bbox[4]:.3f}] "
          f"z=[{bbox[2]:.3f},{bbox[5]:.3f}]")

    # sanity check: area should be roughly (avg chord-span-covered)*bay-length,
    # order-of-magnitude only, since it's a curved surface not flat
    avg_chord_covered = ((1 - SPAR_CHORD_FRAC) * chord_at(0) +
                          (1 - SPAR_CHORD_FRAC) * chord_at(200)) / 2.0
    flat_approx = avg_chord_covered * 200.0
    print(f"(Flat-projection area approx for sanity: {flat_approx:.1f} mm^2 "
          f"-- true curved area should be close to but slightly more than this)")

    n_vols = len(gmsh.model.getEntities(3))
    n_surfs = len(gmsh.model.getEntities(2))
    print(f"\nModel entities: {n_vols} volume(s) (expect 0, this is a shell), "
          f"{n_surfs} surface(s)")

    brep_path = os.path.join(OUT_DIR, "panel_verify.brep")
    step_path = os.path.join(OUT_DIR, "panel_verify.step")
    gmsh.write(brep_path)
    gmsh.write(step_path)
    print(f"\nSaved: {brep_path}\nSaved: {step_path}")

    gmsh.finalize()


if __name__ == "__main__":
    main()
