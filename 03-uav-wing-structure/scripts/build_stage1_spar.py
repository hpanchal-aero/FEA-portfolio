"""
Project 03, Stage 1 -- isolated tapered I-beam spar geometry.

Cross-section: standard I-beam (web + top/bottom caps), strong axis
aligned with the vertical/thickness direction (y), since aero lift and
self-weight both act perpendicular to the chord plane -- flanges resist
this bending, matching real wing-spar orientation.

Taper: root (z=0) height 28.8mm (80% of NACA0012 local thickness at
30% chord, 300mm root chord); tip (z=1000mm) height 14.4mm (80% of
local thickness at 300% chord... i.e. 150mm tip chord). Web thickness
and cap dimensions held constant along the span (only overall height
tapers). Built via addThruSections (loft) between the root and tip
I-beam profile wires.

Geometry construction and validation ONLY. Does not mesh or solve.
"""

import gmsh
import os
import sys

TOL = 1e-6
SPAN = 1000.0          # mm, half-span (z-direction)

# Root / tip overall spar height (80% of local NACA0012 thickness at 30% chord)
H_ROOT = 28.8
H_TIP = 14.4

WEB_T = 2.0             # mm, web thickness (constant along span)
CAP_W = 8.0             # mm, cap (flange) width (constant along span)
CAP_T = 3.0             # mm, cap thickness (constant along span)

RHO = 2780e-12          # t/mm^3, Al 2024-T3 (2780 kg/m^3 x 1e-12; to be confirmed/cited)

OUT_DIR = "mesh/stage1_spar"
os.makedirs(OUT_DIR, exist_ok=True)


def ibeam_profile(H, web_t, cap_w, cap_t):
    """Return the 12 (x,y) vertices of a symmetric I-beam cross-section,
    centered at the origin, strong axis along y."""
    hw, cw, hf = web_t / 2.0, cap_w / 2.0, H / 2.0
    return [
        (-cw,  hf), ( cw,  hf), ( cw,  hf - cap_t), ( hw,  hf - cap_t),
        ( hw, -hf + cap_t), ( cw, -hf + cap_t), ( cw, -hf), (-cw, -hf),
        (-cw, -hf + cap_t), (-hw, -hf + cap_t), (-hw,  hf - cap_t), (-cw,  hf - cap_t),
    ]


def shoelace_area(pts):
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def ibeam_area_analytic(H, web_t, cap_w, cap_t):
    web_h = H - 2 * cap_t
    return web_h * web_t + 2 * cap_w * cap_t


def make_wire(occ, pts, z):
    pt_tags = [occ.addPoint(x, y, z) for x, y in pts]
    line_tags = [occ.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                 for i in range(len(pt_tags))]
    return occ.addCurveLoop(line_tags), pt_tags


def main():
    print("=" * 70)
    print("STAGE 1: TAPERED I-BEAM SPAR")
    print("=" * 70)

    root_pts = ibeam_profile(H_ROOT, WEB_T, CAP_W, CAP_T)
    tip_pts = ibeam_profile(H_TIP, WEB_T, CAP_W, CAP_T)

    root_area_analytic = ibeam_area_analytic(H_ROOT, WEB_T, CAP_W, CAP_T)
    tip_area_analytic = ibeam_area_analytic(H_TIP, WEB_T, CAP_W, CAP_T)
    root_area_shoelace = shoelace_area(root_pts)
    tip_area_shoelace = shoelace_area(tip_pts)

    print(f"Root: H={H_ROOT}mm, area analytic={root_area_analytic:.4f} mm^2, "
          f"shoelace={root_area_shoelace:.4f} mm^2")
    print(f"Tip:  H={H_TIP}mm, area analytic={tip_area_analytic:.4f} mm^2, "
          f"shoelace={tip_area_shoelace:.4f} mm^2")
    if (abs(root_area_analytic - root_area_shoelace) > 1e-6 or
            abs(tip_area_analytic - tip_area_shoelace) > 1e-6):
        print("FATAL: profile area mismatch between analytic and shoelace calc.")
        sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("stage1_spar")
    occ = gmsh.model.occ

    root_loop, _ = make_wire(occ, root_pts, 0.0)
    tip_loop, _ = make_wire(occ, tip_pts, SPAN)
    occ.synchronize()

    loft = occ.addThruSections([root_loop, tip_loop], makeSolid=True, makeRuled=True)
    occ.synchronize()

    vols = gmsh.model.getEntities(3)
    print(f"\nLoft produced {len(vols)} volume(s)")
    if len(vols) != 1:
        print("FATAL: expected 1 volume.")
        gmsh.finalize()
        sys.exit(1)
    vol_tag = vols[0][1]

    actual_vol = occ.getMass(3, vol_tag)
    # expected volume for a LINEAR taper between two similar-shaped polygons
    # scaled only in H (web/cap widths/thicknesses constant): area is linear
    # in H for this profile (web_h*web_t term linear in H, cap term constant),
    # so volume = SPAN * (area_root + area_tip) / 2 is the correct closed form
    # for this specific taper (not a general loft-volume approximation).
    expected_vol = SPAN * (root_area_analytic + tip_area_analytic) / 2.0
    vol_pct = 100 * (actual_vol - expected_vol) / expected_vol
    print(f"Volume: actual={actual_vol:.4f} mm^3, expected(linear-taper closed form)="
          f"{expected_vol:.4f} mm^3 (%diff={vol_pct:+.4f}%)")
    vol_ok = abs(vol_pct) < 0.5

    bbox = gmsh.model.getBoundingBox(3, vol_tag)
    print(f"Bbox: x=[{bbox[0]:.4f},{bbox[3]:.4f}] y=[{bbox[1]:.4f},{bbox[4]:.4f}] "
          f"z=[{bbox[2]:.4f},{bbox[5]:.4f}]")
    bbox_ok = (abs(bbox[2] - 0.0) < TOL and abs(bbox[5] - SPAN) < TOL and
               abs(bbox[1] + H_ROOT / 2) < TOL and abs(bbox[4] - H_ROOT / 2) < TOL)
    print(f"Bbox check: [{'OK' if bbox_ok else 'MISMATCH -- STOP'}]")

    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, t), t) for d, t in surfaces]
    slivers = [(a, t) for a, t in areas if a < 0.5]
    print(f"Surfaces: {len(surfaces)}, slivers (<0.5mm^2): {len(slivers)}")
    sliver_ok = len(slivers) == 0

    boundary = gmsh.model.getBoundary([(3, vol_tag)], oriented=False)
    stray_ok = (len(surfaces) == len(boundary))
    print(f"Cleanliness: {len(surfaces)} surfaces, {len(boundary)} on solid boundary "
          f"[{'OK' if stray_ok else 'MISMATCH -- STOP'}]")

    all_ok = vol_ok and bbox_ok and sliver_ok and stray_ok
    print("\n" + "=" * 70)
    print(f"OVERALL VALIDATION: {'PASS' if all_ok else 'FAIL'}")
    print("=" * 70)

    if all_ok:
        brep_path = os.path.join(OUT_DIR, "stage1_spar.brep")
        step_path = os.path.join(OUT_DIR, "stage1_spar.step")
        gmsh.write(brep_path)
        gmsh.write(step_path)
        mass_g = actual_vol * RHO * 1e6
        print(f"\nSaved: {brep_path}\nSaved: {step_path}")
        print(f"Mass at Al 2024-T3 density (rho={RHO*1e12:.0f} kg/m^3): {mass_g:.2f} g")
    else:
        print("\nGeometry NOT saved -- validation failed.")

    gmsh.finalize()
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
