"""
Stage 6 -- Parametric geometry builder.

Generalizes the exact, verified Stage 5 construction chain
(build_full_bracket_base_frame.py -> build_full_bracket_gusset.py ->
build_full_bracket_toe_fillet.py -> build_full_bracket_hole.py) over
three free design variables:

    t         flange thickness           (Stage 5 baseline: 4.0 mm)
    L_gusset  nominal gusset leg length  (Stage 5 baseline: 20.0 mm)
    r_toe     toe fillet radius          (Stage 5 baseline: 1.0 mm)

Fixed (per Stage 6 sign-off, not swept):
    R3 = 3.0 mm      (interior flange-junction fillet radius)
    HOLE_D = 5.0 mm, at (HOLE_X, HOLE_Y) = (20.0, 20.0)
    L = 60.0 mm (flange length), WIDTH = 40.0 mm (frame width)

TANGENT-POINT FORMULA (derived from and verified against the locked
Stage 5 base-frame script -- reproduces At=(55,2), Bt=(58,5),
Toe=(38,2), Far2=(58,22) exactly at t=4, L_gusset=20, R3=3):

    corner_x = L - t/2
    corner_z = t/2
    At   = (corner_x - R3, corner_z)         -- R3 tangent, Flange A side
    Bt   = (corner_x, corner_z + R3)         -- R3 tangent, Flange B side
    Toe  = (corner_x - L_gusset, corner_z)   -- measured from the ORIGINAL
    Far2 = (corner_x, corner_z + L_gusset)      sharp-corner point, per the
                                                 Stage 5 gusset script's logic

This script performs geometry construction and validity diagnostics
ONLY -- no meshing, solving, or optimization. It is meant to be run
manually across a handful of representative/corner-case design points
first, per the Stage 6 sign-off ("validate geometry before writing
mesh/solve automation"), before any LHS sweep or automation is built.

USAGE:
    python3 build_parametric_bracket.py <t> <L_gusset> <r_toe> [tag]

Exits non-zero (with a printed reason) on ANY geometry validity
failure -- degenerate fillet, non-single volume, sliver surfaces,
hole/gusset overlap, etc. Per the Stage 6 rule, a design point that
fails here is NOT silently skipped or fixed -- it is reported, and
must be resolved or excluded explicitly before being included in any
optimization claim.
"""

import gmsh
import sys
import os

TOL = 1e-6

# ---- Fixed parameters (not swept) ----
L = 60.0
WIDTH = 40.0
R3 = 3.0
HOLE_D = 5.0
HOLE_X = 20.0
HOLE_Y = 20.0
HOLE_R = HOLE_D / 2.0
CUTTER_Z0 = -5.0
CUTTER_LEN = 10.0

# ---- Validity gate thresholds ----
# Effective gusset leg (L_gusset - R3) must comfortably exceed the toe
# fillet radius, by the same order of margin Stage 3's r/h investigation
# showed matters. Conservative factor of 4x, not just >0.
MIN_LEG_TO_TOE_RATIO = 4.0
# Hole must clear the gusset's Toe x-position by at least this margin.
MIN_HOLE_CLEARANCE = 5.0
SLIVER_AREA_THRESHOLD = 1.0  # mm^2, same threshold used throughout Stage 5


def check_design_validity(t, L_gusset, r_toe):
    """Pre-flight checks BEFORE touching OCC. Returns (ok: bool, reasons: list[str])."""
    reasons = []

    corner_x = L - t / 2.0
    corner_z = t / 2.0
    toe_x = corner_x - L_gusset
    effective_leg = L_gusset - R3

    if effective_leg <= 0:
        reasons.append(
            f"Effective gusset leg (L_gusset - R3 = {L_gusset} - {R3} = "
            f"{effective_leg:.3f}) is non-positive. Geometrically invalid."
        )
    elif effective_leg / r_toe < MIN_LEG_TO_TOE_RATIO:
        reasons.append(
            f"Effective gusset leg ({effective_leg:.3f}mm) is less than "
            f"{MIN_LEG_TO_TOE_RATIO}x the toe fillet radius ({r_toe}mm) "
            f"(ratio={effective_leg/r_toe:.2f}). Flagged as a Stage-3-style "
            f"r/h degeneracy risk -- not attempted."
        )

    hole_edge_x = HOLE_X + HOLE_R
    clearance = toe_x - hole_edge_x
    if clearance < MIN_HOLE_CLEARANCE:
        reasons.append(
            f"Hole clearance to gusset toe ({clearance:.3f}mm) is below the "
            f"{MIN_HOLE_CLEARANCE}mm minimum margin (toe_x={toe_x:.3f}, "
            f"hole edge={hole_edge_x:.3f}). Hole may overlap or sit too "
            f"close to the gusset footprint."
        )

    return (len(reasons) == 0), reasons


def build_bracket_geometry(t, L_gusset, r_toe, tag, out_dir):
    """
    Builds the full Stage 6 parametric bracket geometry for one design
    point. Returns (success: bool, brep_path: str or None, mass: float or None,
    diagnostics: dict).
    """
    diagnostics = {"tag": tag, "t": t, "L_gusset": L_gusset, "r_toe": r_toe, "R3": R3}

    ok, reasons = check_design_validity(t, L_gusset, r_toe)
    diagnostics["prevalidity_ok"] = ok
    diagnostics["prevalidity_reasons"] = reasons
    if not ok:
        print(f"\n[{tag}] PRE-VALIDITY CHECK FAILED:")
        for r in reasons:
            print(f"  - {r}")
        return False, None, None, diagnostics

    half_t = t / 2.0
    corner_x = L - half_t
    corner_z = half_t

    at = (corner_x - R3, corner_z)
    bt = (corner_x, corner_z + R3)
    toe = (corner_x - L_gusset, corner_z)
    far2 = (corner_x, corner_z + L_gusset)

    diagnostics["tangent_At"] = at
    diagnostics["tangent_Bt"] = bt
    diagnostics["gusset_Toe"] = toe
    diagnostics["gusset_Far2"] = far2

    gmsh.initialize()
    gmsh.model.add(f"stage6_design_{tag}")
    occ = gmsh.model.occ

    print(f"\n{'='*70}\n[{tag}] t={t}, L_gusset={L_gusset}, r_toe={r_toe}, R3={R3}\n{'='*70}")

    # ---- Step 1: flanges, fuse ----
    flange_a = occ.addBox(0.0, 0.0, -half_t, L, WIDTH, t)
    flange_b = occ.addBox(L - half_t, 0.0, 0.0, t, WIDTH, L)
    occ.synchronize()
    out, _ = occ.fuse([(3, flange_a)], [(3, flange_b)])
    occ.synchronize()
    vol = out[0][1]
    print(f"[{tag}] Step 1 (flanges) OK. bbox={gmsh.model.getBoundingBox(3, vol)}")

    # ---- Step 2: identify reentrant corner edge, apply R3 fillet ----
    all_curves = gmsh.model.getEntities(1)
    reentrant_edge = None
    for dim, curve_tag in all_curves:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, curve_tag)
        x_at = abs(xmin - corner_x) < TOL and abs(xmax - corner_x) < TOL
        z_at = abs(zmin - corner_z) < TOL and abs(zmax - corner_z) < TOL
        spans_width = (ymax - ymin) > (WIDTH - 1.0)
        if x_at and z_at and spans_width:
            reentrant_edge = curve_tag
            break

    if reentrant_edge is None:
        print(f"[{tag}] FATAL: could not identify reentrant corner edge.")
        gmsh.finalize()
        diagnostics["failure_stage"] = "corner_edge_identification"
        return False, None, None, diagnostics

    try:
        filleted = occ.fillet([vol], [reentrant_edge], [R3], True)
        occ.synchronize()
        vol = filleted[0][1]
    except Exception as e:
        print(f"[{tag}] FATAL: R3 fillet failed: {e}")
        gmsh.finalize()
        diagnostics["failure_stage"] = "R3_fillet"
        diagnostics["failure_reason"] = str(e)
        return False, None, None, diagnostics
    print(f"[{tag}] Step 2 (R3={R3}mm interior fillet) OK.")

    # ---- Step 3: gusset (quadrilateral, tangent-point vertices) ----
    p1 = occ.addPoint(toe[0], 0.0, toe[1])
    p2 = occ.addPoint(at[0], 0.0, at[1])
    p3 = occ.addPoint(bt[0], 0.0, bt[1])
    p4 = occ.addPoint(far2[0], 0.0, far2[1])
    l1 = occ.addLine(p1, p2)
    l2 = occ.addLine(p2, p3)
    l3 = occ.addLine(p3, p4)
    l4 = occ.addLine(p4, p1)
    loop = occ.addCurveLoop([l1, l2, l3, l4])
    surf = occ.addPlaneSurface([loop])
    occ.synchronize()

    extruded = occ.extrude([(2, surf)], 0, WIDTH, 0)
    occ.synchronize()
    gusset_vol_tags = [e[1] for e in extruded if e[0] == 3]
    if len(gusset_vol_tags) != 1:
        print(f"[{tag}] FATAL: gusset extrusion produced {len(gusset_vol_tags)} volumes, expected 1.")
        gmsh.finalize()
        diagnostics["failure_stage"] = "gusset_extrusion"
        return False, None, None, diagnostics
    gusset_vol = gusset_vol_tags[0]

    fused, _ = occ.fuse([(3, vol)], [(3, gusset_vol)])
    occ.synchronize()
    if len(fused) != 1:
        print(f"[{tag}] FATAL: gusset fuse produced {len(fused)} volumes, expected 1.")
        gmsh.finalize()
        diagnostics["failure_stage"] = "gusset_fuse"
        return False, None, None, diagnostics
    vol = fused[0][1]
    print(f"[{tag}] Step 3 (gusset, effective leg={L_gusset - R3:.3f}mm) OK.")

    # ---- Step 4: toe fillet ----
    curves = gmsh.model.getEntities(1)
    toe_edge = None
    for dim, curve_tag in curves:
        bb = gmsh.model.getBoundingBox(dim, curve_tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bb
        if (abs(xmin - toe[0]) < 1e-4 and abs(xmax - toe[0]) < 1e-4 and
                abs(zmin - toe[1]) < 1e-4 and abs(zmax - toe[1]) < 1e-4 and
                abs(ymin - 0.0) < 1e-4 and abs(ymax - WIDTH) < 1e-4):
            toe_edge = curve_tag
            break

    if toe_edge is None:
        print(f"[{tag}] FATAL: could not identify toe edge at x={toe[0]}, z={toe[1]}.")
        gmsh.finalize()
        diagnostics["failure_stage"] = "toe_edge_identification"
        return False, None, None, diagnostics

    try:
        filleted = occ.fillet([vol], [toe_edge], [r_toe], True)
        occ.synchronize()
        fillet_vols = [e for e in filleted if e[0] == 3]
        if len(fillet_vols) != 1:
            raise RuntimeError(f"toe fillet produced {len(fillet_vols)} volumes, expected 1")
        vol = fillet_vols[0][1]
    except Exception as e:
        print(f"[{tag}] FATAL: toe fillet (r={r_toe}mm) failed: {e}")
        gmsh.finalize()
        diagnostics["failure_stage"] = "toe_fillet"
        diagnostics["failure_reason"] = str(e)
        return False, None, None, diagnostics
    print(f"[{tag}] Step 4 (toe fillet, r={r_toe}mm) OK.")

    # ---- Step 5: hole cut (fixed geometry) ----
    cutter = occ.addCylinder(HOLE_X, HOLE_Y, CUTTER_Z0, 0, 0, CUTTER_LEN, HOLE_R)
    occ.synchronize()
    cut_result, _ = occ.cut([(3, vol)], [(3, cutter)])
    occ.synchronize()
    cut_vols = [e for e in cut_result if e[0] == 3]
    if len(cut_vols) != 1:
        print(f"[{tag}] FATAL: hole cut produced {len(cut_vols)} volumes, expected 1.")
        gmsh.finalize()
        diagnostics["failure_stage"] = "hole_cut"
        return False, None, None, diagnostics
    vol = cut_vols[0][1]
    print(f"[{tag}] Step 5 (hole, fixed Ø{HOLE_D}mm @ ({HOLE_X},{HOLE_Y})) OK.")

    # ---- Post-build validity diagnostics ----
    mass = occ.getMass(3, vol)
    n_vols = len(gmsh.model.getEntities(3))
    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, s_tag), s_tag) for _, s_tag in surfaces]
    slivers = [(a, s_tag) for a, s_tag in areas if a < SLIVER_AREA_THRESHOLD]

    diagnostics["mass"] = mass
    diagnostics["n_volumes"] = n_vols
    diagnostics["n_surfaces"] = len(surfaces)
    diagnostics["sliver_surfaces"] = slivers
    diagnostics["postvalidity_ok"] = (n_vols == 1 and mass > 0 and len(slivers) == 0)

    print(f"[{tag}] Mass={mass:.4f} mm^3, volumes={n_vols}, surfaces={len(surfaces)}, "
          f"slivers={len(slivers)}")

    if not diagnostics["postvalidity_ok"]:
        print(f"[{tag}] POST-BUILD VALIDITY CHECK FAILED "
              f"(n_volumes={n_vols}, mass={mass}, slivers={slivers})")
        gmsh.finalize()
        return False, None, None, diagnostics

    os.makedirs(out_dir, exist_ok=True)
    brep_path = os.path.join(out_dir, f"stage6_design_{tag}.brep")
    gmsh.write(brep_path)
    print(f"[{tag}] PASSED all validity checks. Saved: {brep_path}")

    gmsh.finalize()
    return True, brep_path, mass, diagnostics


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 build_parametric_bracket.py <t> <L_gusset> <r_toe> [tag]")
        sys.exit(1)

    t_arg = float(sys.argv[1])
    L_gusset_arg = float(sys.argv[2])
    r_toe_arg = float(sys.argv[3])
    tag_arg = sys.argv[4] if len(sys.argv) > 4 else "manual"

    out_dir = "/mesh/stage6_parametric"

    success, brep_path, mass, diag = build_bracket_geometry(
        t_arg, L_gusset_arg, r_toe_arg, tag_arg, out_dir
    )

    print(f"\n{'='*70}\nRESULT: {'PASS' if success else 'FAIL'}\n{'='*70}")
    if success:
        print(f"Mass: {mass:.4f} mm^3")
        print(f"BRep: {brep_path}")
    else:
        print("Design point REJECTED. See diagnostics above.")
        sys.exit(1)
