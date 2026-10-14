#!/usr/bin/env python3
"""
Project 03, Stage 3 geometry -- spar v2 (+ ribs), trimmed to the wing envelope.

WHY THIS EXISTS
  Stage 1/2 spar depth was 80% of local airfoil thickness, leaving a 3.6mm (root)
  / 1.8mm (tip) gap between spar caps and the airfoil outline, so a skin could
  not touch the spar. Option C (user decision): spar depth is resized so the cap
  outer faces sit ON the airfoil outline. This is a spar REDESIGN made to enable
  skin attachment; Stage 1/2 results describe the old 80%-depth spar.

METHOD
  1. Wing envelope: ruled loft of the root and tip NACA0012 polylines, from
     z=-RIB_T/2 to z=SPAN+RIB_T/2 (chord extrapolated linearly, so the root and
     tip ribs keep Stage 2's +/-0.75mm overhang). Because outline points scale
     linearly with chord about x=0, every envelope facet is PLANAR and every rib
     outline is an exact slice of the envelope.
  2. Oversized I-beam loft (root z=0 to tip z=SPAN): cap outer faces exceed the
     envelope by OVERSIZE; inner cap face at h_env(z) - CAP_T (h_env = envelope
     half-height at x=0). Web WEB_T, cap CAP_W x CAP_T as in Stage 1.
  3. Ribs mode: fuse one oversized slab (box) per rib station, one at a time.
  4. Boolean-intersect with the envelope: cap tops, rib rims, and (later) skin
     inner surface all lie on the same envelope facets.

  Ribs are therefore slices of the envelope: thickness 1.5mm, stations, and
  overhang as in Stage 2; rim follows the taper exactly (Stage 2 prisms did not).

The airfoil_outline() below is copied UNCHANGED from build_stage2_ribs.py
(including its non-symmetric trailing-edge closure -- flagged, not fixed here).

USAGE (system Python, not the conda env):
  /usr/bin/python3 scripts/build_stage3_solid.py spar    # spar v2 alone
  /usr/bin/python3 scripts/build_stage3_solid.py ribs    # spar v2 + 6 ribs

Geometry construction and validation ONLY. Does not mesh or solve.
"""
import bisect
import math
import os
import sys

import gmsh

TOL = 1e-3            # mm, bbox tolerance after booleans
SPAN = 1000.0
ROOT_CHORD, TIP_CHORD = 300.0, 150.0
SPAR_CHORD_FRAC = 0.30
import variant_tag as _vt
RIB_SPACING, RIB_STATIONS, RIB_T, RS_SUFFIX = _vt.configuration(SPAN)
N_AIRFOIL_PTS = 40

WEB_T = 2.0
CAP_W = 8.0
CAP_T = 3.0
OVERSIZE = 1.5        # mm, beam cap oversize beyond envelope (trimmed away)

RHO = 2780e-12        # t/mm^3
OLD_SPAR_MASS_G = 220.18

OUT_DIR = "mesh/stage3_solid"
os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# Airfoil helpers (copied unchanged from build_stage2_ribs.py)
# ----------------------------------------------------------------------------
def chord_at(z):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * (z / SPAN)


def naca0012_half_thickness(x_over_c, t=0.12):
    return 5 * t * (0.2969 * math.sqrt(x_over_c) - 0.1260 * x_over_c
                     - 0.3516 * x_over_c**2 + 0.2843 * x_over_c**3
                     - 0.1015 * x_over_c**4)


def airfoil_outline(chord, x_shift):
    xs = [chord * (1 - math.cos(math.pi * i / (2 * (N_AIRFOIL_PTS - 1))))
          for i in range(N_AIRFOIL_PTS)]
    upper = [(x - x_shift, chord * naca0012_half_thickness(x / chord)) for x in xs]
    lower = [(x - x_shift, -chord * naca0012_half_thickness(x / chord)) for x in reversed(xs[:-1])]
    return upper + lower[:-1]


def outline_at(z):
    c = chord_at(z)
    return airfoil_outline(c, SPAR_CHORD_FRAC * c)


# ----------------------------------------------------------------------------
# Numerical helpers (pure Python, for independent expected-volume checks)
# ----------------------------------------------------------------------------
def shoelace(pts):
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def upper_xy(z):
    up = outline_at(z)[:N_AIRFOIL_PTS]      # upper surface, LE -> TE, x ascending
    return [p[0] for p in up], [p[1] for p in up]


def interp_upper(xs, ys, x):
    i = bisect.bisect_right(xs, x) - 1
    i = max(0, min(i, len(xs) - 2))
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[i] + t * (ys[i + 1] - ys[i])


def h_env(z):
    """Envelope half-height at x=0 (spar axis)."""
    xs, ys = upper_xy(z)
    return interp_upper(xs, ys, 0.0)


def env_area(z):
    return shoelace(outline_at(z))


def spar_area(z, n=400):
    """Cross-section area of spar v2 at z (I-section trimmed to the envelope)."""
    xs, ys = upper_xy(z)
    h = interp_upper(xs, ys, 0.0)
    y_in = h - CAP_T
    web = WEB_T * 2.0 * y_in
    hw = CAP_W / 2.0
    dx = CAP_W / n
    vals = [interp_upper(xs, ys, -hw + i * dx) - y_in for i in range(n + 1)]
    integral = dx * (sum(vals) - 0.5 * (vals[0] + vals[-1]))
    return web + 2.0 * integral


def simpson(f, a, b, n):
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3.0


def ymax_in_band(z, half_w):
    xs, ys = upper_xy(z)
    cands = [interp_upper(xs, ys, -half_w), interp_upper(xs, ys, half_w)]
    cands += [y for x, y in zip(xs, ys) if -half_w < x < half_w]
    return max(cands)


# ----------------------------------------------------------------------------
# Gmsh helpers
# ----------------------------------------------------------------------------
def ibeam_profile(H, web_t, cap_w, cap_t):
    hw, cw, hf = web_t / 2.0, cap_w / 2.0, H / 2.0
    return [
        (-cw,  hf), ( cw,  hf), ( cw,  hf - cap_t), ( hw,  hf - cap_t),
        ( hw, -hf + cap_t), ( cw, -hf + cap_t), ( cw, -hf), (-cw, -hf),
        (-cw, -hf + cap_t), (-hw, -hf + cap_t), (-hw,  hf - cap_t), (-cw,  hf - cap_t),
    ]


def make_loop(occ, pts, z):
    pt_tags = [occ.addPoint(x, y, z) for x, y in pts]
    line_tags = [occ.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                 for i in range(len(pt_tags))]
    return occ.addCurveLoop(line_tags)


def loft_volume(occ, loop_a, loop_b):
    res = occ.addThruSections([loop_a, loop_b], makeSolid=True, makeRuled=True)
    occ.synchronize()
    return [t for d, t in res if d == 3][0]


def fail(msg):
    print(f"ERROR: {msg} -- STOP.")
    gmsh.finalize()
    sys.exit(1)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "ribs"
    if mode not in ("spar", "ribs"):
        print("usage: build_stage3_solid.py [spar|ribs]")
        sys.exit(2)

    print("=" * 70)
    print(f"STAGE 3 SOLID (Option C, spar depth = airfoil depth): mode = {mode}")
    print("=" * 70)

    z_lo, z_hi = -RIB_T / 2.0, SPAN + RIB_T / 2.0
    h0, h1 = h_env(0.0), h_env(SPAN)
    print(f"Envelope half-height at spar axis: root {h0:.4f} mm, tip {h1:.4f} mm")
    print(f"Spar v2 web height (between caps): root {2*(h0-CAP_T):.4f} mm, "
          f"tip {2*(h1-CAP_T):.4f} mm; overall depth root {2*h0:.4f}, tip {2*h1:.4f}")

    gmsh.initialize()
    gmsh.model.add("stage3_solid")
    occ = gmsh.model.occ

    # 1. envelope
    env_lo = make_loop(occ, outline_at(z_lo), z_lo)
    env_hi = make_loop(occ, outline_at(z_hi), z_hi)
    env_tag = loft_volume(occ, env_lo, env_hi)
    env_vol = occ.getMass(3, env_tag)
    env_expected = simpson(env_area, z_lo, z_hi, 2)   # A ~ c^2 -> exact for quadratic
    pct = 100 * (env_vol - env_expected) / env_expected
    print(f"\nEnvelope volume: {env_vol:.4f} (expected {env_expected:.4f}, %diff={pct:+.5f}%)")
    if abs(pct) > 0.05:
        fail("envelope volume mismatch")

    # 2. oversized beam
    def beam_pts(z):
        h = h_env(z)
        return ibeam_profile(2.0 * (h + OVERSIZE), WEB_T, CAP_W, CAP_T + OVERSIZE)
    root_pts, tip_pts = beam_pts(0.0), beam_pts(SPAN)
    beam_tag = loft_volume(occ, make_loop(occ, root_pts, 0.0),
                           make_loop(occ, tip_pts, SPAN))
    beam_vol = occ.getMass(3, beam_tag)
    beam_expected = SPAN * (shoelace(root_pts) + shoelace(tip_pts)) / 2.0
    pct = 100 * (beam_vol - beam_expected) / beam_expected
    print(f"Oversized beam volume: {beam_vol:.4f} (expected {beam_expected:.4f}, %diff={pct:+.5f}%)")
    if abs(pct) > 0.05:
        fail("oversized beam volume mismatch")

    # 3. slabs (ribs mode), fused one at a time
    current = beam_tag
    if mode == "ribs":
        for zs in RIB_STATIONS:
            slab = occ.addBox(-100.0, -25.0, zs - RIB_T / 2.0, 330.0, 50.0, RIB_T)
            occ.synchronize()
            fused, _ = occ.fuse([(3, current)], [(3, slab)])
            occ.synchronize()
            vols = [t for d, t in fused if d == 3]
            if len(vols) != 1:
                fail(f"fuse at z={zs:.0f} gave {len(vols)} volumes")
            current = vols[0]
            print(f"  fused slab at z={zs:.0f} -> 1 volume")

    # 4. trim to envelope
    out, _ = occ.intersect([(3, current)], [(3, env_tag)])
    occ.synchronize()
    out_vols = [t for d, t in out if d == 3]
    all_vols = gmsh.model.getEntities(3)
    print(f"\nAfter trim: {len(out_vols)} result volume(s), {len(all_vols)} volume(s) in model")
    if len(out_vols) != 1 or len(all_vols) != 1:
        fail("trim did not leave exactly 1 volume")
    vol_tag = out_vols[0]
    actual_vol = occ.getMass(3, vol_tag)

    # expected volume (independent numerical integration)
    if mode == "spar":
        expected_vol = simpson(spar_area, 0.0, SPAN, 400)
    else:
        expected_vol = 0.0
        for i, zs in enumerate(RIB_STATIONS):
            a, b = zs - RIB_T / 2.0, zs + RIB_T / 2.0
            expected_vol += simpson(env_area, a, b, 2)
            if i + 1 < len(RIB_STATIONS):
                a2, b2 = zs + RIB_T / 2.0, RIB_STATIONS[i + 1] - RIB_T / 2.0
                expected_vol += simpson(spar_area, a2, b2, 100)
    vpct = 100 * (actual_vol - expected_vol) / expected_vol
    print(f"Volume: actual={actual_vol:.4f} mm^3, expected(numerical)={expected_vol:.4f} "
          f"(%diff={vpct:+.5f}%)")
    vol_ok = abs(vpct) < 0.05

    spar_vol_num = simpson(spar_area, 0.0, SPAN, 400)
    print(f"Spar v2 alone (numerical): {spar_vol_num:.2f} mm^3 = "
          f"{spar_vol_num*RHO*1e6:.2f} g  (old spar: {OLD_SPAR_MASS_G} g)")

    # bbox
    # OCC getBoundingBox() is LOOSE for trimmed faces (it returned the untrimmed
    # cutter/beam extents), so the tight bbox comes from the nodes of a coarse 2D
    # surface mesh instead: nodes lie on the true trimmed faces, and all extremes
    # of this planar-faceted geometry are at vertices.
    loose = gmsh.model.getBoundingBox(3, vol_tag)
    print(f"(OCC loose bbox, for the record: x=[{loose[0]:.3f},{loose[3]:.3f}] "
          f"y=[{loose[1]:.3f},{loose[4]:.3f}] z=[{loose[2]:.3f},{loose[5]:.3f}])")
    gmsh.option.setNumber("Mesh.MeshSizeMax", 10.0)
    gmsh.model.mesh.generate(2)
    _, ncoord, _ = gmsh.model.mesh.getNodes()
    nx, ny, nz = ncoord[0::3], ncoord[1::3], ncoord[2::3]
    bbox = (min(nx), min(ny), min(nz), max(nx), max(ny), max(nz))
    print(f"Bbox (from {len(nx)} surface-mesh nodes): x=[{bbox[0]:.4f},{bbox[3]:.4f}] "
          f"y=[{bbox[1]:.4f},{bbox[4]:.4f}] z=[{bbox[2]:.4f},{bbox[5]:.4f}]")

    # independent centroid check: OCC z-centroid vs numerical integral of z*A(z)
    if mode == "spar":
        segs = [(0.0, SPAN, spar_area, 400)]
    else:
        segs = []
        for i, zs in enumerate(RIB_STATIONS):
            segs.append((zs - RIB_T / 2.0, zs + RIB_T / 2.0, env_area, 2))
            if i + 1 < len(RIB_STATIONS):
                segs.append((zs + RIB_T / 2.0, RIB_STATIONS[i + 1] - RIB_T / 2.0,
                             spar_area, 100))
    zc_expected = sum(simpson(lambda z, f=f: z * f(z), a, b, n)
                      for a, b, f, n in segs) / expected_vol
    cx, cy, cz = occ.getCenterOfMass(3, vol_tag)
    print(f"Centroid z: OCC={cz:.4f} mm, numerical={zc_expected:.4f} mm "
          f"(diff={cz - zc_expected:+.5f} mm); centroid x={cx:.4f}, y={cy:.5f}")
    cent_ok = abs(cz - zc_expected) < 5e-3
    print(f"Centroid check: [{'OK' if cent_ok else 'MISMATCH'}]")
    vol_ok = vol_ok and cent_ok
    if mode == "spar":
        ym = ymax_in_band(0.0, CAP_W / 2.0)
        exp = (-CAP_W / 2.0, -ym, 0.0, CAP_W / 2.0, ym, SPAN)
    else:
        ol = outline_at(z_lo)
        exp = (min(p[0] for p in ol), min(p[1] for p in ol), z_lo,
               max(p[0] for p in ol), max(p[1] for p in ol), z_hi)
    print(f"Expected: x=[{exp[0]:.4f},{exp[3]:.4f}] y=[{exp[1]:.4f},{exp[4]:.4f}] "
          f"z=[{exp[2]:.4f},{exp[5]:.4f}]")
    bbox_ok = all(abs(a - b) < TOL for a, b in zip(bbox, exp))
    print(f"Bbox check: [{'OK' if bbox_ok else 'MISMATCH'}]")

    # cleanliness
    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, t), t) for d, t in surfaces]
    slivers = [(a, t) for a, t in areas if a < 0.5]
    boundary = gmsh.model.getBoundary([(3, vol_tag)], oriented=False)
    stray_ok = len(surfaces) == len(boundary)
    sliver_ok = len(slivers) == 0
    print(f"Surfaces: {len(surfaces)}, on solid boundary: {len(boundary)}, "
          f"slivers (<0.5mm^2): {len(slivers)}")
    if slivers:
        print("  smallest faces:", sorted(slivers)[:5])
    print(f"Cleanliness: [{'OK' if stray_ok else 'MISMATCH'}]")

    all_ok = vol_ok and bbox_ok and sliver_ok and stray_ok
    print("\n" + "=" * 70)
    print(f"OVERALL VALIDATION: {'PASS' if all_ok else 'FAIL'}")
    print("=" * 70)

    if all_ok:
        name = "stage3_spar_v2" if mode == "spar" else "stage3_spar_ribs_v2" + RS_SUFFIX
        brep_path = os.path.join(OUT_DIR, name + ".brep")
        step_path = os.path.join(OUT_DIR, name + ".step")
        gmsh.write(brep_path)
        gmsh.write(step_path)
        print(f"\nSaved: {brep_path}\nSaved: {step_path}")
        print(f"Mass at Al 2024-T3 density: {actual_vol*RHO*1e6:.2f} g")
    else:
        print("\nGeometry NOT saved -- validation failed.")

    gmsh.finalize()
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
