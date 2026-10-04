#!/usr/bin/env python3
"""
Project 03, Stage 3 -- skin imprint geometry.

Loads the validated spar v2 + ribs solid (mesh/stage3_solid/stage3_spar_ribs_v2.step),
rebuilds the wing envelope (same functions as build_stage3_solid.py), copies its 78
LATERAL faces (skin inner surface = airfoil envelope; end caps excluded), and
FRAGMENTS solid + skin faces so that:
  * skin pieces lying on the solid (rib rims, spar cap tops) become faces SHARED with
    the solid boundary ("shared" pieces -> shared nodes later),
  * the remaining skin pieces are free faces ("free" pieces -> shell only).

Independent checks (pure Python, no OCC):
  1. total skin area == lateral area of the envelope loft (planar quads)
  2. shared skin area == 6 rib-rim slabs + 2 x cap-top strips between ribs
  3. solid volume unchanged (128,582.1216 mm^3), still 1 volume

Writes mesh/stage3_skin/stage3_skin_geom.{step,brep} and stage3_skin_faces.json
(centroid, area, shared flag per skin piece) for the Netgen meshing stage.
Geometry ONLY. Run with system Python:
  /usr/bin/python3 scripts/build_stage3_skin_geom.py
"""
import json
import math
import os
import sys

import gmsh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_stage3_solid as b

STEP_IN = "mesh/stage3_solid/stage3_spar_ribs_v2.step"
OUT_DIR = "mesh/stage3_skin"
os.makedirs(OUT_DIR, exist_ok=True)
VOL_EXPECTED = 128582.1216
AREA_TOL_PCT = 0.1
SPAN = b.SPAN
RIB_T = b.RIB_T


def fail(msg):
    print(f"ERROR: {msg} -- STOP.")
    gmsh.finalize()
    sys.exit(1)


# ---------- independent area expectations (pure Python) ----------
def lateral_area(za, zb):
    """Area of the ruled lateral surface between z=za and z=zb: sum of planar quads."""
    pa = [(x, y, za) for x, y in b.outline_at(za)]
    pb = [(x, y, zb) for x, y in b.outline_at(zb)]
    n = len(pa)
    tot = 0.0
    for i in range(n):
        j = (i + 1) % n
        a, bb, c, d = pa[i], pa[j], pb[j], pb[i]
        d1 = [c[k] - a[k] for k in range(3)]
        d2 = [d[k] - bb[k] for k in range(3)]
        cx = (d1[1] * d2[2] - d1[2] * d2[1],
              d1[2] * d2[0] - d1[0] * d2[2],
              d1[0] * d2[1] - d1[1] * d2[0])
        tot += 0.5 * math.sqrt(cx[0] ** 2 + cx[1] ** 2 + cx[2] ** 2)
    return tot


def window_len(z, hw):
    xs, ys = b.upper_xy(z)
    L = 0.0
    for i in range(len(xs) - 1):
        x0, x1 = xs[i], xs[i + 1]
        lo, hi = max(min(x0, x1), -hw), min(max(x0, x1), hw)
        if hi <= lo:
            continue
        L += (hi - lo) / abs(x1 - x0) * math.hypot(x1 - x0, ys[i + 1] - ys[i])
    return L


def cap_strip_area(za, zb, hw, n=200):
    """One cap-top strip (|x|<=hw) between za and zb, incl. small slant factor."""
    def f(z):
        e = 1e-3
        dh = (b.h_env(min(z + e, SPAN)) - b.h_env(max(z - e, 0.0))) / (min(z + e, SPAN) - max(z - e, 0.0))
        return window_len(z, hw) * math.sqrt(1.0 + dh * dh)
    return b.simpson(f, za, zb, n)


def main():
    print("=" * 70)
    print("STAGE 3 SKIN IMPRINT GEOMETRY")
    print("=" * 70)
    if not os.path.isfile(STEP_IN):
        print(f"ERROR: {STEP_IN} not found."); sys.exit(1)

    gmsh.initialize()
    gmsh.model.add("stage3_skin_geom")
    occ = gmsh.model.occ

    tags = occ.importShapes(STEP_IN)
    solids = [t for d, t in tags if d == 3]
    if len(solids) != 1:
        fail(f"expected 1 solid in {STEP_IN}, got {len(solids)}")
    occ.synchronize()
    v0 = occ.getMass(3, solids[0])
    print(f"Imported solid volume: {v0:.4f} mm^3 (expect {VOL_EXPECTED})")
    if abs(v0 - VOL_EXPECTED) > 0.01:
        fail("imported solid volume mismatch")

    # envelope (same construction as build_stage3_solid.py)
    z_lo, z_hi = -RIB_T / 2.0, SPAN + RIB_T / 2.0
    env_lo = b.make_loop(occ, b.outline_at(z_lo), z_lo)
    env_hi = b.make_loop(occ, b.outline_at(z_hi), z_hi)
    env_tag = b.loft_volume(occ, env_lo, env_hi)
    bnd = gmsh.model.getBoundary([(3, env_tag)], oriented=False)
    lateral = []
    for d, t in bnd:
        bb = gmsh.model.getBoundingBox(2, t)
        if abs(bb[5] - bb[2]) < 1e-6:
            continue                     # end cap (planar, constant z)
        lateral.append(t)
    print(f"Envelope boundary faces: {len(bnd)}, lateral: {len(lateral)} "
          f"(polyline segments: {len(b.outline_at(0.0))})")
    if len(lateral) != len(b.outline_at(0.0)):
        fail("lateral face count != number of polyline segments")

    copies = occ.copy([(2, t) for t in lateral])
    occ.remove([(3, env_tag)], recursive=True)
    occ.synchronize()

    lat_area = sum(occ.getMass(2, t) for d, t in copies)
    lat_expected = lateral_area(z_lo, z_hi)
    pct = 100 * (lat_area - lat_expected) / lat_expected
    print(f"Skin lateral area (before imprint): {lat_area:.4f} mm^2, "
          f"independent {lat_expected:.4f} (%diff={pct:+.5f}%)")
    if abs(pct) > AREA_TOL_PCT:
        fail("skin area mismatch")

    # imprint
    out, outmap = occ.fragment([(3, solids[0])], copies)
    occ.synchronize()
    pieces = set()
    for m in outmap[1:]:
        for d, t in m:
            if d != 2:
                fail(f"unexpected dim {d} among skin pieces")
            pieces.add(t)

    vols = gmsh.model.getEntities(3)
    if len(vols) != 1:
        fail(f"expected 1 volume after fragment, got {len(vols)}")
    vol_tag = vols[0][1]
    v1 = occ.getMass(3, vol_tag)
    vpct = 100 * (v1 - VOL_EXPECTED) / VOL_EXPECTED
    print(f"\nAfter fragment: {len(vols)} volume, {v1:.4f} mm^3 ({vpct:+.5f}% vs validated)")
    if abs(vpct) > 0.001:
        fail("solid volume changed by fragment")

    solid_bnd = {t for d, t in gmsh.model.getBoundary([(3, vol_tag)], oriented=False)}
    shared = sorted(pieces & solid_bnd)
    free = sorted(pieces - solid_bnd)
    all_surfs = {t for d, t in gmsh.model.getEntities(2)}
    stray = all_surfs - solid_bnd - pieces
    print(f"Skin pieces: {len(pieces)}  (shared with solid: {len(shared)}, free: {len(free)})")
    print(f"Solid boundary faces: {len(solid_bnd)}; other (non-skin) solid faces: "
          f"{len(solid_bnd - pieces)}; stray surfaces not in solid/skin: {len(stray)}")
    if stray:
        fail(f"stray surfaces present: {sorted(stray)[:10]}")

    a_shared = sum(occ.getMass(2, t) for t in shared)
    a_free = sum(occ.getMass(2, t) for t in free)
    a_total = a_shared + a_free
    pct = 100 * (a_total - lat_area) / lat_area
    print(f"\nSkin area after imprint: shared {a_shared:.4f} + free {a_free:.4f} = "
          f"{a_total:.4f} mm^2 (%diff vs before {pct:+.5f}%)")
    if abs(pct) > 0.001:
        fail("skin area not conserved by fragment")

    # expected shared area: 6 rib-rim slabs + upper/lower cap-top strips between ribs
    rim = sum(lateral_area(zs - RIB_T / 2.0, zs + RIB_T / 2.0) for zs in b.RIB_STATIONS)
    caps = 0.0
    for k in range(len(b.RIB_STATIONS) - 1):
        za = b.RIB_STATIONS[k] + RIB_T / 2.0
        zb = b.RIB_STATIONS[k + 1] - RIB_T / 2.0
        caps += 2.0 * cap_strip_area(za, zb, b.CAP_W / 2.0)
    exp_shared = rim + caps
    spct = 100 * (a_shared - exp_shared) / exp_shared
    print(f"Expected shared area (independent): rim slabs {rim:.4f} + cap strips {caps:.4f} = "
          f"{exp_shared:.4f} mm^2 (%diff={spct:+.4f}%)")
    shared_ok = abs(spct) < AREA_TOL_PCT
    print(f"Shared-area check: [{'OK' if shared_ok else 'MISMATCH'}]")

    areas = sorted((occ.getMass(2, t), t, t in solid_bnd) for t in pieces)
    print("Smallest 5 skin pieces (area, tag, shared):", [(round(a, 3), t, s) for a, t, s in areas[:5]])

    faces_json = []
    for t in sorted(pieces):
        c = occ.getCenterOfMass(2, t)
        faces_json.append({"tag": t, "shared": t in solid_bnd, "area": occ.getMass(2, t),
                           "cx": c[0], "cy": c[1], "cz": c[2]})

    all_ok = shared_ok
    print("\n" + "=" * 70)
    print(f"OVERALL VALIDATION: {'PASS' if all_ok else 'FAIL'}")
    print("=" * 70)
    if all_ok:
        gmsh.write(os.path.join(OUT_DIR, "stage3_skin_geom.brep"))
        gmsh.write(os.path.join(OUT_DIR, "stage3_skin_geom.step"))
        with open(os.path.join(OUT_DIR, "stage3_skin_faces.json"), "w") as f:
            json.dump(faces_json, f, indent=1)
        print(f"Saved geometry and {len(faces_json)}-piece face list to {OUT_DIR}/")
    else:
        print("Geometry NOT saved -- validation failed.")
    gmsh.finalize()
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
