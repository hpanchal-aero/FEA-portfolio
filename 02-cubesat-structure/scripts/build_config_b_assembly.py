"""
Project 02, Configuration B -- pocketed-panel assembly geometry.

Same rail/panel layout as Configuration A, but each of the 4 panels gets a
through-cut pocket (one of 4 shapes) at the given removal level, centered at
z=50, with the 5mm perimeter frame retained -- matching the isolated-panel
design in build_config_b_panels.py (dimension formulas reused from there).

Usage: /usr/bin/python3 scripts/build_config_b_assembly.py <shape> <level_pct>
  shape: rect | circle | cross | grid
  level_pct: one of 10, 15, 20, 25, 30

Geometry construction and validation ONLY. Does not mesh or solve.
"""

import gmsh
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_config_b_panels import rect_dims, circle_dims, cross_dims, grid_dims  # noqa: E402

TOL = 1e-6

ENV = 50.0
RAIL = 8.5
HEIGHT = 100.0
PANEL_T = 2.0
RHO = 2700e-9
FRAME_W = 5.0

RAIL_INNER = ENV - RAIL          # 41.5
PANEL_LEN = 2 * RAIL_INNER       # 83.0  (matches PANEL_W in build_config_b_panels)
PANEL_AREA = PANEL_LEN * HEIGHT  # 8300  (matches PANEL_AREA there)
W_INT = PANEL_LEN - 2 * FRAME_W  # 73.0
H_INT = HEIGHT - 2 * FRAME_W     # 90.0

EXPECTED_RAIL_VOL_EACH = RAIL * RAIL * HEIGHT
EXPECTED_PANEL_VOL_EACH = PANEL_LEN * PANEL_T * HEIGHT
EXPECTED_XSEC_AREA = 4 * (RAIL * RAIL) + 4 * (PANEL_LEN * PANEL_T)  # 953, at z away from any pocket
SECTION_Z = 2.0    # inside the retained 5mm frame margin by construction for every shape
                    # (frame width FRAME_W=5mm is never cut into, so z<FRAME_W and z>HEIGHT-FRAME_W
                    # are always clear -- unlike z=10, which the cross shape's full-length vertical
                    # arm intersects at every removal level; see build_assembly_cross_10pct.log)

OUT_DIR = "mesh/config_b_study"
os.makedirs(OUT_DIR, exist_ok=True)


def cutter_local(occ, shape, A_target, big):
    """Build pocket cutter(s) in LOCAL (u=width, v=thru-thickness, w=height)
    coordinates, centered at the origin. Returns list of 3D tags."""
    if shape == "rect":
        w, h = rect_dims(A_target)
        if w is None:
            return None
        return [occ.addBox(-w / 2, -big / 2, -h / 2, w, big, h)]
    if shape == "circle":
        r = circle_dims(A_target)
        return [occ.addCylinder(0, -big / 2, 0, 0, big, 0, r)]
    if shape == "cross":
        w_arm = cross_dims(A_target)
        if w_arm is None:
            return None
        a1 = occ.addBox(-W_INT / 2, -big / 2, -w_arm / 2, W_INT, big, w_arm)
        a2 = occ.addBox(-w_arm / 2, -big / 2, -H_INT / 2, w_arm, big, H_INT)
        return [a1, a2]
    if shape == "grid":
        RIB = 8.0
        side = (A_target / 4) ** 0.5
        footprint = 2 * side + RIB
        if footprint > min(W_INT, H_INT):
            return None
        offset = side / 2 + RIB / 2
        tags = []
        for su in (-1, 1):
            for sw in (-1, 1):
                cu, cw = su * offset, sw * offset
                tags.append(occ.addBox(cu - side / 2, -big / 2, cw - side / 2, side, big, side))
        return tags
    raise ValueError(shape)


def place_cutters_global(occ, local_tags, orientation, center_uw, y_or_x_center):
    """Move/rotate a local (u,v,w)-centered cutter set into global coords for
    one of the 4 panel orientations, then return their (dim,tag) pairs."""
    out = []
    for t in local_tags:
        if orientation == "y_normal":
            # local u->x, v->y, w->z; just translate to (0, y_or_x_center, 50)
            occ.translate([(3, t)], center_uw[0], y_or_x_center, center_uw[1])
        else:  # "x_normal": local u->y, v->x, w->z -> rotate u-axis(x) into y-axis(y), then translate
            occ.rotate([(3, t)], 0, 0, 0, 0, 0, 1, 1.5707963267948966)  # +90 deg about z: x->y, y->-x
            occ.translate([(3, t)], y_or_x_center, center_uw[0], center_uw[1])
        out.append((3, t))
    return out


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("rect", "circle", "cross", "grid"):
        print("Usage: build_config_b_assembly.py <rect|circle|cross|grid> <level_pct>")
        sys.exit(1)
    shape = sys.argv[1]
    level = float(sys.argv[2])
    A_target = level / 100.0 * PANEL_AREA
    tag_name = f"config_b_assembly_{shape}_{int(level)}pct"
    removed_vol_each = A_target * PANEL_T
    expected_total_vol = (4 * EXPECTED_RAIL_VOL_EACH
                           + 4 * (EXPECTED_PANEL_VOL_EACH - removed_vol_each))

    gmsh.initialize()
    gmsh.model.add(tag_name)
    occ = gmsh.model.occ

    print("=" * 70)
    print(f"CONFIG B ASSEMBLY: shape={shape}, level={level}%, "
          f"target removed area/panel={A_target:.4f} mm^2")
    print("=" * 70)

    print("\nSTEP 1: rails (unchanged from Config A)")
    rail_positions = [(-ENV, -ENV), (RAIL_INNER, -ENV), (-ENV, RAIL_INNER), (RAIL_INNER, RAIL_INNER)]
    rail_tags = []
    for i, (x0, y0) in enumerate(rail_positions, start=1):
        t = occ.addBox(x0, y0, 0.0, RAIL, RAIL, HEIGHT)
        occ.synchronize()
        bbox = gmsh.model.getBoundingBox(3, t)
        ok = (abs(bbox[3] - bbox[0] - RAIL) < TOL and abs(bbox[4] - bbox[1] - RAIL) < TOL
              and abs(bbox[5] - bbox[2] - HEIGHT) < TOL)
        print(f"  Rail {i}: [{'OK' if ok else 'MISMATCH -- STOP'}]")
        if not ok:
            gmsh.finalize(); sys.exit(1)
        rail_tags.append(t)

    print("\nSTEP 2: panels with pocket")
    panel_defs = [
        ("y=-ENV", -RAIL_INNER, -ENV, PANEL_LEN, PANEL_T, "y_normal", -ENV + PANEL_T / 2),
        ("y=+ENV", -RAIL_INNER, ENV - PANEL_T, PANEL_LEN, PANEL_T, "y_normal", ENV - PANEL_T / 2),
        ("x=-ENV", -ENV, -RAIL_INNER, PANEL_T, PANEL_LEN, "x_normal", -ENV + PANEL_T / 2),
        ("x=+ENV", ENV - PANEL_T, -RAIL_INNER, PANEL_T, PANEL_LEN, "x_normal", ENV - PANEL_T / 2),
    ]
    panel_tags = []
    for label, x0, y0, dx, dy, orient, thru_center in panel_defs:
        panel = occ.addBox(x0, y0, 0.0, dx, dy, HEIGHT)
        occ.synchronize()
        bbox = gmsh.model.getBoundingBox(3, panel)
        ok = (abs(bbox[3] - bbox[0] - dx) < TOL and abs(bbox[4] - bbox[1] - dy) < TOL
              and abs(bbox[5] - bbox[2] - HEIGHT) < TOL)
        print(f"  Panel [{label}] pre-cut: [{'OK' if ok else 'MISMATCH -- STOP'}]")
        if not ok:
            gmsh.finalize(); sys.exit(1)

        local = cutter_local(occ, shape, A_target, PANEL_T * 4)
        if local is None:
            print(f"  [{tag_name}] INFEASIBLE: pocket does not fit at this level.")
            gmsh.finalize(); sys.exit(1)
        occ.synchronize()
        # Verify no cutter extends into the retained frame margin (local w-axis maps to
        # global z before translation to z-center 50). Bbox is in LOCAL coords here.
        for lt in local:
            lbb = gmsh.model.getBoundingBox(3, lt)
            w_half_extent = max(abs(lbb[2]), abs(lbb[5]))  # local w is bbox index 2/5 (z of local box)
            if w_half_extent > (HEIGHT / 2 - FRAME_W) + TOL:
                print(f"  [{tag_name}] FATAL: cutter extends into the frame margin "
                      f"(half-extent {w_half_extent:.4f} > {HEIGHT/2 - FRAME_W:.4f}). "
                      f"SECTION_Z={SECTION_Z} would be inside the pocket for this shape.")
                gmsh.finalize(); sys.exit(1)
        cutters = place_cutters_global(occ, local, orient, (0.0, 50.0), thru_center)
        occ.synchronize()

        result, _ = occ.cut([(3, panel)], cutters)
        occ.synchronize()
        if len(result) != 1:
            print(f"  [{tag_name}] FATAL: pocket cut on panel [{label}] gave {len(result)} volumes.")
            gmsh.finalize(); sys.exit(1)
        cut_panel = result[0][1]
        post_vol = occ.getMass(3, cut_panel)
        removed_vol = EXPECTED_PANEL_VOL_EACH - post_vol
        pct = 100 * (removed_vol - removed_vol_each) / removed_vol_each
        print(f"    post-cut volume={post_vol:.4f} mm^3, removed={removed_vol:.4f} "
              f"(target {removed_vol_each:.4f}, %diff={pct:+.4f}%)")
        if abs(pct) > 1.0:
            print("  MISMATCH -- STOP")
            gmsh.finalize(); sys.exit(1)
        panel_tags.append(cut_panel)

    print("\nSTEP 3: rail/panel overlap check (pre-fuse)")
    rail_copies = [occ.copy([(3, t)])[0][1] for t in rail_tags]
    panel_copies = [occ.copy([(3, t)])[0][1] for t in panel_tags]
    try:
        inter, _ = occ.intersect([(3, t) for t in rail_copies], [(3, t) for t in panel_copies],
                                  removeObject=True, removeTool=True)
        occ.synchronize()
        overlap = sum(occ.getMass(3, t) for d, t in inter) if inter else 0.0
        if inter:
            occ.remove(inter, recursive=True); occ.synchronize()
    except Exception:
        overlap = 0.0
    print(f"  Overlap volume: {overlap:.6f} mm^3 [{'OK' if overlap < 1e-6 else 'WARNING'}]")

    print("\nSTEP 4: fuse")
    all_tags = rail_tags + panel_tags
    fused, _ = occ.fuse([(3, all_tags[0])], [(3, t) for t in all_tags[1:]])
    occ.synchronize()
    if len(fused) != 1:
        print(f"FATAL: fuse gave {len(fused)} volumes, expected 1.")
        gmsh.finalize(); sys.exit(1)
    vol_tag = fused[0][1]

    print("\nSTEP 5: validation")
    actual_vol = occ.getMass(3, vol_tag)
    vol_pct = 100 * (actual_vol - expected_total_vol) / expected_total_vol
    print(f"  Volume: actual={actual_vol:.4f}, expected={expected_total_vol:.4f} "
          f"(%diff={vol_pct:+.4f}%) [{'OK' if abs(vol_pct) < 0.1 else 'MISMATCH -- STOP'}]")
    if abs(vol_pct) >= 0.1:
        gmsh.finalize(); sys.exit(1)

    bbox = gmsh.model.getBoundingBox(3, vol_tag)
    bbox_ok = (abs(bbox[0] + ENV) < TOL and abs(bbox[3] - ENV) < TOL and
               abs(bbox[1] + ENV) < TOL and abs(bbox[4] - ENV) < TOL and
               abs(bbox[2]) < TOL and abs(bbox[5] - HEIGHT) < TOL)
    print(f"  Bbox: {'OK' if bbox_ok else 'MISMATCH -- STOP'}")
    if not bbox_ok:
        gmsh.finalize(); sys.exit(1)

    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, t), t) for d, t in surfaces]
    slivers = [(a, t) for a, t in areas if a < 1.0]
    print(f"  Surfaces: {len(surfaces)}, slivers: {len(slivers)}")
    sliver_ok = len(slivers) == 0

    print(f"\n  Sectional check at z={SECTION_Z} (outside pocket for all levels, expected {EXPECTED_XSEC_AREA}):")
    section_result = []
    try:
        plane = occ.addDisk(0, 0, SECTION_Z, 200, 200)
        occ.synchronize()
        solid_copy = occ.copy([(3, vol_tag)])[0][1]
        section_result, _ = occ.intersect([(2, plane)], [(3, solid_copy)],
                                           removeObject=True, removeTool=True)
        occ.synchronize()
        sec_area = sum(occ.getMass(2, t) for d, t in section_result if d == 2)
        sec_pct = 100 * (sec_area - EXPECTED_XSEC_AREA) / EXPECTED_XSEC_AREA
        print(f"    area={sec_area:.4f} mm^2 (%diff={sec_pct:+.4f}%)")
        xsec_ok = abs(sec_pct) < 0.5
    except Exception as e:
        print(f"    FAILED: {e}"); xsec_ok = False
    print(f"  Sectional check: [{'OK' if xsec_ok else 'MISMATCH -- STOP'}]")

    if section_result:
        occ.remove(section_result, recursive=True); occ.synchronize()

    n_vols = len(gmsh.model.getEntities(3))
    all_faces = len(gmsh.model.getEntities(2))
    solid_faces = None  # gmsh has no "faces of this solid" query without a boundary call
    boundary = gmsh.model.getBoundary([(3, vol_tag)], oriented=False)
    solid_faces = len(boundary)
    stray_ok = (n_vols == 1 and all_faces == solid_faces)
    print(f"\n  Model cleanliness: {n_vols} volume(s), {all_faces} surface(s) in model, "
          f"{solid_faces} surface(s) on the solid boundary [{'OK, no stray faces' if stray_ok else 'MISMATCH -- STOP'}]")

    all_ok = sliver_ok and xsec_ok and stray_ok
    print("\n" + "=" * 70)
    print(f"OVERALL VALIDATION: {'PASS' if all_ok else 'FAIL'}")
    print("=" * 70)

    if all_ok:
        brep_path = os.path.join(OUT_DIR, f"{tag_name}.brep")
        step_path = os.path.join(OUT_DIR, f"{tag_name}.step")
        gmsh.write(brep_path)
        gmsh.write(step_path)
        print(f"\nSaved: {brep_path}\nSaved: {step_path}")
        print(f"Mass at Al 6061-T6 density: {actual_vol * RHO * 1000:.2f} g")
    else:
        print("\nGeometry NOT saved -- validation failed.")

    gmsh.finalize()
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
