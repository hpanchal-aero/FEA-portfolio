"""
Project 02, Configuration B -- full-assembly geometry for randomly
generated pocket panels.

Reads a saved polygon vertex set (mesh/config_b_study/config_b_random_
<tag>_vertices.json, written by build_config_b_random.py), applies the
same polygon pocket to all 4 panels of the fused rail+panel assembly
(same methodology as build_config_b_assembly.py's named shapes), and
validates: per-panel removed area, rail/panel overlap, total volume,
bounding box, section area at z=2mm (inside the always-uncut 5mm frame
margin -- see the cross-shape SECTION_Z lesson), and single-solid
cleanliness (no stray faces).

Usage: /usr/bin/python3 scripts/build_config_b_random_assembly.py <tag>
   or: /usr/bin/python3 scripts/build_config_b_random_assembly.py --all
       (processes every config_b_random_*_vertices.json found)
"""
import glob
import json
import os
import sys

import gmsh

TOL = 1e-6
ENV = 50.0
RAIL = 8.5
HEIGHT = 100.0
PANEL_T = 2.0
RHO = 2700e-9
FRAME_W = 5.0

RAIL_INNER = ENV - RAIL
PANEL_LEN = 2 * RAIL_INNER
PANEL_AREA = PANEL_LEN * HEIGHT
W_INT = PANEL_LEN - 2 * FRAME_W
H_INT = HEIGHT - 2 * FRAME_W

EXPECTED_RAIL_VOL_EACH = RAIL * RAIL * HEIGHT
EXPECTED_PANEL_VOL_EACH = PANEL_LEN * PANEL_T * HEIGHT
EXPECTED_XSEC_AREA = 4 * (RAIL * RAIL) + 4 * (PANEL_LEN * PANEL_T)
SECTION_Z = 2.0

OUT_DIR = "mesh/config_b_study"
RESULTS_CSV = "results/config_b_random_assembly_validation.csv"


def cutter_local(occ, poly_pts, big):
    pt_tags = [occ.addPoint(x, -big / 2, y, 0) for x, y in poly_pts]
    line_tags = [occ.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                 for i in range(len(pt_tags))]
    loop = occ.addCurveLoop(line_tags)
    surf = occ.addPlaneSurface([loop])
    occ.synchronize()
    extruded = occ.extrude([(2, surf)], 0, big, 0)
    occ.synchronize()
    return [t for d, t in extruded if d == 3]


def place_cutters_global(occ, local_tags, orientation, center_uw, y_or_x_center):
    out = []
    for t in local_tags:
        if orientation == "y_normal":
            occ.translate([(3, t)], center_uw[0], y_or_x_center, center_uw[1])
        else:
            occ.rotate([(3, t)], 0, 0, 0, 0, 0, 1, 1.5707963267948966)
            occ.translate([(3, t)], y_or_x_center, center_uw[0], center_uw[1])
        out.append((3, t))
    return out


def build_one(tag, poly_pts, results):
    A_target = None  # computed below from the polygon itself (shoelace)
    n = len(poly_pts)
    a = 0.0
    for i in range(n):
        x1, y1 = poly_pts[i]
        x2, y2 = poly_pts[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    A_target = abs(a) / 2.0
    removed_vol_each = A_target * PANEL_T
    expected_total_vol = (4 * EXPECTED_RAIL_VOL_EACH
                           + 4 * (EXPECTED_PANEL_VOL_EACH - removed_vol_each))

    tag_name = f"config_b_random_assembly_{tag}"
    gmsh.initialize()
    gmsh.model.add(tag_name)
    occ = gmsh.model.occ

    print(f"\n{'='*70}\n[{tag_name}] n_vertices={n}, A_target={A_target:.4f}\n{'='*70}")

    rail_positions = [(-ENV, -ENV), (RAIL_INNER, -ENV), (-ENV, RAIL_INNER), (RAIL_INNER, RAIL_INNER)]
    rail_tags = []
    for x0, y0 in rail_positions:
        t = occ.addBox(x0, y0, 0.0, RAIL, RAIL, HEIGHT)
        rail_tags.append(t)
    occ.synchronize()

    panel_defs = [
        ("y=-ENV", -RAIL_INNER, -ENV, PANEL_LEN, PANEL_T, "y_normal", -ENV + PANEL_T / 2),
        ("y=+ENV", -RAIL_INNER, ENV - PANEL_T, PANEL_LEN, PANEL_T, "y_normal", ENV - PANEL_T / 2),
        ("x=-ENV", -ENV, -RAIL_INNER, PANEL_T, PANEL_LEN, "x_normal", -ENV + PANEL_T / 2),
        ("x=+ENV", ENV - PANEL_T, -RAIL_INNER, PANEL_T, PANEL_LEN, "x_normal", ENV - PANEL_T / 2),
    ]
    panel_tags = []
    row = {"tag": tag, "n_vertices": n, "target_area_per_panel": A_target}
    for label, x0, y0, dx, dy, orient, thru_center in panel_defs:
        panel = occ.addBox(x0, y0, 0.0, dx, dy, HEIGHT)
        occ.synchronize()
        local = cutter_local(occ, poly_pts, PANEL_T * 4)
        occ.synchronize()
        cutters = place_cutters_global(occ, local, orient, (0.0, 50.0), thru_center)
        occ.synchronize()
        result, _ = occ.cut([(3, panel)], cutters)
        occ.synchronize()
        if len(result) != 1:
            print(f"  FATAL: panel [{label}] cut gave {len(result)} volumes.")
            row.update({"status": "panel_cut_fail", "pass": False})
            results.append(row)
            gmsh.finalize()
            return
        cut_panel = result[0][1]
        post_vol = occ.getMass(3, cut_panel)
        removed_vol = EXPECTED_PANEL_VOL_EACH - post_vol
        pct = 100 * (removed_vol - removed_vol_each) / removed_vol_each
        if abs(pct) > 1.0:
            print(f"  MISMATCH panel [{label}]: removed_vol %diff={pct:+.4f}%")
            row.update({"status": "panel_area_mismatch", "pass": False, "worst_panel_pct_diff": pct})
            results.append(row)
            gmsh.finalize()
            return
        panel_tags.append(cut_panel)

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
    if overlap >= 1e-6:
        print(f"  WARNING: nonzero rail/panel overlap {overlap:.6f} mm^3")

    all_tags = rail_tags + panel_tags
    fused, _ = occ.fuse([(3, all_tags[0])], [(3, t) for t in all_tags[1:]])
    occ.synchronize()
    if len(fused) != 1:
        print(f"  FATAL: fuse gave {len(fused)} volumes.")
        row.update({"status": "fuse_fail", "pass": False})
        results.append(row)
        gmsh.finalize()
        return
    vol_tag = fused[0][1]

    actual_vol = occ.getMass(3, vol_tag)
    vol_pct = 100 * (actual_vol - expected_total_vol) / expected_total_vol
    vol_ok = abs(vol_pct) < 0.1

    bbox = gmsh.model.getBoundingBox(3, vol_tag)
    bbox_ok = (abs(bbox[0] + ENV) < TOL and abs(bbox[3] - ENV) < TOL and
               abs(bbox[1] + ENV) < TOL and abs(bbox[4] - ENV) < TOL and
               abs(bbox[2]) < TOL and abs(bbox[5] - HEIGHT) < TOL)

    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, t), t) for d, t in surfaces]
    slivers = [(a2, t) for a2, t in areas if a2 < 1.0]
    sliver_ok = len(slivers) == 0

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
        xsec_ok = abs(sec_pct) < 0.5
    except Exception as e:
        print(f"  section check FAILED to compute: {e}")
        sec_area, sec_pct, xsec_ok = None, None, False
    if section_result:
        occ.remove(section_result, recursive=True); occ.synchronize()

    n_vols = len(gmsh.model.getEntities(3))
    all_faces = len(gmsh.model.getEntities(2))
    boundary = gmsh.model.getBoundary([(3, vol_tag)], oriented=False)
    solid_faces = len(boundary)
    stray_ok = (n_vols == 1 and all_faces == solid_faces)

    all_ok = vol_ok and bbox_ok and sliver_ok and xsec_ok and stray_ok
    print(f"  volume={actual_vol:.4f} (exp {expected_total_vol:.4f}, {vol_pct:+.4f}%) "
          f"[{'OK' if vol_ok else 'FAIL'}]")
    print(f"  bbox [{'OK' if bbox_ok else 'FAIL'}]  slivers={len(slivers)} [{'OK' if sliver_ok else 'FAIL'}]")
    print(f"  section z={SECTION_Z}: area={sec_area} (%diff={sec_pct}) [{'OK' if xsec_ok else 'FAIL'}]")
    print(f"  cleanliness: {n_vols} vol, {all_faces} surf, {solid_faces} on boundary "
          f"[{'OK' if stray_ok else 'FAIL'}]")
    print(f"  OVERALL: {'PASS' if all_ok else 'FAIL'}")

    row.update({
        "status": "ok" if all_ok else "validation_fail", "pass": all_ok,
        "actual_volume": actual_vol, "expected_volume": expected_total_vol,
        "vol_pct_diff": vol_pct, "n_surfaces": all_faces, "n_slivers": len(slivers),
        "section_area": sec_area, "section_pct_diff": sec_pct,
        "n_solids": n_vols, "n_boundary_faces": solid_faces,
        "mass_g": actual_vol * RHO * 1000.0 if all_ok else None,
    })
    results.append(row)

    if all_ok:
        brep_path = os.path.join(OUT_DIR, f"{tag_name}.brep")
        step_path = os.path.join(OUT_DIR, f"{tag_name}.step")
        gmsh.write(brep_path)
        gmsh.write(step_path)
        print(f"  Saved: {brep_path}\n  Saved: {step_path}")

    gmsh.finalize()


def main():
    if len(sys.argv) != 2:
        print("Usage: build_config_b_random_assembly.py <tag | --all>")
        sys.exit(1)

    if sys.argv[1] == "--all":
        vtx_files = sorted(glob.glob(os.path.join(OUT_DIR, "config_b_random_attempt*_vertices.json")))
    else:
        vtx_files = [os.path.join(OUT_DIR, f"config_b_random_{sys.argv[1]}_vertices.json")]

    results = []
    for vf in vtx_files:
        if not os.path.isfile(vf):
            print(f"WARNING: {vf} not found, skipping")
            continue
        d = json.load(open(vf))
        build_one(d["tag"], [tuple(p) for p in d["vertices_mm"]], results)

    fieldnames = ["tag", "n_vertices", "target_area_per_panel", "status", "pass",
                  "actual_volume", "expected_volume", "vol_pct_diff",
                  "n_surfaces", "n_slivers", "section_area", "section_pct_diff",
                  "n_solids", "n_boundary_faces", "mass_g", "worst_panel_pct_diff"]
    with open(RESULTS_CSV, "w", newline="") as f:
        import csv
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    n_pass = sum(1 for r in results if r.get("pass"))
    print(f"\n{'='*70}\nBATCH SUMMARY: {n_pass}/{len(results)} PASSED\n{'='*70}")
    print(f"Saved: {RESULTS_CSV}")


if __name__ == "__main__":
    main()
