"""
Project 02, Configuration B -- lightened panel geometry generation
and validation. 4 pocket shapes x 5 removal levels = 20 geometries,
each a through-cut pocket in an isolated 83x100x2mm panel with a
5mm perimeter frame retained at full thickness.

Fairness constraint: at a given removal level, all 4 shapes target
the SAME removed area (not just "looks similar") -- verified
computationally per geometry, not assumed from the input formula.

Geometry construction and validation ONLY. Does not mesh or solve.
"""

import gmsh
import math
import csv
import os
import sys

PANEL_W = 83.0   # x, mm
PANEL_H = 100.0  # z, mm  (matches panel height in the assembled frame)
PANEL_T = 2.0    # y, mm (thickness)
FRAME_W = 5.0    # mm, perimeter frame width retained at full thickness
RIB = 8.0        # mm, inter-pocket rib width for the grid shape

W_INT = PANEL_W - 2 * FRAME_W   # 73.0
H_INT = PANEL_H - 2 * FRAME_W   # 90.0
PANEL_AREA = PANEL_W * PANEL_H  # 8300

LEVELS = [10, 15, 20, 25, 30]
SLIVER_THRESHOLD = 1.0  # mm^2

OUT_DIR = "mesh/config_b_study"
os.makedirs(OUT_DIR, exist_ok=True)
RESULTS_CSV = "results/config_b_geometry_validation.csv"
os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)


def rect_dims(A):
    w = math.sqrt(A * W_INT / H_INT)
    h = math.sqrt(A * H_INT / W_INT)
    return w, h


def circle_dims(A):
    r = math.sqrt(A / math.pi)
    return r


def cross_dims(A):
    b = -(W_INT + H_INT)
    c = A
    disc = b * b - 4 * c
    if disc < 0:
        return None
    w_arm = (-b - math.sqrt(disc)) / 2
    return w_arm


def grid_dims(A):
    side = math.sqrt(A / 4)
    footprint = 2 * side + RIB
    return side, footprint


def build_rect_pocket(occ, cut_from, A):
    w, h = rect_dims(A)
    cutter = occ.addBox(-w/2, -PANEL_T, -h/2, w, PANEL_T * 3, h)
    occ.synchronize()
    result, _ = occ.cut([(3, cut_from)], [(3, cutter)])
    occ.synchronize()
    return result, {"rect_w": w, "rect_h": h}


def build_circle_pocket(occ, cut_from, A):
    r = circle_dims(A)
    cutter = occ.addCylinder(0, -PANEL_T, 0, 0, PANEL_T * 3, 0, r)
    occ.synchronize()
    result, _ = occ.cut([(3, cut_from)], [(3, cutter)])
    occ.synchronize()
    return result, {"circ_r": r}


def build_cross_pocket(occ, cut_from, A):
    w_arm = cross_dims(A)
    if w_arm is None:
        return None, {"cross_w": None}
    arm1 = occ.addBox(-W_INT/2, -PANEL_T, -w_arm/2, W_INT, PANEL_T * 3, w_arm)
    arm2 = occ.addBox(-w_arm/2, -PANEL_T, -H_INT/2, w_arm, PANEL_T * 3, H_INT)
    occ.synchronize()
    fused, _ = occ.fuse([(3, arm1)], [(3, arm2)])
    occ.synchronize()
    cutter_tag = fused[0][1]
    result, _ = occ.cut([(3, cut_from)], [(3, cutter_tag)])
    occ.synchronize()
    return result, {"cross_w": w_arm}


def build_grid_pocket(occ, cut_from, A):
    side, footprint = grid_dims(A)
    if footprint > min(W_INT, H_INT):
        return None, {"grid_side": None}
    offset = side / 2 + RIB / 2
    cutters = []
    for sx in (-1, 1):
        for sz in (-1, 1):
            cx = sx * offset
            cz = sz * offset
            c = occ.addBox(cx - side/2, -PANEL_T, cz - side/2, side, PANEL_T * 3, side)
            cutters.append(c)
    occ.synchronize()
    result = [(3, cut_from)]
    for c in cutters:
        result, _ = occ.cut(result, [(3, c)])
        occ.synchronize()
        cut_from = result[0][1]
    return result, {"grid_side": side, "grid_footprint": footprint}


SHAPE_BUILDERS = {
    "rect": build_rect_pocket,
    "circle": build_circle_pocket,
    "cross": build_cross_pocket,
    "grid": build_grid_pocket,
}


def build_and_validate(shape, level, results):
    A_target = level / 100 * PANEL_AREA
    tag = f"{shape}_{level}pct"
    print(f"\n{'='*70}\n[{tag}] target removed area = {A_target:.2f} mm^2\n{'='*70}")

    gmsh.initialize()
    gmsh.model.add(tag)
    occ = gmsh.model.occ

    panel = occ.addBox(-PANEL_W/2, -PANEL_T/2, -PANEL_H/2, PANEL_W, PANEL_T, PANEL_H)
    occ.synchronize()
    panel_vol_pre = occ.getMass(3, panel)
    print(f"  Pre-cut panel volume: {panel_vol_pre:.4f} mm^3 (expect {PANEL_W*PANEL_T*PANEL_H:.4f})")

    builder = SHAPE_BUILDERS[shape]
    result, dims = builder(occ, panel, A_target)

    row = {"shape": shape, "level_pct": level, "target_area": A_target, **dims}

    if result is None:
        print(f"  [{tag}] INFEASIBLE at this level -- geometry not built.")
        row.update({"status": "infeasible", "pass": False})
        results.append(row)
        gmsh.finalize()
        return

    vols = gmsh.model.getEntities(3)
    if len(vols) != 1:
        print(f"  [{tag}] FATAL: cut produced {len(vols)} volumes, expected 1.")
        row.update({"status": "multi_volume_fail", "pass": False})
        results.append(row)
        gmsh.finalize()
        return

    vol_tag = vols[0][1]
    post_vol = occ.getMass(3, vol_tag)
    removed_vol = panel_vol_pre - post_vol
    removed_area = removed_vol / PANEL_T  # through-cut, so removed_vol = removed_area * thickness
    area_pct_diff = 100 * (removed_area - A_target) / A_target
    print(f"  Post-cut volume: {post_vol:.4f} mm^3, removed_vol={removed_vol:.4f}, "
          f"removed_area={removed_area:.4f} mm^2 (target {A_target:.4f}, %diff={area_pct_diff:+.4f}%)")
    area_ok = abs(area_pct_diff) < 1.0

    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, t), t) for d, t in surfaces]
    slivers = [(a, t) for a, t in areas if a < SLIVER_THRESHOLD]
    print(f"  Surfaces: {len(surfaces)}, slivers (<{SLIVER_THRESHOLD}mm^2): {len(slivers)}")
    sliver_ok = len(slivers) == 0

    bbox = gmsh.model.getBoundingBox(3, vol_tag)
    bbox_ok = (abs(bbox[0] + PANEL_W/2) < 1e-3 and abs(bbox[3] - PANEL_W/2) < 1e-3 and
               abs(bbox[2] + PANEL_H/2) < 1e-3 and abs(bbox[5] - PANEL_H/2) < 1e-3)
    print(f"  Bbox check (outer envelope preserved): [{'OK' if bbox_ok else 'FAIL'}]")

    all_ok = area_ok and sliver_ok and bbox_ok
    row.update({
        "status": "ok" if all_ok else "validation_fail",
        "pass": all_ok,
        "post_vol": post_vol,
        "removed_area": removed_area,
        "area_pct_diff": area_pct_diff,
        "n_surfaces": len(surfaces),
        "n_slivers": len(slivers),
    })
    results.append(row)

    if all_ok:
        out_path = os.path.join(OUT_DIR, f"config_b_{tag}.brep")
        gmsh.write(out_path)
        print(f"  [{tag}] PASSED. Saved: {out_path}")
    else:
        print(f"  [{tag}] FAILED validation.")

    gmsh.finalize()


def main():
    results = []
    for shape in SHAPE_BUILDERS:
        for level in LEVELS:
            build_and_validate(shape, level, results)

    fieldnames = ["shape", "level_pct", "target_area", "status", "pass",
                  "post_vol", "removed_area", "area_pct_diff",
                  "n_surfaces", "n_slivers",
                  "rect_w", "rect_h", "circ_r", "cross_w", "grid_side", "grid_footprint"]
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    n_pass = sum(1 for r in results if r.get("pass"))
    print(f"\n{'='*70}\nBATCH SUMMARY: {n_pass}/{len(results)} PASSED\n{'='*70}")
    print(f"{'shape':<10}{'level':<8}{'status':<20}{'area_%diff':<12}{'slivers'}")
    for r in results:
        print(f"{r['shape']:<10}{r['level_pct']:<8}{r.get('status',''):<20}"
              f"{r.get('area_pct_diff', 0):<12.4f}{r.get('n_slivers','')}")
    print(f"\nSaved: {RESULTS_CSV}")


if __name__ == "__main__":
    main()
