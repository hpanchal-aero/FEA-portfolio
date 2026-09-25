"""
Project 02, Configuration B -- randomly generated pocket panel geometry.

FINAL CONFIG (see development-log.md for the exploration that led here):
N_POINTS_MIN=5, N_POINTS_MAX=12 (cycled across attempts), N_ATTEMPTS_MAX=150,
no cap on candidates kept (223 valid candidates were produced from 150
attempts in the actual run used for this project).

Generates candidate irregular polygons across N_ATTEMPTS_MAX attempts,
cycling the point count from N_POINTS_MIN to N_POINTS_MAX (5..12). Each
attempt places that many random points inside the panel's interior region
and connects them two ways -- convex hull, and angular/star order around
the centroid -- both guaranteed simple/non-self-intersecting by
construction. No cap on the number of valid results kept; every candidate
that passes all checks is built and validated.

Validation per candidate: exact target area via uniform scaling from the
centroid (shoelace formula), minimum edge length >= MIN_FEATURE_MM, fits
within the panel's interior region, no duplicate (rounded-vertex) shape.
Each accepted polygon is then built as a through-cut pocket panel
(matching build_config_b_panels.py's isolated-panel methodology: box cut
by an extruded polygon prism) and validated the same way as the 4 named
shapes (area match, no slivers, correct bbox, single-solid cleanliness).

Usage: /usr/bin/python3 scripts/build_config_b_random.py

Fixed random seed for reproducibility.
"""
import csv
import math
import os
import random

import gmsh

N_POINTS_MIN = 5
N_POINTS_MAX = 12
N_ATTEMPTS_MAX = 150
SEED = 2

PANEL_W = 83.0
PANEL_H = 100.0
PANEL_T = 2.0
FRAME_W = 5.0
LEVEL_PCT = 10
TARGET_AREA = LEVEL_PCT / 100.0 * PANEL_W * PANEL_H  # 830.0 mm^2

W_INT = PANEL_W - 2 * FRAME_W  # 73.0
H_INT = PANEL_H - 2 * FRAME_W  # 90.0

MIN_FEATURE_MM = 3.0
SLIVER_THRESHOLD = 1.0

OUT_DIR = "mesh/config_b_study"
os.makedirs(OUT_DIR, exist_ok=True)
RESULTS_CSV = "results/config_b_random_geometry_validation.csv"
os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)


def polygon_area(pts):
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def convex_hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return None

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = lower[:-1] + upper[:-1]
    return hull if len(hull) >= 3 else None


def angular_order(pts):
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    return sorted(pts, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))


def segments_intersect(p1, p2, p3, p4):
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
    return (ccw(p1, p3, p4) != ccw(p2, p3, p4)) and (ccw(p1, p2, p3) != ccw(p1, p2, p4))


def is_simple_polygon(pts):
    n = len(pts)
    if n < 3:
        return False
    edges = [(pts[i], pts[(i + 1) % n]) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue
            if segments_intersect(*edges[i], *edges[j]):
                return False
    return True


def min_edge_length(pts):
    n = len(pts)
    return min(math.dist(pts[i], pts[(i + 1) % n]) for i in range(n))


def scale_to_area(pts, target_area):
    area = polygon_area(pts)
    if area < 1e-6:
        return None
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    k = math.sqrt(target_area / area)
    return [(cx + k * (x - cx), cy + k * (y - cy)) for x, y in pts]


def fits_in_frame(pts):
    return all(abs(x) <= W_INT / 2 and abs(y) <= H_INT / 2 for x, y in pts)


def generate_candidates():
    random.seed(SEED)
    candidates = []
    seen = set()
    attempts = 0
    n_range = list(range(N_POINTS_MIN, N_POINTS_MAX + 1))
    while attempts < N_ATTEMPTS_MAX:
        n_pts = n_range[attempts % len(n_range)]
        attempts += 1
        raw_pts = [
            (random.uniform(-W_INT / 2, W_INT / 2), random.uniform(-H_INT / 2, H_INT / 2))
            for _ in range(n_pts)
        ]
        for method_name, order_fn in (("hull", convex_hull), ("star", angular_order)):
            ordered = order_fn(raw_pts)
            if ordered is None or not is_simple_polygon(ordered):
                continue
            scaled = scale_to_area(ordered, TARGET_AREA)
            if scaled is None:
                continue
            if not is_simple_polygon(scaled):
                continue
            if not fits_in_frame(scaled):
                continue
            if min_edge_length(scaled) < MIN_FEATURE_MM:
                continue
            key = tuple(round(x, 3) for pt in scaled for x in pt)
            if key in seen:
                continue
            seen.add(key)
            candidates.append((f"attempt{attempts}_{method_name}", scaled))
    return candidates, attempts


def build_and_validate(tag, poly_pts, results):
    print(f"\n{'='*70}\n[{tag}] N_vertices={len(poly_pts)}\n{'='*70}")
    print(f"  vertices (mm): {[(round(x,2), round(y,2)) for x,y in poly_pts]}")
    area = polygon_area(poly_pts)
    min_edge = min_edge_length(poly_pts)
    print(f"  area={area:.4f} mm^2 (target {TARGET_AREA}), min_edge={min_edge:.4f} mm")

    gmsh.initialize()
    gmsh.model.add(tag)
    occ = gmsh.model.occ

    panel = occ.addBox(-PANEL_W / 2, -PANEL_T / 2, -PANEL_H / 2, PANEL_W, PANEL_T, PANEL_H)
    occ.synchronize()
    panel_vol_pre = occ.getMass(3, panel)

    pt_tags = [occ.addPoint(x, -PANEL_T * 2, y) for x, y in poly_pts]
    line_tags = [occ.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)]) for i in range(len(pt_tags))]
    loop = occ.addCurveLoop(line_tags)
    surf = occ.addPlaneSurface([loop])
    occ.synchronize()
    extruded = occ.extrude([(2, surf)], 0, PANEL_T * 4, 0)
    occ.synchronize()
    cutter_vol = [t for d, t in extruded if d == 3]

    result, _ = occ.cut([(3, panel)], [(3, t) for t in cutter_vol])
    occ.synchronize()

    row = {"tag": tag, "n_vertices": len(poly_pts), "target_area": TARGET_AREA,
           "input_area": area, "min_edge_mm": min_edge}

    if len(result) != 1:
        print(f"  [{tag}] FATAL: cut produced {len(result)} volumes, expected 1.")
        row.update({"status": "multi_volume_fail", "pass": False})
        results.append(row)
        gmsh.finalize()
        return

    vol_tag = result[0][1]
    post_vol = occ.getMass(3, vol_tag)
    removed_vol = panel_vol_pre - post_vol
    removed_area = removed_vol / PANEL_T
    area_pct_diff = 100 * (removed_area - TARGET_AREA) / TARGET_AREA
    print(f"  post-cut volume={post_vol:.4f}, removed_area={removed_area:.4f} "
          f"(%diff={area_pct_diff:+.4f}%)")
    area_ok = abs(area_pct_diff) < 1.0

    surfaces = gmsh.model.getEntities(2)
    areas = [(occ.getMass(2, t), t) for d, t in surfaces]
    slivers = [(a, t) for a, t in areas if a < SLIVER_THRESHOLD]
    sliver_ok = len(slivers) == 0
    print(f"  surfaces={len(surfaces)}, slivers={len(slivers)}")

    bbox = gmsh.model.getBoundingBox(3, vol_tag)
    bbox_ok = (abs(bbox[0] + PANEL_W / 2) < 1e-3 and abs(bbox[3] - PANEL_W / 2) < 1e-3 and
               abs(bbox[2] + PANEL_H / 2) < 1e-3 and abs(bbox[5] - PANEL_H / 2) < 1e-3)
    print(f"  bbox: [{'OK' if bbox_ok else 'FAIL'}]")

    n_vols = len(gmsh.model.getEntities(3))
    boundary = gmsh.model.getBoundary([(3, vol_tag)], oriented=False)
    stray_ok = (n_vols == 1 and len(surfaces) == len(boundary))
    print(f"  cleanliness: {n_vols} volume(s), {len(surfaces)} surf in model, "
          f"{len(boundary)} on solid boundary [{'OK' if stray_ok else 'FAIL'}]")

    all_ok = area_ok and sliver_ok and bbox_ok and stray_ok
    row.update({"status": "ok" if all_ok else "validation_fail", "pass": all_ok,
                "post_vol": post_vol, "removed_area": removed_area,
                "area_pct_diff": area_pct_diff, "n_surfaces": len(surfaces),
                "n_slivers": len(slivers)})
    results.append(row)

    if all_ok:
        out_path = os.path.join(OUT_DIR, f"config_b_random_{tag}.brep")
        gmsh.write(out_path)
        vtx_path = os.path.join(OUT_DIR, f"config_b_random_{tag}_vertices.json")
        with open(vtx_path, "w") as vf:
            import json
            json.dump({"tag": tag, "n_vertices": len(poly_pts), "vertices_mm": poly_pts}, vf)
        print(f"  [{tag}] PASSED. Saved: {out_path} and {vtx_path}")
    else:
        print(f"  [{tag}] FAILED validation -- not saved.")

    gmsh.finalize()


def main():
    candidates, attempts = generate_candidates()
    print(f"Generated {len(candidates)} valid candidate polygons from {attempts} attempts "
          f"(no target cap -- all valid candidates kept)")

    results = []
    for tag, pts in candidates:
        build_and_validate(tag, pts, results)

    fieldnames = ["tag", "n_vertices", "target_area", "input_area", "min_edge_mm",
                  "status", "pass", "post_vol", "removed_area", "area_pct_diff",
                  "n_surfaces", "n_slivers"]
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    n_pass = sum(1 for r in results if r.get("pass"))
    print(f"\n{'='*70}\nBATCH SUMMARY: {n_pass}/{len(results)} PASSED\n{'='*70}")
    print(f"Saved: {RESULTS_CSV}")


if __name__ == "__main__":
    main()
