"""
Stage 6 -- generic Netgen mesh builder for one parametric design point.

Generalizes mesh_netgen_level1_v2.py (Stage 5) to work across the
Stage 6 design space, where the three feature-face areas are NOT
fixed constants -- they depend on the design point's (t, L_gusset,
r_toe). Rather than hardcoded expected areas or fragile bbox-location
heuristics, this script computes the expected area of each feature
face ANALYTICALLY from the design parameters, then matches faces
against those computed values (with the same "abort on mismatch,
don't guess" discipline as Stage 5).

DERIVATIONS (verified exact against Stage 5's hardcoded values at the
baseline design t=4, L_gusset=20, r_toe=1, R3=3):

  R3 fillet area = R3 * (pi/2) * WIDTH
      The Flange A / Flange B junction is always a true 90 degree
      corner regardless of t, so this is a CONSTANT across the whole
      sweep (188.4956 mm^2 at R3=3, WIDTH=40).

  Toe fillet area = r_toe * (pi/4) * WIDTH
      The angle at the Toe vertex is always exactly 45 degrees,
      independent of L_gusset -- because both Toe and Far2 are
      offset by the SAME L_gusset from the corner point, so the
      Toe->Far2 edge always has dx=dz=L_gusset (a 45 degree line).
      Scales linearly with r_toe (31.4159 mm^2 at r_toe=1).

  Hole bore area = 2*pi*HOLE_R*t
      The hole passes through Flange A's full thickness t.
      Scales linearly with t (62.8319 mm^2 at t=4).

USAGE (single design point):
    from mesh_stage6_design import build_mesh
    build_mesh(design_id, t, L_gusset, r_toe, hole_h, r3_h, toe_h,
               brep_path, inp_out)
"""

from netgen.occ import OCCGeometry
import numpy as np
import os
import sys
import math

WIDTH = 40.0
R3 = 3.0
HOLE_D = 5.0
HOLE_R = HOLE_D / 2.0

E_MOD = 71700.0
NU = 0.33
RHO = 2810e-9
FORCE_X = -235.44
REF_NODE_COORD = (60.0, 20.0, 63.0)  # arbitrary coupling reference point, unrelated to geometry
MAXH_GLOBAL = 6.0

AREA_TOL = 0.05  # mm^2, tolerance for analytic-area face matching
CEILING = 1_180_000  # documented direct-solver equation ceiling


def expected_areas(t, L_gusset, r_toe):
    r3_area = R3 * (math.pi / 2.0) * WIDTH
    toe_area = r_toe * (math.pi / 4.0) * WIDTH
    hole_area = 2.0 * math.pi * HOLE_R * t
    return {"R3": r3_area, "toe": toe_area, "hole": hole_area}


def identify_faces(shape, t, L_gusset, r_toe):
    """
    Returns dict {"R3": idx, "toe": idx, "hole": idx} or raises with a
    clear message if any category has zero or more than one match.
    """
    targets = expected_areas(t, L_gusset, r_toe)
    faces = list(shape.faces)
    areas = [f.mass for f in faces]

    found = {}
    for label, target_area in targets.items():
        matches = [i for i, a in enumerate(areas) if abs(a - target_area) < AREA_TOL]
        if len(matches) == 0:
            raise RuntimeError(
                f"Face identification FAILED for '{label}': expected area "
                f"{target_area:.4f} mm^2, no face matched within tolerance "
                f"{AREA_TOL}. All face areas: {[f'{a:.4f}' for a in areas]}"
            )
        if len(matches) > 1:
            raise RuntimeError(
                f"Face identification AMBIGUOUS for '{label}': expected area "
                f"{target_area:.4f} mm^2, {len(matches)} faces matched: "
                f"{matches} with areas {[areas[m] for m in matches]}"
            )
        found[label] = matches[0]

    return found, targets


def build_mesh(design_id, t, L_gusset, r_toe, hole_h, r3_h, toe_h, brep_path, inp_out):
    """
    Returns a dict of diagnostics (always), with 'success': bool.
    Never raises -- all failures are captured and returned so a batch
    runner can log every design point without one failure halting the
    whole sweep.
    """
    diag = {
        "design_id": design_id, "t": t, "L_gusset": L_gusset, "r_toe": r_toe,
        "hole_h": hole_h, "r3_h": r3_h, "toe_h": toe_h,
        "success": False,
    }

    print(f"\n{'='*70}\n[{design_id}] mesh: t={t}, L_gusset={L_gusset}, r_toe={r_toe} | "
          f"sizing: hole={hole_h}, R3={r3_h}, toe={toe_h}\n{'='*70}")

    if not os.path.exists(brep_path):
        diag["failure_stage"] = "brep_missing"
        diag["failure_reason"] = f"BRep not found: {brep_path}"
        print(f"[{design_id}] ERROR: {diag['failure_reason']}")
        return diag

    try:
        geo = OCCGeometry(brep_path)
        shape = geo.shape
    except Exception as e:
        diag["failure_stage"] = "brep_load"
        diag["failure_reason"] = str(e)
        print(f"[{design_id}] ERROR loading BRep: {e}")
        return diag

    try:
        found_faces, target_areas = identify_faces(shape, t, L_gusset, r_toe)
    except RuntimeError as e:
        diag["failure_stage"] = "face_identification"
        diag["failure_reason"] = str(e)
        print(f"[{design_id}] ERROR: {e}")
        return diag

    print(f"[{design_id}] Face IDs: hole={found_faces['hole']} (expect area "
          f"{target_areas['hole']:.4f}), R3={found_faces['R3']} (expect "
          f"{target_areas['R3']:.4f}), toe={found_faces['toe']} (expect "
          f"{target_areas['toe']:.4f})")

    geo.SetFaceMeshsize(found_faces["hole"], hole_h)
    geo.SetFaceMeshsize(found_faces["R3"], r3_h)
    geo.SetFaceMeshsize(found_faces["toe"], toe_h)

    try:
        mesh = geo.GenerateMesh(maxh=MAXH_GLOBAL)
    except Exception as e:
        diag["failure_stage"] = "mesh_generation"
        diag["failure_reason"] = str(e)
        print(f"[{design_id}] ERROR during mesh generation: {e}")
        return diag

    points = mesh.Points()
    n_pts = len(points)
    coords = np.zeros((n_pts + 1, 3))
    for i, p in enumerate(points, start=1):
        pnt = p.p
        coords[i] = [pnt[0], pnt[1], pnt[2]]

    elements = mesh.Elements3D()
    tets = []
    for el in elements:
        pts = [v.nr for v in el.vertices]
        if len(pts) != 4:
            continue
        tets.append(pts)
    tets = np.array(tets)
    tets = tets[:, [0, 1, 3, 2]]  # Netgen->CalculiX orientation fix, confirmed Stage 5

    def signed_volume(p0, p1, p2, p3):
        return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

    vols = np.array([
        signed_volume(coords[tt[0]], coords[tt[1]], coords[tt[2]], coords[tt[3]])
        for tt in tets
    ])
    n_bad = int((vols <= 0).sum())
    diag["n_linear_tets"] = len(tets)
    diag["min_linear_volume"] = float(vols.min()) if len(vols) else None
    diag["mean_linear_volume"] = float(vols.mean()) if len(vols) else None
    diag["n_bad_volume_elements"] = n_bad

    print(f"[{design_id}] Linear tets: {len(tets)}, min vol={vols.min():.6f}, "
          f"mean vol={vols.mean():.6f}, bad={n_bad}")

    if n_bad > 0:
        diag["failure_stage"] = "negative_volume_elements"
        diag["failure_reason"] = f"{n_bad} elements with non-positive linear volume"
        print(f"[{design_id}] ERROR: {diag['failure_reason']}")
        return diag

    # Promote to C3D10 (straight-sided, deterministic midside nodes)
    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    midside_cache = {}
    next_node_tag = n_pts + 1
    extra_coords = {}
    quad_elems = []
    for tt in tets:
        full_conn = list(tt)
        for (a, b) in edge_pairs:
            na, nb = tt[a], tt[b]
            key = (min(na, nb), max(na, nb))
            if key not in midside_cache:
                mid = (coords[na] + coords[nb]) / 2.0
                extra_coords[next_node_tag] = mid
                midside_cache[key] = next_node_tag
                next_node_tag += 1
            full_conn.append(midside_cache[key])
        quad_elems.append(full_conn)

    total_nodes = n_pts + len(midside_cache)
    diag["n_nodes"] = total_nodes
    diag["n_elements"] = len(quad_elems)

    all_node_tags = list(range(1, n_pts + 1)) + list(extra_coords.keys())

    def get_coord(tag):
        return coords[tag] if tag <= n_pts else extra_coords[tag]

    root_nodes = [tg for tg in all_node_tags if abs(get_coord(tg)[0] - 0.0) < 1e-3]
    tip_nodes = [tg for tg in all_node_tags if abs(get_coord(tg)[2] - 60.0) < 1e-3]
    diag["n_root_nodes"] = len(root_nodes)
    diag["n_tip_nodes"] = len(tip_nodes)

    if len(root_nodes) == 0 or len(tip_nodes) == 0:
        diag["failure_stage"] = "boundary_node_identification"
        diag["failure_reason"] = (
            f"root_nodes={len(root_nodes)}, tip_nodes={len(tip_nodes)} (both must be > 0)"
        )
        print(f"[{design_id}] ERROR: {diag['failure_reason']}")
        return diag

    est_equations = 3 * total_nodes
    diag["est_equations"] = est_equations
    diag["pct_of_ceiling"] = 100 * est_equations / CEILING
    print(f"[{design_id}] Nodes={total_nodes}, Elements={len(quad_elems)}, "
          f"Equations~{est_equations} ({diag['pct_of_ceiling']:.1f}% of ceiling)")

    if est_equations > CEILING:
        diag["failure_stage"] = "exceeds_solver_ceiling"
        diag["failure_reason"] = (
            f"Estimated {est_equations} equations exceeds documented "
            f"{CEILING} direct-solver ceiling"
        )
        print(f"[{design_id}] ERROR: {diag['failure_reason']}")
        return diag

    os.makedirs(os.path.dirname(inp_out), exist_ok=True)
    ref_tag = total_nodes + 1
    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Stage 6 design {design_id} -- t={t}, L_gusset={L_gusset}, r_toe={r_toe}\n")
        f.write("Al 7075-T6, F=235.44N in -x at Flange B tip, root encastre\n")

        f.write("*NODE\n")
        for tag in range(1, n_pts + 1):
            c = coords[tag]
            f.write(f"{tag}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")
        for tag, c in extra_coords.items():
            f.write(f"{tag}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")
        f.write(f"{ref_tag}, {REF_NODE_COORD[0]}, {REF_NODE_COORD[1]}, {REF_NODE_COORD[2]}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=EBRACKET\n")
        for i, conn in enumerate(quad_elems, start=1):
            f.write(f"{i}, " + ", ".join(str(n) for n in conn) + "\n")

        f.write("*NSET, NSET=NROOT\n")
        for i in range(0, len(root_nodes), 10):
            f.write(", ".join(str(n) for n in root_nodes[i:i + 10]) + "\n")
        f.write("*NSET, NSET=NTIP\n")
        for i in range(0, len(tip_nodes), 10):
            f.write(", ".join(str(n) for n in tip_nodes[i:i + 10]) + "\n")

        f.write("*SURFACE, NAME=STIP, TYPE=NODE\nNTIP\n")
        f.write("*MATERIAL, NAME=AL7075T6\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=EBRACKET, MATERIAL=AL7075T6\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNROOT, 1, 3, 0.0\n")
        f.write(f"*COUPLING, CONSTRAINT NAME=CTIP, REF NODE={ref_tag}, SURFACE=STIP\n")
        f.write("*KINEMATIC\n1, 3\n")
        f.write(f"*CLOAD\n{ref_tag}, 1, {FORCE_X}\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")

    diag["success"] = True
    diag["inp_path"] = inp_out
    print(f"[{design_id}] PASSED mesh quality checks. Wrote: {inp_out}")
    return diag


if __name__ == "__main__":
    if len(sys.argv) < 8:
        print("Usage: python3 mesh_stage6_design.py <design_id> <t> <L_gusset> <r_toe> "
              "<hole_h> <r3_h> <toe_h> [brep_path] [inp_out]")
        sys.exit(1)

    design_id = sys.argv[1]
    t = float(sys.argv[2])
    L_gusset = float(sys.argv[3])
    r_toe = float(sys.argv[4])
    hole_h = float(sys.argv[5])
    r3_h = float(sys.argv[6])
    toe_h = float(sys.argv[7])
    brep_path = sys.argv[8] if len(sys.argv) > 8 else f"mesh/stage6_parametric/stage6_design_{design_id}.brep"
    inp_out = sys.argv[9] if len(sys.argv) > 9 else f"simulation/stage6_parametric/stage6_{design_id}_mesh.inp"

    result = build_mesh(design_id, t, L_gusset, r_toe, hole_h, r3_h, toe_h, brep_path, inp_out)
    print(f"\n{'='*70}\nRESULT: {'PASS' if result['success'] else 'FAIL'}\n{'='*70}")
    if not result["success"]:
        print(f"Failure stage: {result.get('failure_stage')}")
        print(f"Failure reason: {result.get('failure_reason')}")
        sys.exit(1)
