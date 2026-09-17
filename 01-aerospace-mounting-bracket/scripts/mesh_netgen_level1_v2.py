"""
Stage 5 -- Netgen mesh, CORRECTED refinement method.

Original point-based RestrictH() was confirmed too weak: Level 1->2
target size dropped 40% at the toe but actual mean edge length near
the toe only dropped 5% (0.83->0.79mm). Netgen's default grading rate
limits a single point's influence.

Fix: use SetFaceMeshsize() on the actual geometric surfaces (identified
by area fingerprint, same method used throughout this project with
Gmsh) -- a face-level size constraint is a much stronger, more direct
control than a point constraint.

Face indices confirmed by area match against this exact BRep:
  face[2]  = toe fillet   (area 31.4159)
  face[11] = R3 fillet    (area 188.4956)
  face[14] = hole bore    (area 62.8319)

This REPLACES both mesh_netgen.py (Level 1) and mesh_netgen_level2.py --
both are superseded, not reused, since the refinement mechanism itself
was invalid in both.
"""

from netgen.occ import OCCGeometry
import numpy as np
import os
import sys

BREP_FILE = "mesh/full_bracket_study/full_bracket_geometry_final.brep"

E_MOD = 71700.0
NU = 0.33
RHO = 2810e-9
FORCE_X = -235.44
REF_NODE_COORD = (60.0, 20.0, 63.0)
MAXH_GLOBAL = 6.0

TOE_FACE_IDX = 2
R3_FACE_IDX = 11
HOLE_FACE_IDX = 14

EXPECTED_AREAS = {TOE_FACE_IDX: 31.4159, R3_FACE_IDX: 188.4956, HOLE_FACE_IDX: 62.8319}

def build_mesh(level_name, hole_h, r3_h, toe_h, inp_out):
    print(f"\n=== {level_name}: hole={hole_h}, R3={r3_h}, toe={toe_h} ===")

    geo = OCCGeometry(BREP_FILE)
    shape = geo.shape
    faces = list(shape.faces)

    for idx, expected_area in EXPECTED_AREAS.items():
        area = faces[idx].mass
        ok = abs(area - expected_area) < 0.1
        print(f"  face[{idx}]: area={area:.4f} (expected {expected_area}) [{'OK' if ok else 'MISMATCH -- STOP'}]")
        if not ok:
            print("ERROR: face index/area mismatch. Aborting.")
            sys.exit(1)

    geo.SetFaceMeshsize(HOLE_FACE_IDX, hole_h)
    geo.SetFaceMeshsize(R3_FACE_IDX, r3_h)
    geo.SetFaceMeshsize(TOE_FACE_IDX, toe_h)

    print("Generating mesh...")
    mesh = geo.GenerateMesh(maxh=MAXH_GLOBAL)
    print("Mesh generated.")

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
    tets = tets[:, [0, 1, 3, 2]]  # Netgen->CalculiX orientation fix, confirmed
    print(f"Points: {n_pts}, Linear tets: {len(tets)}")

    def signed_volume(p0, p1, p2, p3):
        return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

    vols = np.array([
        signed_volume(coords[t[0]], coords[t[1]], coords[t[2]], coords[t[3]])
        for t in tets
    ])
    n_bad = (vols <= 0).sum()
    print(f"Linear volume check: min={vols.min():.6f}, mean={vols.mean():.6f}, bad={n_bad}")
    if n_bad > 0:
        print("ERROR: non-positive volume elements. Aborting.")
        sys.exit(1)

    # Sanity: verify actual mean edge length near toe reflects the target
    toe_pt = np.array([38.0, 20.0, 2.0])
    edge_lens = []
    for t in tets:
        centroid = np.mean([coords[n] for n in t], axis=0)
        if np.linalg.norm(centroid - toe_pt) < 3.0:
            c = [coords[n] for n in t]
            for i in range(4):
                for j in range(i + 1, 4):
                    edge_lens.append(np.linalg.norm(c[i] - c[j]))
    edge_lens = np.array(edge_lens)
    print(f"Near-toe (< 3mm) edges: {len(edge_lens)}, mean length={edge_lens.mean():.4f}mm "
          f"(target size: {toe_h}mm)")

    print("Promoting to C3D10 (straight-sided, deterministic midside nodes)...")
    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    midside_cache = {}
    next_node_tag = n_pts + 1
    extra_coords = {}
    quad_elems = []
    for t in tets:
        full_conn = list(t)
        for (a, b) in edge_pairs:
            na, nb = t[a], t[b]
            key = (min(na, nb), max(na, nb))
            if key not in midside_cache:
                mid = (coords[na] + coords[nb]) / 2.0
                extra_coords[next_node_tag] = mid
                midside_cache[key] = next_node_tag
                next_node_tag += 1
            full_conn.append(midside_cache[key])
        quad_elems.append(full_conn)

    total_nodes = n_pts + len(midside_cache)
    print(f"Total nodes: {total_nodes}, Total elements: {len(quad_elems)}")

    all_node_tags = list(range(1, n_pts + 1)) + list(extra_coords.keys())
    def get_coord(tag):
        return coords[tag] if tag <= n_pts else extra_coords[tag]

    root_nodes = [t for t in all_node_tags if abs(get_coord(t)[0] - 0.0) < 1e-3]
    tip_nodes = [t for t in all_node_tags if abs(get_coord(t)[2] - 60.0) < 1e-3]
    print(f"Root nodes: {len(root_nodes)}, Tip nodes: {len(tip_nodes)}")
    if len(root_nodes) == 0 or len(tip_nodes) == 0:
        print("ERROR: root/tip node set empty.")
        sys.exit(1)

    est_equations = 3 * total_nodes
    ceiling = 1_180_000
    print(f"Equation estimate: {est_equations} ({100*est_equations/ceiling:.1f}% of ceiling)")

    os.makedirs(os.path.dirname(inp_out), exist_ok=True)
    ref_tag = total_nodes + 1
    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Stage 5 -- Netgen mesh {level_name}, straight C3D10, face-based sizing\n")
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
            f.write(", ".join(str(n) for n in root_nodes[i:i+10]) + "\n")
        f.write("*NSET, NSET=NTIP\n")
        for i in range(0, len(tip_nodes), 10):
            f.write(", ".join(str(n) for n in tip_nodes[i:i+10]) + "\n")

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

    print(f"Wrote: {inp_out}")
    return total_nodes, len(quad_elems)

if __name__ == "__main__":
    import sys
    level = sys.argv[1] if len(sys.argv) > 1 else "1"

    if level == "1":
        build_mesh("Level 1", hole_h=0.5, r3_h=0.5, toe_h=0.25,
                   inp_out="simulation/stage5_netgen_L1_solve.inp")
    elif level == "2":
        build_mesh("Level 2", hole_h=0.35, r3_h=0.35, toe_h=0.15,
                   inp_out="simulation/stage5_netgen_L2_solve.inp")
    elif level == "3":
        build_mesh("Level 3", hole_h=0.25, r3_h=0.25, toe_h=0.10,
                   inp_out="simulation/stage5_netgen_L3_solve.inp")
    else:
        print(f"Unknown level '{level}', use 1, 2, or 3")
