"""
Project 02, Sub-Stage 2b -- isolated rail verification.

8.5x8.5x100mm rail, pure axial compression, base fixed, uniform
pressure 0.1643 MPa applied at the free top face (z=100).

Analytical reference: sigma = P/A, delta = PL/(AE)
Euler buckling closed-form check (K=2, fixed-free) included.
No eigenvalue buckling FEA.

Netgen-based straight-sided C3D10 pipeline, consistent with the
toolchain locked in Project 01 Stage 5 onward.
"""

from netgen.occ import OCCGeometry, Box, Pnt
import numpy as np
import math
import os
import sys

SIDE = 8.5     # x and y, mm (rail cross-section)
L = 100.0      # z, mm (height)

E_MOD = 68900.0
NU = 0.33
RHO = 2700e-9

PRESSURE = 0.1643  # MPa
AREA = SIDE * SIDE
FORCE = PRESSURE * AREA

MESH_H = 1.5   # mm, comfortably under the 8.5mm cross-section
TOL = 1e-3

OUT_DIR = "mesh/rail_2b_study"
os.makedirs(OUT_DIR, exist_ok=True)
INP_OUT = "simulation/stage2b_rail_solve.inp"
os.makedirs(os.path.dirname(INP_OUT), exist_ok=True)


def main():
    print("=" * 70)
    print("SUB-STAGE 2b: ISOLATED RAIL")
    print("=" * 70)
    print(f"Geometry: {SIDE}x{SIDE}x{L}mm, Area={AREA}mm^2")
    print(f"Applied pressure: {PRESSURE} MPa -> Force = {FORCE:.4f} N")

    sigma_theory = PRESSURE
    delta_theory = FORCE * L / (AREA * E_MOD)
    print(f"\nAnalytical: sigma = {sigma_theory:.4f} MPa, delta = {delta_theory:.6f} mm")

    I = SIDE ** 4 / 12
    K = 2.0  # fixed-free
    P_cr = (math.pi ** 2) * E_MOD * I / (K * L) ** 2
    print(f"Euler buckling (K={K}, fixed-free): I={I:.4f} mm^4, P_cr={P_cr:.2f} N")
    print(f"Applied/critical load ratio = {FORCE / P_cr:.4f}")

    print("\nBuilding geometry...")
    box = Box(Pnt(-SIDE / 2, -SIDE / 2, 0.0), Pnt(SIDE / 2, SIDE / 2, L))
    geo = OCCGeometry(box)

    print(f"Meshing (target size {MESH_H}mm)...")
    mesh = geo.GenerateMesh(maxh=MESH_H)
    print("Mesh generated successfully.")

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
    tets = tets[:, [0, 1, 3, 2]]

    def signed_volume(p0, p1, p2, p3):
        return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

    vols = np.array([
        signed_volume(coords[t[0]], coords[t[1]], coords[t[2]], coords[t[3]])
        for t in tets
    ])
    n_bad = (vols <= 0).sum()
    print(f"Linear tets: {len(tets)}, min vol={vols.min():.6f}, bad={n_bad}")
    if n_bad > 0:
        print("ERROR: non-positive volume elements. Aborting.")
        sys.exit(1)

    print("Promoting to C3D10 (straight-sided)...")
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

    bottom_nodes = [t for t in all_node_tags if abs(get_coord(t)[2] - 0.0) < TOL]
    top_nodes = [t for t in all_node_tags if abs(get_coord(t)[2] - L) < TOL]
    print(f"Bottom (z=0) nodes: {len(bottom_nodes)}, Top (z={L}) nodes: {len(top_nodes)}")
    if len(bottom_nodes) == 0 or len(top_nodes) == 0:
        print("ERROR: bottom/top node set empty.")
        sys.exit(1)

    est_equations = 3 * total_nodes
    ceiling = 1_180_000
    print(f"Equation estimate: {est_equations} ({100 * est_equations / ceiling:.1f}% of ceiling)")

    tet_face_defs = [
        (0, 1, 2, 4, 5, 6),
        (0, 1, 3, 4, 8, 7),
        (1, 2, 3, 5, 9, 8),
        (0, 2, 3, 6, 9, 7),
    ]
    dload_faces = []
    for eid, conn in enumerate(quad_elems, start=1):
        for face_num, idxs in enumerate(tet_face_defs, start=1):
            corner_idxs = idxs[:3]
            if all(abs(get_coord(conn[ci])[2] - L) < TOL for ci in corner_idxs):
                dload_faces.append((eid, face_num))
    print(f"Top-face element faces for pressure load: {len(dload_faces)}")
    if len(dload_faces) == 0:
        print("ERROR: no top-face element faces found.")
        sys.exit(1)

    with open(INP_OUT, "w") as f:
        f.write("*HEADING\n")
        f.write("Project 02, Sub-Stage 2b -- isolated rail, axial compression\n")

        f.write("*NODE\n")
        for tag in range(1, n_pts + 1):
            c = coords[tag]
            f.write(f"{tag}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")
        for tag, c in extra_coords.items():
            f.write(f"{tag}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=ERAIL\n")
        for i, conn in enumerate(quad_elems, start=1):
            f.write(f"{i}, " + ", ".join(str(n) for n in conn) + "\n")

        f.write("*NSET, NSET=NBOTTOM\n")
        for i in range(0, len(bottom_nodes), 10):
            f.write(", ".join(str(n) for n in bottom_nodes[i:i + 10]) + "\n")
        f.write("*NSET, NSET=NTOP\n")
        for i in range(0, len(top_nodes), 10):
            f.write(", ".join(str(n) for n in top_nodes[i:i + 10]) + "\n")

        f.write("*MATERIAL, NAME=AL6061T6\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=ERAIL, MATERIAL=AL6061T6\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNBOTTOM, 1, 3, 0.0\n")
        f.write("*DSLOAD\n")
        for eid, face_num in dload_faces:
            f.write(f"{eid}, P{face_num}, {PRESSURE}\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")

    print(f"Wrote: {INP_OUT}")


if __name__ == "__main__":
    main()
