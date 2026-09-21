"""
Project 02, Configuration A -- full-assembly mesh generation and CalculiX
input writer. Mesh + .inp ONLY: does not launch the solver.

Usage:  /usr/bin/python3 scripts/mesh_config_a.py <maxh_mm>

Geometry: mesh/config_a_study/config_a_clean.step (single solid, 18 faces),
produced and validated by scripts/build_config_a_base_frame.py.

Load/BC (locked Project 02 assumption): bottom face z=0 fully fixed;
uniform pressure 0.1643 MPa on the top face z=100 (compressive). This is
the equal-strain / equal-stress load-sharing assumption applied over the
whole top face (sum of all 8 members' top faces = 953 mm^2), giving a
total applied force of 0.1643*953 = 156.58 N, consistent with
1.33 kg x 12 g x 9.81 = 156.6 N.

Netgen straight-sided C3D10 pipeline, same as Sub-Stages 2a/2b.
"""

import math
import os
import sys
import numpy as np
from netgen.occ import OCCGeometry

STEP_FILE = "mesh/config_a_study/config_a_clean.step"

L = 100.0                # z-height, mm
E_MOD = 68900.0          # MPa
NU = 0.33
RHO = 2700e-9            # t/mm^3
PRESSURE = 0.1643        # MPa

EXPECTED_VOLUME = 95300.0     # mm^3 (validated geometry)
EXPECTED_TOP_AREA = 953.0     # mm^2 (validated section area)
EXPECTED_FORCE = PRESSURE * EXPECTED_TOP_AREA

TOL = 1e-3
REL_TOL = 1e-6

EQ_CEILING = 1_180_000        # documented direct-solver equation ceiling
ELEM_CEILING = 230_000        # empirical WSL hardware ceiling (Project 01)


def main():
    if len(sys.argv) != 2:
        print("Usage: mesh_config_a.py <maxh_mm>")
        sys.exit(1)
    maxh = float(sys.argv[1])
    label = f"{maxh:.1f}".replace(".", "p")
    inp_out = f"simulation/config_a_h{label}.inp"
    os.makedirs("simulation", exist_ok=True)

    print("=" * 70)
    print(f"CONFIG A FULL ASSEMBLY -- mesh generation, maxh = {maxh} mm")
    print("=" * 70)

    geo = OCCGeometry(STEP_FILE)
    n_faces = len(list(geo.shape.faces))
    print(f"Geometry: {STEP_FILE}, faces = {n_faces}")
    if n_faces != 18:
        print("ERROR: expected 18 faces (clean single solid). Aborting.")
        sys.exit(1)

    mesh = geo.GenerateMesh(maxh=maxh)
    print("Mesh generated successfully.")

    points = mesh.Points()
    n_pts = len(points)
    coords = np.zeros((n_pts + 1, 3))          # row 0 is an unused dummy
    for i, p in enumerate(points, start=1):
        pnt = p.p
        coords[i] = [pnt[0], pnt[1], pnt[2]]

    tets = []
    for el in mesh.Elements3D():
        pts = [v.nr for v in el.vertices]
        if len(pts) != 4:
            continue
        tets.append(pts)
    tets = np.array(tets)
    tets = tets[:, [0, 1, 3, 2]]               # orientation for CalculiX

    # ---- vectorized volume + quality ------------------------------------
    P = coords[tets]                            # (n, 4, 3)
    v1 = P[:, 1] - P[:, 0]
    v2 = P[:, 2] - P[:, 0]
    v3 = P[:, 3] - P[:, 0]
    vols = np.einsum("ij,ij->i", np.cross(v1, v2), v3) / 6.0
    n_bad = int((vols <= 0).sum())
    total_vol = vols.sum()
    print(f"\nLinear tets: {len(tets)}, min vol = {vols.min():.6f}, bad = {n_bad}")
    if n_bad > 0:
        print("ERROR: non-positive volume elements. Aborting.")
        sys.exit(1)

    vol_err = (total_vol - EXPECTED_VOLUME) / EXPECTED_VOLUME
    print(f"Mesh volume: {total_vol:.4f} mm^3 vs {EXPECTED_VOLUME:.4f} "
          f"(rel err {vol_err:+.2e})")

    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    sumsq = np.zeros(len(tets))
    for a, b in edge_pairs:
        d = P[:, a] - P[:, b]
        sumsq += np.einsum("ij,ij->i", d, d)
    quality = 12.0 * (3.0 * vols) ** (2.0 / 3.0) / sumsq
    print(f"Tet quality q=12(3V)^(2/3)/sum(l^2) (1 = regular): "
          f"min={quality.min():.4f}, mean={quality.mean():.4f}, "
          f"n(q<0.1)={int((quality < 0.1).sum())}, "
          f"n(q<0.05)={int((quality < 0.05).sum())}")

    # ---- promote to straight-sided C3D10 ---------------------------------
    print("\nPromoting to C3D10 (straight-sided)...")
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
    n_elem = len(quad_elems)
    all_coords = np.vstack(
        [coords, np.array([extra_coords[k] for k in sorted(extra_coords)])]
    )
    assert all_coords.shape[0] == total_nodes + 1

    # ---- node sets (slice [1:] to exclude the dummy row 0) -----------------
    z_all = all_coords[1:, 2]
    bottom_nodes = (np.nonzero(np.abs(z_all - 0.0) < TOL)[0] + 1).tolist()
    top_nodes = (np.nonzero(np.abs(z_all - L) < TOL)[0] + 1).tolist()
    print(f"Bottom (z=0) nodes: {len(bottom_nodes)}, "
          f"Top (z={L}) nodes: {len(top_nodes)}")
    if not bottom_nodes or not top_nodes:
        print("ERROR: bottom/top node set empty.")
        sys.exit(1)

    # ---- top-face element faces, independent area check --------------------
    tet_faces = [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)]  # P1..P4
    dload_faces = []
    top_area = 0.0
    for fnum, idx in enumerate(tet_faces, start=1):
        idx = list(idx)
        mask = np.all(np.abs(P[:, idx, 2] - L) < TOL, axis=1)
        sel = np.nonzero(mask)[0]
        for e in sel:
            dload_faces.append((int(e) + 1, fnum))
        if len(sel):
            a0 = P[sel][:, idx[0]]
            a1 = P[sel][:, idx[1]]
            a2 = P[sel][:, idx[2]]
            top_area += 0.5 * np.linalg.norm(np.cross(a1 - a0, a2 - a0), axis=1).sum()
    area_err = (top_area - EXPECTED_TOP_AREA) / EXPECTED_TOP_AREA
    applied_force = PRESSURE * top_area
    print(f"Top-face element faces: {len(dload_faces)}, "
          f"area = {top_area:.4f} mm^2 vs {EXPECTED_TOP_AREA} "
          f"(rel err {area_err:+.2e})")
    print(f"Applied force = {applied_force:.4f} N "
          f"(expected {EXPECTED_FORCE:.4f} N)")

    # ---- resource report / gate --------------------------------------------
    est_eq = 3 * total_nodes
    print("\n" + "-" * 70)
    print(f"C3D10 nodes: {total_nodes}")
    print(f"C3D10 elements: {n_elem}  "
          f"({100 * n_elem / ELEM_CEILING:.1f}% of {ELEM_CEILING} empirical element ceiling)")
    print(f"Equations: {est_eq}  "
          f"({100 * est_eq / EQ_CEILING:.1f}% of {EQ_CEILING} equation ceiling)")
    print("-" * 70)

    gate_ok = (abs(vol_err) < REL_TOL and abs(area_err) < REL_TOL
               and len(dload_faces) > 0)
    print(f"MESH CHECKS: {'PASS' if gate_ok else 'FAIL'} "
          f"(positive volumes, volume match, top-area match)")
    if not gate_ok:
        print("ERROR: mesh checks failed. .inp NOT written.")
        sys.exit(1)

    # ---- write .inp ---------------------------------------------------------
    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Project 02, Config A full assembly, maxh={maxh} mm, axial compression\n")
        f.write("** Load assumption: uniform pressure 0.1643 MPa on top face (z=100),\n")
        f.write("** equal-stress load sharing across rails + panels; total = 0.1643*953 N.\n")
        f.write("** Bottom face (z=0) fully fixed (rigid deployer contact idealization).\n")

        f.write("*NODE\n")
        for tag in range(1, total_nodes + 1):
            c = all_coords[tag]
            f.write(f"{tag}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=EALL\n")
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
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=AL6061T6\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNBOTTOM, 1, 3, 0.0\n")
        f.write("*DSLOAD\n")
        for eid, face_num in dload_faces:
            f.write(f"{eid}, P{face_num}, {PRESSURE}\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")

    print(f"Wrote: {inp_out}")
    print("NOTE: solver NOT launched (mesh + .inp only).")


if __name__ == "__main__":
    main()
