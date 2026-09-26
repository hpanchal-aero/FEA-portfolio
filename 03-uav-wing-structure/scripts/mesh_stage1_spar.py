"""
Project 03, Stage 1 -- mesh, load, and .inp writer for the isolated spar.

Load: self-weight (gravity, *DLOAD GRAV) + an elliptical spanwise line
load applied as nodal forces on the TOP cap's outer edge nodes (a direct-
to-spar approximation, since no skin exists yet to distribute the load
over -- to be refined in Stage 4 once skin exists).

PLACEHOLDER total load: assumes total aircraft mass 5 kg, n=1 (straight
and level), giving 24.5 N total half-span lift -- NOT a real design load,
explicitly flagged for revision once an actual aircraft weight/lift
requirement is provided.

Boundary condition: root face (z=0) fully fixed (cantilever).

Usage: /usr/bin/python3 scripts/mesh_stage1_spar.py <maxh_mm>
"""
import math
import os
import sys
import numpy as np
from netgen.occ import OCCGeometry

STEP_FILE = "mesh/stage1_spar/stage1_spar.step"
SPAN = 1000.0
H_ROOT, H_TIP = 28.8, 14.4
CAP_T = 3.0

E_MOD = 73100.0     # MPa, Al 2024-T3 (to be cited)
NU = 0.33
RHO = 2780e-12      # t/mm^3 (2780 kg/m^3 x 1e-12)
G = 9810.0          # mm/s^2

TOTAL_LIFT_N = 24.5   # PLACEHOLDER -- see docstring
TOL = 1e-3

EQ_CEILING = 1_180_000
ELEM_CEILING = 230_000


def main():
    if len(sys.argv) != 2:
        print("Usage: mesh_stage1_spar.py <maxh_mm>"); sys.exit(1)
    maxh = float(sys.argv[1])
    hlabel = f"{maxh:.1f}".replace(".", "p")
    tag = f"stage1_spar_h{hlabel}"
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    print("=" * 70)
    print(f"STAGE 1 SPAR: maxh={maxh}mm")
    print("=" * 70)

    geo = OCCGeometry(STEP_FILE)
    n_faces = len(list(geo.shape.faces))
    if len(list(geo.shape.solids)) != 1:
        print("ERROR: expected 1 solid."); sys.exit(1)
    print(f"Geometry: {n_faces} faces")

    mesh = geo.GenerateMesh(maxh=maxh)
    points = mesh.Points()
    n_pts = len(points)
    coords = np.zeros((n_pts + 1, 3))
    for i, p in enumerate(points, start=1):
        coords[i] = [p.p[0], p.p[1], p.p[2]]

    tets = []
    for el in mesh.Elements3D():
        pts = [v.nr for v in el.vertices]
        if len(pts) == 4:
            tets.append(pts)
    tets = np.array(tets)
    tets = tets[:, [0, 1, 3, 2]]

    P = coords[tets]
    v1, v2, v3 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]
    vols = np.einsum("ij,ij->i", np.cross(v1, v2), v3) / 6.0
    n_bad = int((vols <= 0).sum())
    total_vol = vols.sum()
    print(f"Linear tets: {len(tets)}, min vol={vols.min():.6f}, bad={n_bad}, "
          f"mesh volume={total_vol:.4f} mm^3")
    if n_bad > 0:
        print("ERROR: non-positive volume elements."); sys.exit(1)

    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    midside_cache, extra_coords, quad_elems = {}, {}, []
    next_tag = n_pts + 1
    for t in tets:
        conn = list(t)
        for a, b in edge_pairs:
            na, nb = t[a], t[b]
            key = (min(na, nb), max(na, nb))
            if key not in midside_cache:
                extra_coords[next_tag] = (coords[na] + coords[nb]) / 2.0
                midside_cache[key] = next_tag
                next_tag += 1
            conn.append(midside_cache[key])
        quad_elems.append(conn)

    total_nodes = n_pts + len(midside_cache)
    n_elem = len(quad_elems)
    all_coords = np.vstack([coords, np.array([extra_coords[k] for k in sorted(extra_coords)])])

    z_all = all_coords[1:, 2]
    y_all = all_coords[1:, 1]
    root_nodes = (np.nonzero(np.abs(z_all - 0.0) < TOL)[0] + 1).tolist()
    print(f"Root (z=0) nodes: {len(root_nodes)}")
    if not root_nodes:
        print("ERROR: root node set empty."); sys.exit(1)

    # top-cap outer-edge nodes: y close to local H(z)/2 - small tol, all z
    # H(z) tapers linearly from H_ROOT to H_TIP
    def H_at(z):
        return H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    # NOTE: this selects the entire top cap OUTER FACE (full CAP_W width x
    # full span), not a thin edge line -- a reasonable Stage-1 approximation
    # for where an eventual skin's aero pressure would act, given no skin
    # exists yet to distribute the load over properly (see Stage 4).
    top_face_nodes = []
    for idx in range(1, total_nodes + 1):
        y, z = all_coords[idx, 1], all_coords[idx, 2]
        if abs(y - H_at(z) / 2.0) < 0.15:
            top_face_nodes.append(idx)
    print(f"Top-cap outer-face nodes (load application): {len(top_face_nodes)}")
    if not top_face_nodes:
        print("ERROR: no top-edge nodes found for load application."); sys.exit(1)

    # elliptical line load q(z) = q0*sqrt(1-(z/b)^2); total = q0*pi*b/4 = TOTAL_LIFT_N
    q0 = TOTAL_LIFT_N / (math.pi * SPAN / 4.0)
    # distribute total force among top_face_nodes weighted by local q(z) and
    # an approximate tributary length (SPAN / n_nodes_per_unit_z, simplified
    # as equal split among nodes at similar z, weighted by q(z) at that z)
    z_vals = np.array([all_coords[i, 2] for i in top_face_nodes])
    q_vals = q0 * np.sqrt(np.clip(1.0 - (z_vals / SPAN) ** 2, 0.0, None))
    weights = q_vals / q_vals.sum() if q_vals.sum() > 0 else np.ones_like(q_vals) / len(q_vals)
    nodal_forces = weights * TOTAL_LIFT_N  # N, applied in -y (lift direction... using +y here as "up")

    est_eq = 3 * total_nodes
    print(f"\nC3D10 nodes: {total_nodes}  ({100*n_elem/ELEM_CEILING:.1f}% of element ceiling)")
    print(f"Equations: {est_eq}  ({100*est_eq/EQ_CEILING:.1f}% of ceiling)")
    print(f"Total applied lift: {TOTAL_LIFT_N} N (PLACEHOLDER load magnitude)")

    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Project 03 Stage 1, spar alone, maxh={maxh}mm\n")
        f.write("** PLACEHOLDER load: 24.5N total elliptical lift (5kg aircraft, n=1)\n")
        f.write("** + self-weight (gravity). Root (z=0) fixed. To be refined Stage 4.\n")

        f.write("*NODE\n")
        for tg in range(1, total_nodes + 1):
            c = all_coords[tg]
            f.write(f"{tg}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=EALL\n")
        for i, conn in enumerate(quad_elems, start=1):
            f.write(f"{i}, " + ", ".join(str(n) for n in conn) + "\n")

        f.write("*NSET, NSET=NROOT\n")
        for i in range(0, len(root_nodes), 10):
            f.write(", ".join(str(n) for n in root_nodes[i:i+10]) + "\n")

        f.write("*MATERIAL, NAME=AL2024T3\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=AL2024T3\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNROOT, 1, 6, 0.0\n")
        f.write("*DLOAD\nEALL, GRAV, 9810.0, 0.0, -1.0, 0.0\n")
        f.write("*CLOAD\n")
        for nd, force in zip(top_face_nodes, nodal_forces):
            f.write(f"{nd}, 2, {force:.8f}\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")

    print(f"\nWrote: {inp_out}")


if __name__ == "__main__":
    main()
