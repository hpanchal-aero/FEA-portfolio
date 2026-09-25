"""
Project 02, Configuration B -- mesh generation and CalculiX input writer for
a pocketed-panel assembly. Mesh + .inp ONLY: does not launch the solver.

Usage:  /usr/bin/python3 scripts/mesh_config_b.py <shape> <level_pct> <maxh_mm>

Geometry: mesh/config_b_study/config_b_assembly_<shape>_<level>pct.step
(single solid, validated by scripts/build_config_b_assembly.py).

Load/BC: identical assumption to Config A (equal-strain load sharing,
total 156.5779 N, bottom face z=0 fixed), but applied via the *EQUATION
tied-top-displacement mechanism verified on Config A in this session
(scripts/make_equation_load_case.py), rather than uniform *DSLOAD pressure --
required here because panels are no longer prismatic, so a uniform pressure
would silently assume equal stress rather than let equal strain emerge.
A reference node at (0,0,150) carries the total force; every NTOP node has
Uz tied to it; the reference node's Ux,Uy are fixed to avoid a singular system.

Netgen straight-sided C3D10 pipeline, same as Config A / Sub-Stages 2a/2b.
"""

import os
import sys
import numpy as np
from netgen.occ import OCCGeometry

L = 100.0
E_MOD = 68900.0
NU = 0.33
RHO = 2700e-9
FORCE_TOTAL = 156.5779   # N, same total as Config A (0.1643 MPa x 953 mm^2)
REF_XYZ = (0.0, 0.0, 150.0)

TOL = 1e-3
REL_TOL_VOL = 1e-3   # looser than Config A: pocket geometry volume is not
                      # independently re-derived here; build script already
                      # validated the source geometry's volume exactly.

EQ_CEILING = 1_180_000
ELEM_CEILING = 230_000


def main():
    if len(sys.argv) == 4 and sys.argv[1] != "--random":
        shape, level, maxh = sys.argv[1], sys.argv[2], float(sys.argv[3])
        step_file = f"mesh/config_b_study/config_b_assembly_{shape}_{level}pct.step"
        hlabel = f"{maxh:.1f}".replace(".", "p")
        tag = f"config_b_{shape}_{level}pct_h{hlabel}"
        label_for_print = f"shape={shape} level={level}% maxh={maxh}mm"
    elif len(sys.argv) == 4 and sys.argv[1] == "--random":
        rand_tag, maxh = sys.argv[2], float(sys.argv[3])
        step_file = f"mesh/config_b_study/config_b_random_assembly_{rand_tag}.step"
        hlabel = f"{maxh:.1f}".replace(".", "p")
        tag = f"config_b_random_{rand_tag}_h{hlabel}"
        label_for_print = f"random_tag={rand_tag} maxh={maxh}mm"
    else:
        print("Usage: mesh_config_b.py <shape> <level_pct> <maxh_mm>")
        print("   or: mesh_config_b.py --random <tag> <maxh_mm>")
        sys.exit(1)
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    print("=" * 70)
    print(f"CONFIG B: {label_for_print}")
    print("=" * 70)

    if not os.path.isfile(step_file):
        print(f"ERROR: {step_file} not found."); sys.exit(1)
    geo = OCCGeometry(step_file)
    n_faces = len(list(geo.shape.faces))
    n_solids = len(list(geo.shape.solids))
    print(f"Geometry: {step_file}, solids={n_solids}, faces={n_faces}")
    if n_solids != 1:
        print("ERROR: expected 1 solid."); sys.exit(1)

    mesh = geo.GenerateMesh(maxh=maxh)
    print("Mesh generated successfully.")

    points = mesh.Points()
    n_pts = len(points)
    coords = np.zeros((n_pts + 1, 3))
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
    tets = tets[:, [0, 1, 3, 2]]

    P = coords[tets]
    v1, v2, v3 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]
    vols = np.einsum("ij,ij->i", np.cross(v1, v2), v3) / 6.0
    n_bad = int((vols <= 0).sum())
    total_vol = vols.sum()
    print(f"Linear tets: {len(tets)}, min vol={vols.min():.6f}, bad={n_bad}, "
          f"mesh volume={total_vol:.4f} mm^3")
    if n_bad > 0:
        print("ERROR: non-positive volume elements. Aborting."); sys.exit(1)

    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    midside_cache, extra_coords, quad_elems = {}, {}, []
    next_node_tag = n_pts + 1
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

    z_all = all_coords[1:, 2]
    bottom_nodes = (np.nonzero(np.abs(z_all - 0.0) < TOL)[0] + 1).tolist()
    top_nodes = (np.nonzero(np.abs(z_all - L) < TOL)[0] + 1).tolist()
    print(f"Bottom (z=0) nodes: {len(bottom_nodes)}, Top (z={L}) nodes: {len(top_nodes)}")
    if not bottom_nodes or not top_nodes:
        print("ERROR: bottom/top node set empty."); sys.exit(1)

    ref_id = total_nodes + 1
    est_eq = 3 * total_nodes + len(top_nodes)  # +1 eq per tied node (rough; *EQUATION
                                                # adds constraint eqs, not extra DOFs)
    print(f"\nC3D10 nodes: {total_nodes} (+1 reference node)")
    print(f"C3D10 elements: {n_elem}  ({100*n_elem/ELEM_CEILING:.1f}% of {ELEM_CEILING} element ceiling)")
    print(f"Equations (approx, incl. ties): {est_eq}  ({100*est_eq/EQ_CEILING:.1f}% of {EQ_CEILING} ceiling)")

    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Project 02, Config B, {label_for_print}\n")
        f.write("** Load: total 156.5779 N via *EQUATION-tied top-face Uz to a reference\n")
        f.write("** node (equal-strain load sharing emerges from stiffness, not imposed\n")
        f.write("** as uniform pressure -- required for non-prismatic pocketed panels).\n")
        f.write("** Bottom face (z=0) fully fixed.\n")

        f.write("*NODE\n")
        for tg in range(1, total_nodes + 1):
            c = all_coords[tg]
            f.write(f"{tg}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")
        f.write(f"{ref_id}, {REF_XYZ[0]:.10f}, {REF_XYZ[1]:.10f}, {REF_XYZ[2]:.10f}\n")

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

        f.write("*EQUATION\n")
        for k in top_nodes:
            f.write(f"2\n{k}, 3, 1., {ref_id}, 3, -1.\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\n")
        f.write(f"{ref_id}, 1, 2, 0.0\n")
        f.write("NBOTTOM, 1, 3, 0.0\n")
        f.write("*CLOAD\n")
        f.write(f"{ref_id}, 3, {-FORCE_TOTAL}\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")

    print(f"\nWrote: {inp_out}")
    print("NOTE: solver NOT launched (mesh + .inp only).")


if __name__ == "__main__":
    main()
