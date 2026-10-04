"""
Stage 3 panel-verify: mesh the curved root-bay panel (shell, S6),
apply real structural BCs, and a PLACEHOLDER unit compressive load at
the z=200mm (interior rib) edge, for a *BUCKLE run.

BCs (tied to the real structure, not idealized textbook conditions):
  x=0   (spar edge):          UX=UY=UZ=0 (fully restrained -- spar is
                               much stiffer than the skin)
  z=0   (root rib edge):      UX=UY=UZ=0 (matches root fixity elsewhere)
  z=200 (interior rib edge):  UX=UY=0, UZ FREE (load applied here)
  trailing edge:               free (no trailing-edge structure in this
                               project's spec)

Load: PLACEHOLDER unit total force in -Z at the z=200 edge, split
EQUALLY among that edge's nodes (not consistent Simpson weighting --
the mesh here is Netgen's own triangulation, not a structured grid, so
exact boundary-segment lengths aren't readily available; acceptable
for a placeholder that gets replaced by the real derived load).

S6 promotion: edge midside nodes placed at the LINEAR midpoint between
corner nodes -- a small geometric approximation on curved geometry
(the true midpoint of a curved edge sits slightly off this), analogous
to the straight-sided-tet approximation used for C3D10 elsewhere in
this portfolio, but for a curved shell instead of a volume.

Usage: /usr/bin/python3 scripts/mesh_stage3_panel_verify.py <maxh_mm>
"""
import os
import sys
import numpy as np
from netgen.occ import OCCGeometry

STEP_FILE = "mesh/stage3_panel_verify/panel_verify.step"
THICKNESS = 1.0
E_MOD, NU, RHO = 73100.0, 0.33, 2780e-12
PLACEHOLDER_LOAD_N = 1000.0  # total, applied -Z at the z=200 edge
TOL = 1e-2


def main():
    if len(sys.argv) != 2:
        print("Usage: mesh_stage3_panel_verify.py <maxh_mm>"); sys.exit(1)
    maxh = float(sys.argv[1])
    hlabel = f"{maxh:.1f}".replace(".", "p")
    tag = f"stage3_panel_verify_h{hlabel}"
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    print("=" * 70)
    print(f"STAGE 3 PANEL VERIFY MESH: maxh={maxh}mm")
    print("=" * 70)

    geo = OCCGeometry(STEP_FILE)
    if len(list(geo.shape.faces)) != 1 or len(list(geo.shape.solids)) != 0:
        print("ERROR: expected 1 face, 0 solids."); sys.exit(1)

    mesh = geo.GenerateMesh(maxh=maxh)
    points = mesh.Points()
    n_pts = len(points)
    coords = np.zeros((n_pts + 1, 3))
    for i, p in enumerate(points, start=1):
        coords[i] = [p.p[0], p.p[1], p.p[2]]

    tris = []
    for el in mesh.Elements2D():
        v = [pt.nr for pt in el.vertices]
        if len(v) != 3:
            print(f"ERROR: non-triangular element found ({len(v)} nodes)."); sys.exit(1)
        tris.append(v)
    tris = np.array(tris)
    print(f"Linear triangles: {len(tris)}, nodes: {n_pts}")

    edge_pairs = [(0, 1), (1, 2), (0, 2)]
    midside_cache, extra_coords, s6_elems = {}, {}, []
    next_tag = n_pts + 1
    for t in tris:
        conn = list(t)
        for a, b in edge_pairs:
            na, nb = t[a], t[b]
            key = (min(na, nb), max(na, nb))
            if key not in midside_cache:
                extra_coords[next_tag] = (coords[na] + coords[nb]) / 2.0
                midside_cache[key] = next_tag
                next_tag += 1
            conn.append(midside_cache[key])
        s6_elems.append(conn)

    total_nodes = n_pts + len(midside_cache)
    all_coords = np.vstack([coords, np.array([extra_coords[k] for k in sorted(extra_coords)])])
    print(f"S6 nodes: {total_nodes}, S6 elements: {len(s6_elems)}")

    x_all, z_all = all_coords[1:, 0], all_coords[1:, 2]

    n_x0 = (np.nonzero(np.abs(x_all - 0.0) < TOL)[0] + 1).tolist()
    n_z0 = (np.nonzero(np.abs(z_all - 0.0) < TOL)[0] + 1).tolist()
    n_z200 = (np.nonzero(np.abs(z_all - 200.0) < TOL)[0] + 1).tolist()
    print(f"x=0 (spar) nodes: {len(n_x0)}")
    print(f"z=0 (root rib) nodes: {len(n_z0)}")
    print(f"z=200 (interior rib, loaded) nodes: {len(n_z200)}")
    for name, nset in [("x=0", n_x0), ("z=0", n_z0), ("z=200", n_z200)]:
        if not nset:
            print(f"ERROR: {name} node set empty."); sys.exit(1)

    # Real derived load (scripts/derive_stage3_panel_load.py): the
    # spar's verified bending moment M(z) applied at the panel's own
    # surface height y(x,z), NOT a placeholder anymore.
    M_Z200 = 5487.3885  # N*mm, from derive_stage3_panel_load.py at z=200mm
    I_Z200 = 7657.3140   # mm^4, same
    edge_len_estimate = None  # computed below from node spread

    def sigma_zz(x, y):
        return -M_Z200 * y / I_Z200

    xs_edge = all_coords[n_z200, 0]
    ys_edge = all_coords[n_z200, 1]
    order = np.argsort(xs_edge)
    n_z200_sorted = [n_z200[i] for i in order]
    xs_sorted = xs_edge[order]
    # tributary width per node: half-distance to each neighbor, endpoints
    # get half the distance to their single neighbor (trapezoidal rule)
    trib = np.zeros(len(xs_sorted))
    trib[1:-1] = (xs_sorted[2:] - xs_sorted[:-2]) / 2.0
    trib[0] = (xs_sorted[1] - xs_sorted[0]) / 2.0
    trib[-1] = (xs_sorted[-1] - xs_sorted[-2]) / 2.0

    node_force = {}
    for nd, x, y, w in zip(n_z200_sorted, xs_sorted, ys_edge[order], trib):
        node_force[nd] = sigma_zz(x, y) * THICKNESS * w  # N, negative = compressive
    total_applied = sum(node_force.values())
    print(f"Derived load (real, from Stage 2 M(z)/I(z)): "
          f"{total_applied:.4f} N total on z=200 edge "
          f"(reference integral: -1387.9948 N)")
    if abs((total_applied - (-1387.9948)) / 1387.9948) > 0.05:
        print("WARNING: nodal-sum total differs from the reference integral "
              "by >5% -- check tributary-width estimate before trusting this run.")

    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Stage 3 panel-verify (curved root bay), maxh={maxh}mm, S6\n")
        f.write("** PLACEHOLDER load (equal split), real BCs tied to structure\n")

        f.write("*NODE\n")
        for tg in range(1, total_nodes + 1):
            c = all_coords[tg]
            f.write(f"{tg}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")

        f.write("*ELEMENT, TYPE=S6, ELSET=EALL\n")
        for i, conn in enumerate(s6_elems, start=1):
            f.write(f"{i}, " + ", ".join(str(v) for v in conn) + "\n")

        for name, nset in [("NX0", n_x0), ("NZ0", n_z0), ("NZ200", n_z200)]:
            f.write(f"*NSET, NSET={name}\n")
            for i in range(0, len(nset), 10):
                f.write(", ".join(str(v) for v in nset[i:i+10]) + "\n")

        f.write("*MATERIAL, NAME=AL2024T3\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SHELL SECTION, ELSET=EALL, MATERIAL=AL2024T3\n")
        f.write(f"{THICKNESS}\n")

        f.write("*BOUNDARY\n")
        f.write("NX0, 1, 3, 0.0\n")
        f.write("NZ0, 1, 3, 0.0\n")
        f.write("NZ200, 1, 2, 0.0\n")

        f.write("*STEP\n*BUCKLE\n4\n")
        f.write("*CLOAD\n")
        for nd, F in node_force.items():
            f.write(f"{nd}, 3, {F:.8f}\n")
        f.write("*NODE FILE\nU\n*END STEP\n")

    print(f"\nWrote: {inp_out}")


if __name__ == "__main__":
    main()

