"""
Project 03, Stage 2 -- mesh, load, and .inp writer for the spar + ribs
assembly (placeholder 200mm spacing, 6 ribs).

Load: identical to Stage 1 -- self-weight (gravity, *DLOAD GRAV) + the
same 24.5N PLACEHOLDER elliptical spanwise lift applied as nodal forces
on the spar's top-cap outer face (still a direct-to-spar approximation;
still to be refined in Stage 4 once skin exists). Ribs do not carry any
of this load directly in this pass -- they only stiffen the assembly.

Boundary condition: root face fixed (cantilever) -- UNLIKE STAGE 1, this
is NOT at z=0. The root rib is centered on z=0 and extends RIB_T/2 in
both directions (see build_stage2_ribs.py), so z=0 sits inside solid
material, not on a boundary face. The true root face is at the
assembly's actual minimum z, computed from the mesh itself (not
hardcoded), so this stays correct if RIB_T changes later.

Usage: /usr/bin/python3 scripts/mesh_stage2_ribs.py <maxh_mm>
"""
import math
import os
import sys
import numpy as np
from netgen.occ import OCCGeometry

STEP_FILE = "mesh/stage2_ribs/stage2_ribs.step"
SPAN = 1000.0
H_ROOT, H_TIP = 28.8, 14.4
CAP_T = 3.0
RIB_T = 1.5

E_MOD = 73100.0     # MPa, Al 2024-T3 (to be cited)
NU = 0.33
RHO = 2780e-12      # t/mm^3 (2780 kg/m^3 x 1e-12)
G = 9810.0          # mm/s^2

TOTAL_LIFT_N = 24.5   # PLACEHOLDER -- see docstring
TOL = 1e-3

EQ_CEILING = 1_180_000
ELEM_CEILING = 230_000

EXPECTED_GEOM_VOLUME = 117951.9990  # mm^3, from build_stage2_ribs.py validation


def main():
    if len(sys.argv) != 2:
        print("Usage: mesh_stage2_ribs.py <maxh_mm>"); sys.exit(1)
    maxh = float(sys.argv[1])
    hlabel = f"{maxh:.1f}".replace(".", "p")
    tag = f"stage2_ribs_h{hlabel}"
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    print("=" * 70)
    print(f"STAGE 2 SPAR+RIBS: maxh={maxh}mm")
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

    vol_pct_diff = 100 * (total_vol - EXPECTED_GEOM_VOLUME) / EXPECTED_GEOM_VOLUME
    print(f"Mesh volume vs known geometry volume ({EXPECTED_GEOM_VOLUME:.4f} mm^3): "
          f"{vol_pct_diff:+.4f}%")
    if abs(vol_pct_diff) > 0.5:
        print("ERROR: meshed volume disagrees with known geometry volume by >0.5% -- STOP.")
        sys.exit(1)

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

    # Root face: the assembly's TRUE minimum z, computed from the mesh
    # itself -- NOT hardcoded to 0.0 as in Stage 1. Expected to land near
    # -RIB_T/2 (-0.75mm) given the root rib's centered-extrusion overhang,
    # but computed directly rather than assumed.
    z_min = z_all.min()
    root_nodes = (np.nonzero(np.abs(z_all - z_min) < TOL)[0] + 1).tolist()
    expected_z_min = -RIB_T / 2.0
    print(f"Root face detected at z_min={z_min:.6f}mm "
          f"(expected ~{expected_z_min:.4f}mm from RIB_T/2 overhang, "
          f"diff={z_min - expected_z_min:+.6f}mm)")
    print(f"Root nodes: {len(root_nodes)}")
    if not root_nodes:
        print("ERROR: root node set empty."); sys.exit(1)

    # top-cap outer-edge nodes: unchanged from Stage 1 -- spar geometry
    # itself is untouched by adding ribs, so the same H(z)-based selection
    # still isolates the spar's top cap outer face. Rib surfaces are
    # curved (airfoil profile), so they are not expected to systematically
    # coincide with y=H(z)/2 except by coincidence at their fuse lines --
    # checked via the printed count below, not assumed silently.
    def H_at(z):
        return H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    # Restricted to the spar's actual cap width in x (CAP_W=8.0mm,
    # i.e. |x|<=4.0mm, with a small margin) IN ADDITION to the y
    # proximity test -- without this, rib-surface points (which also
    # pass through y=H(z)/2 somewhere along their curved airfoil
    # outline, since ribs are taller than the spar's own flange at
    # most stations) contaminate the selection. Confirmed via direct
    # inspection: ~7-8% of nodes selected by the y-only test had
    # |x|>4.5mm (outside the spar cap), and this contamination varied
    # between mesh densities -- directly causing the large (~5%)
    # mesh-to-mesh tip-deflection discrepancy and equilibrium drift
    # observed before this fix.
    CAP_W = 8.0
    x_tol = CAP_W / 2.0 + 0.5  # mm, small margin beyond the exact cap edge
    top_face_nodes = []
    for idx in range(1, total_nodes + 1):
        x, y, z = all_coords[idx, 0], all_coords[idx, 1], all_coords[idx, 2]
        if abs(y - H_at(z) / 2.0) < 0.15 and abs(x) < x_tol:
            top_face_nodes.append(idx)
    print(f"Top-cap outer-face nodes (load application): {len(top_face_nodes)}")
    if not top_face_nodes:
        print("ERROR: no top-edge nodes found for load application."); sys.exit(1)

    # Elliptical lift q(z) = q0*sqrt(1-(z/b)^2), total = TOTAL_LIFT_N.
    # Load is assigned SLAB BY SLAB in z: each 5mm slab receives exactly
    # the analytic integral of q(z) over that slab, shared equally by the
    # load nodes inside it. The previous per-node weighting (force
    # proportional to q at each node) over-weighted regions where the mesh
    # is locally refined around ribs (2.0mm mesh: -4.3% load moment vs
    # intended, 10 load nodes on the fixed root face); slab assignment is
    # independent of node density. Nodes outside 0<=z<=SPAN (root/tip rib
    # overhang) are excluded.
    SLAB = 5.0
    q0 = TOTAL_LIFT_N / (math.pi * SPAN / 4.0)

    def cum_q(zz):
        u = min(max(zz / SPAN, 0.0), 1.0)
        return q0 * (SPAN / 2.0) * (u * math.sqrt(1.0 - u * u) + math.asin(u))

    top_face_nodes = [nd for nd in top_face_nodes if 0.0 <= all_coords[nd, 2] <= SPAN]
    z_vals = np.array([all_coords[i, 2] for i in top_face_nodes])
    n_slabs = int(round(SPAN / SLAB))
    slab_idx = np.minimum((z_vals / SLAB).astype(int), n_slabs - 1)
    counts = np.bincount(slab_idx, minlength=n_slabs)
    if (counts == 0).any():
        print(f"ERROR: {int((counts == 0).sum())} empty load slab(s) -- STOP.")
        sys.exit(1)
    nodal_forces = np.zeros(len(top_face_nodes))
    for k in range(n_slabs):
        slab_load = cum_q((k + 1) * SLAB) - cum_q(k * SLAB)
        nodal_forces[slab_idx == k] = slab_load / counts[k]
    print(f"Load nodes after slab assignment: {len(top_face_nodes)}, "
          f"total force {nodal_forces.sum():.6f} N")

    est_eq = 3 * total_nodes
    print(f"\nC3D10 nodes: {total_nodes}  ({100*n_elem/ELEM_CEILING:.1f}% of element ceiling)")
    print(f"Equations: {est_eq}  ({100*est_eq/EQ_CEILING:.1f}% of ceiling)")
    print(f"Total applied lift: {TOTAL_LIFT_N} N (PLACEHOLDER load magnitude)")

    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Project 03 Stage 2, spar+ribs (200mm placeholder spacing), maxh={maxh}mm\n")
        f.write("** PLACEHOLDER load: 24.5N total elliptical lift (5kg aircraft, n=1)\n")
        f.write("** + self-weight (gravity). Root (true assembly z_min) fixed.\n")
        f.write(f"** Root face at z={z_min:.6f}mm (root rib center at z=0, "
                f"RIB_T/2={RIB_T/2:.2f}mm overhang).\n")

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
