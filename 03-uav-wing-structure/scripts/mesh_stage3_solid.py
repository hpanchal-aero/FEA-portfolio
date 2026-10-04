#!/usr/bin/env python3
"""
Project 03, Stage 3 baseline (no skin) -- mesh, load, and .inp writer for
spar v2 (Option C: spar depth = airfoil depth) alone, or spar v2 + 6 ribs.

Derived from mesh_stage2_ribs.py; the old Stage 1/2 scripts are untouched.
Changes: STEP file per mode; expected volume = values validated in
build_stage3_solid.py; top-cap outer face now sits on the airfoil envelope, so
load nodes are selected by |y - h_env(z)| < 0.15 and |x| <= CAP_W/2 (cap edges
are now geometric boundaries); hard element-ceiling stop; load-moment check.

Load: identical placeholder to Stages 1/2 -- self-weight (*DLOAD GRAV) + 24.5 N
PLACEHOLDER elliptical lift as nodal forces on the top-cap outer face, assigned
slab-by-slab (5 mm slabs, exact analytic integral per slab, split equally among
the slab's load nodes -- independent of mesh density).

Root: fixed at the mesh's own z_min (0 for spar, -RIB_T/2 for ribs).

Usage (system Python):
  /usr/bin/python3 scripts/mesh_stage3_solid.py <spar|ribs> <maxh_mm>
"""
import math
import os
import sys
import numpy as np
from netgen.occ import OCCGeometry

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_stage3_solid import h_env   # envelope half-height at x=0 (spar axis)

SPAN = 1000.0
CAP_W = 8.0
RIB_T = 1.5

E_MOD = 73100.0     # MPa, Al 2024-T3 (to be cited)
NU = 0.33
RHO = 2780e-12      # t/mm^3 (2780 kg/m^3 x 1e-12)

TOTAL_LIFT_N = 24.5   # PLACEHOLDER -- not a sourced design requirement
TOL = 1e-3

EQ_CEILING = 1_180_000
ELEM_CEILING = 230_000

MODES = {
    "spar": dict(step="mesh/stage3_solid/stage3_spar_v2.step",
                 vol=89910.4455, z_min=0.0, name="stage3_spar_v2"),
    "ribs": dict(step="mesh/stage3_solid/stage3_spar_ribs_v2.step",
                 vol=128582.1216, z_min=-RIB_T / 2.0, name="stage3_ribs_v2"),
}


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in MODES:
        print("Usage: mesh_stage3_solid.py <spar|ribs> <maxh_mm>"); sys.exit(1)
    mode = sys.argv[1]
    cfg = MODES[mode]
    maxh = float(sys.argv[2])
    hlabel = f"{maxh:.1f}".replace(".", "p")
    tag = f"{cfg['name']}_h{hlabel}"
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    print("=" * 70)
    print(f"STAGE 3 BASELINE (no skin): mode={mode}, maxh={maxh}mm")
    print("=" * 70)

    if not os.path.isfile(cfg["step"]):
        print(f"ERROR: {cfg['step']} not found -- run build_stage3_solid.py {mode}."); sys.exit(1)

    geo = OCCGeometry(cfg["step"])
    if len(list(geo.shape.solids)) != 1:
        print("ERROR: expected 1 solid."); sys.exit(1)
    print(f"Geometry: {len(list(geo.shape.faces))} faces")

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

    n_elem = len(tets)
    n_nodes_est = n_pts + int(6 * n_elem / 5.5)   # rough pre-check only
    print(f"Linear tets: {n_elem}, linear nodes: {n_pts}")
    if n_elem > ELEM_CEILING:
        print(f"ERROR: {n_elem} elements exceeds ceiling {ELEM_CEILING} -- STOP (WSL crash history)."); sys.exit(1)

    P = coords[tets]
    v1, v2, v3 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]
    vols = np.einsum("ij,ij->i", np.cross(v1, v2), v3) / 6.0
    n_bad = int((vols <= 0).sum())
    total_vol = vols.sum()
    print(f"Min tet vol={vols.min():.6f}, bad={n_bad}, mesh volume={total_vol:.4f} mm^3")
    if n_bad > 0:
        print("ERROR: non-positive volume elements."); sys.exit(1)
    vol_pct = 100 * (total_vol - cfg["vol"]) / cfg["vol"]
    print(f"Mesh volume vs validated geometry volume ({cfg['vol']:.4f}): {vol_pct:+.4f}%")
    if abs(vol_pct) > 0.5:
        print("ERROR: meshed volume disagrees with geometry by >0.5% -- STOP."); sys.exit(1)

    # promote linear tets to C3D10 (straight-sided midside nodes)
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
    all_coords = np.vstack([coords, np.array([extra_coords[k] for k in sorted(extra_coords)])])
    est_eq = 3 * total_nodes
    print(f"\nC3D10 nodes: {total_nodes}, elements: {n_elem} "
          f"({100*n_elem/ELEM_CEILING:.1f}% of element ceiling)")
    print(f"Equations: {est_eq} ({100*est_eq/EQ_CEILING:.1f}% of ceiling)")
    if est_eq > EQ_CEILING:
        print("ERROR: equation count exceeds ceiling -- STOP."); sys.exit(1)

    z_all = all_coords[1:, 2]
    z_min = z_all.min()
    root_nodes = (np.nonzero(np.abs(z_all - z_min) < TOL)[0] + 1).tolist()
    print(f"Root face at z_min={z_min:.6f} mm (expected {cfg['z_min']:.4f}, "
          f"diff={z_min - cfg['z_min']:+.6f}); root nodes: {len(root_nodes)}")
    if not root_nodes or abs(z_min - cfg["z_min"]) > TOL:
        print("ERROR: root face not where expected -- STOP."); sys.exit(1)

    # load nodes: top-cap outer face (on the airfoil envelope), |x| within cap width
    x_tol = CAP_W / 2.0 + TOL
    top_face_nodes = []
    for idx in range(1, total_nodes + 1):
        x, y, z = all_coords[idx]
        if abs(x) <= x_tol and abs(y - h_env(z)) < 0.15:
            top_face_nodes.append(idx)
    print(f"Top-cap outer-face nodes (before z filter): {len(top_face_nodes)}")
    top_face_nodes = [nd for nd in top_face_nodes if 0.0 <= all_coords[nd, 2] <= SPAN]
    if not top_face_nodes:
        print("ERROR: no load nodes found."); sys.exit(1)

    SLAB = 5.0
    q0 = TOTAL_LIFT_N / (math.pi * SPAN / 4.0)

    def cum_q(zz):
        u = min(max(zz / SPAN, 0.0), 1.0)
        return q0 * (SPAN / 2.0) * (u * math.sqrt(1.0 - u * u) + math.asin(u))

    z_vals = np.array([all_coords[i, 2] for i in top_face_nodes])
    n_slabs = int(round(SPAN / SLAB))
    slab_idx = np.minimum((z_vals / SLAB).astype(int), n_slabs - 1)
    counts = np.bincount(slab_idx, minlength=n_slabs)
    if (counts == 0).any():
        print(f"ERROR: {int((counts == 0).sum())} empty load slab(s) -- STOP."); sys.exit(1)
    nodal_forces = np.zeros(len(top_face_nodes))
    for k in range(n_slabs):
        nodal_forces[slab_idx == k] = (cum_q((k + 1) * SLAB) - cum_q(k * SLAB)) / counts[k]

    on_root = sum(1 for nd in top_face_nodes if abs(all_coords[nd, 2] - z_min) < TOL)
    M_app = float(np.sum(nodal_forces * z_vals))
    M_int = q0 * SPAN ** 2 / 3.0
    print(f"Load nodes: {len(top_face_nodes)} (on fixed root face: {on_root}); "
          f"total force {nodal_forces.sum():.6f} N (intended {TOTAL_LIFT_N})")
    print(f"Applied load moment about z=0: {M_app:.4f} N*mm vs intended {M_int:.4f} "
          f"({100*(M_app-M_int)/M_int:+.4f}%)")
    if abs(M_app - M_int) / M_int > 0.005:
        print("ERROR: load moment off by >0.5% -- STOP."); sys.exit(1)

    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Project 03 Stage 3 baseline (no skin), {mode}, spar v2, maxh={maxh}mm\n")
        f.write("** PLACEHOLDER load: 24.5N total elliptical lift (5kg aircraft, n=1)\n")
        f.write("** + self-weight (gravity). Root (true z_min) fixed.\n")
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
