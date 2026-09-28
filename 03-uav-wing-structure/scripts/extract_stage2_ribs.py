"""
Project 03, Stage 2 -- extract CalculiX results for the spar+ribs
assembly.

Verification target: numerical tapered Euler-Bernoulli beam theory for
the spar with the ribs treated as point WEIGHTS (zero bending
stiffness). Ribs are 1.5mm plates normal to the span axis, so they add
essentially no spanwise bending stiffness; their effect on tip
deflection in this load case comes from added self-weight.

Tip metric: deflection ON THE SPAR AXIS (|x|<4.5mm). The ribs' off-axis
weight twists the low-torsion open I-section (Uy is linear in x across
the tip), and the tip-band mean over all nodes is dominated by rib
nodes and depends on how nodes fall across the chord.

Equilibrium: CalculiX's reported reaction (FORC) at constrained nodes
does NOT include external load applied directly on those nodes
(established empirically on four runs: with that load added back the
residual is <1e-4 N). The corrected expectation therefore adds the
gravity share (consistent nodal load for straight-sided C3D10:
-V/20 per corner node, +V/5 per midside node) and any CLOAD on the
fixed nodes. Both the naive and corrected errors are reported.

Usage: python3 scripts/extract_stage2_ribs.py <job_path_without_ext>
"""
import json
import os
import sys
import numpy as np

SPAN = 1000.0
H_ROOT, H_TIP = 28.8, 14.4
WEB_T, CAP_W, CAP_T = 2.0, 8.0, 3.0
RIB_T = 1.5
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
TOTAL_LIFT_N = 24.5
TOL = 1e-3

ASSEMBLY_VOLUME_MM3 = 117951.9990  # geometry (build_stage2_ribs.py)

# Rib volumes and overlap with the spar (boolean intersection); overlaps
# sum to 594.0000 so 79200 + 39345.999 - 594.000 = 117951.999.
RIB_Z = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
RIB_VOL = [11083.3800, 8977.5378, 7093.3632, 5430.8562, 3990.0168, 2770.8450]
RIB_OVERLAP = [70.1919, 131.7600, 123.1200, 114.4800, 105.8400, 48.6081]

STAGE1_TIP_DEFLECTION_MM = 3.841  # FEA, 2.0mm mesh, spar alone

COEF = np.array([-1/20.0] * 4 + [1/5.0] * 6)


def parse_inp(path):
    nodes, elems, nroot, cload = [], [], [], []
    mode = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                if u.startswith("*NODE") and not u.startswith("*NODE FILE"):
                    mode = "node"
                elif u.startswith("*ELEMENT"):
                    mode = "elem"
                elif u.startswith("*NSET") and "NROOT" in u:
                    mode = "nroot"
                elif u.startswith("*CLOAD"):
                    mode = "cload"
                else:
                    mode = None
                continue
            p = [t.strip() for t in s.split(",") if t.strip()]
            if mode == "node":
                nodes.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                elems.append([int(v) for v in p[1:11]])
            elif mode == "nroot":
                nroot.extend(int(v) for v in p)
            elif mode == "cload":
                cload.append((int(p[0]), int(p[1]), float(p[2])))
    return np.array(nodes), np.array(elems), np.array(nroot), cload


def read_frd_block(path, name, ncomp, n_nodes):
    data = None
    in_block = False
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                in_block = (line[5:13].strip() == name)
                if in_block:
                    data = np.zeros((n_nodes + 1, ncomp))
                continue
            if in_block:
                if line.startswith(" -3"):
                    in_block = False
                elif line.startswith(" -1"):
                    node = int(line[3:13])
                    for k in range(ncomp):
                        data[node, k] = float(line[13 + 12*k:25 + 12*k])
    if data is None:
        print(f"ERROR: block {name} not found"); sys.exit(1)
    return data


def I_of_z(z):
    H = H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    return (1/12.0) * (CAP_W * H**3 - (CAP_W - WEB_T) * (H - 2*CAP_T)**3)


def A_of_z(z):
    H = H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    return (H - 2*CAP_T) * WEB_T + 2 * CAP_W * CAP_T


def beam_tip(rib_weights_on):
    """Tapered cantilever tip deflection by numerical double integration.
    Ribs (if on) enter only as point weights, zero bending stiffness."""
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = TOTAL_LIFT_N / (np.pi * SPAN / 4.0)
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN)**2, 0, None)) - RHO * A_of_z(z) * G

    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q[i] + q[i+1]) * dz
    if rib_weights_on:
        for zr, vol, ov in zip(RIB_Z, RIB_VOL, RIB_OVERLAP):
            P = (vol - ov) * RHO * G
            V[z < zr - 1e-9] -= P
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dz
    curv = M / (E_MOD * I_of_z(z))
    th = np.zeros(n); v = np.zeros(n)
    for i in range(1, n):
        th[i] = th[i-1] + 0.5 * (curv[i-1] + curv[i]) * dz
    for i in range(1, n):
        v[i] = v[i-1] + 0.5 * (th[i-1] + th[i]) * dz
    return v[-1]


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = (s[:, i] for i in range(6))
    return np.sqrt(0.5*((sxx-syy)**2+(syy-szz)**2+(szz-sxx)**2) + 3*(sxy**2+syz**2+szx**2))


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_stage2_ribs.py <job>"); sys.exit(1)
    job = sys.argv[1]
    inp, frd = job + ".inp", job + ".frd"
    xyz, elems, root, cload = parse_inp(inp)
    n = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    print(f"Job: {job}\nNodes: {n}, elements: {len(elems)}")

    U = read_frd_block(frd, "DISP", 3, n)
    RF = read_frd_block(frd, "FORC", 3, n)
    S = read_frd_block(frd, "STRESS", 6, n)
    vm = von_mises(S)
    z = coords[:, 2]

    z_min = z[1:].min()
    zr = coords[root, 2]
    print(f"Fixed set (NROOT): {len(root)} nodes, z range [{zr.min():.6f}, {zr.max():.6f}] mm "
          f"(mesh z_min={z_min:.6f}, expected ~{-RIB_T/2:.4f})")

    P = coords[elems[:, :4]]
    vol = np.abs(np.einsum("ij,ij->i",
                           np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                           P[:, 3] - P[:, 0])) / 6.0
    total_vol = vol.sum()
    weight = RHO * G * total_vol
    print(f"Mesh volume {total_vol:.4f} mm^3 (geometry {ASSEMBLY_VOLUME_MM3:.4f}, "
          f"diff {100*(total_vol - ASSEMBLY_VOLUME_MM3)/ASSEMBLY_VOLUME_MM3:+.4f}%)")

    print("\n--- Equilibrium ---")
    is_root = np.zeros(n + 1, dtype=bool)
    is_root[root] = True
    grav_share_vol = (vol * (is_root[elems] * COEF).sum(axis=1)).sum()
    F_grav_root = -RHO * G * grav_share_vol
    F_lift_total = sum(f for nd, dof, f in cload if dof == 2)
    F_lift_root = sum(f for nd, dof, f in cload if dof == 2 and is_root[nd])
    expected_naive = -(F_lift_total - weight)
    expected_corr = expected_naive + F_grav_root + F_lift_root
    r = RF[root].sum(axis=0)
    print(f"Sum RF fixed nodes: X={r[0]:.6f} Y={r[1]:.6f} Z={r[2]:.6f} N")
    print(f"Applied lift {F_lift_total:.6f} N (+Y), self-weight {weight:.6f} N (-Y)")
    print(f"External load applied directly on fixed nodes: gravity share "
          f"{F_grav_root:+.6f} N, CLOAD {F_lift_root:+.6f} N "
          f"(not included in reported RF)")
    eq_naive = 100 * (r[1] - expected_naive) / abs(expected_naive)
    eq_corr = 100 * (r[1] - expected_corr) / abs(expected_corr)
    print(f"Expected reaction, naive:            {expected_naive:+.6f} N  -> error {eq_naive:+.4f}%")
    print(f"Expected reaction, fixed-node load excluded: {expected_corr:+.6f} N  -> error {eq_corr:+.4f}%  (PRIMARY)")

    print("\n--- Tip deflection (Y) ---")
    tip_band = np.nonzero(np.abs(z[1:] - SPAN) < 2.0)[0] + 1
    axis = tip_band[np.abs(coords[tip_band, 0]) < 4.5]
    x_b, uy_b = coords[tip_band, 0], U[tip_band, 1]
    slope, intercept = np.polyfit(x_b, uy_b, 1)
    uy_axis = U[axis, 1].mean()
    print(f"Tip-band nodes: {len(tip_band)}; spar-axis nodes (|x|<4.5): {len(axis)}")
    print(f"Spar-axis FEA Uy (PRIMARY metric): {uy_axis:.6f} mm")
    print(f"Regression Uy = a + b*x: a={intercept:.6f} mm, b={slope:+.6f} mm/mm "
          f"(twist-like, from off-axis rib weight)")
    print(f"(Tip-band mean over all nodes, NOT used: {uy_b.mean():.6f} mm)")

    v_spar = beam_tip(False)
    v_ribw = beam_tip(True)
    pct_ref = 100 * (uy_axis - v_ribw) / v_ribw
    print(f"\nReference (tapered beam, ribs as point weights only): {v_ribw:.6f} mm")
    print(f"FEA spar-axis vs reference: {pct_ref:+.4f}%")
    print(f"(Reference spar alone, no ribs: {v_spar:.6f} mm; added rib weight alone "
          f"changes it by {100*(v_ribw - v_spar)/v_spar:+.2f}%)")

    delta = uy_axis - STAGE1_TIP_DEFLECTION_MM
    print(f"\nTotal change vs Stage 1 FEA ({STAGE1_TIP_DEFLECTION_MM:.4f} mm): "
          f"{delta:+.4f} mm ({100*delta/STAGE1_TIP_DEFLECTION_MM:+.2f}%) -- dominated by "
          f"added rib self-weight, NOT a stiffness gain")

    print("\n--- Peak von Mises (local corner peak; convergence NOT established) ---")
    order = np.argsort(vm[1:])[::-1][:5] + 1
    for nd in order:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.4f} MPa at ({c[0]:7.3f},{c[1]:7.3f},{c[2]:7.3f})")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_results.json"
    res = {
        "job": job, "nodes": n, "elements": int(len(elems)),
        "root_reaction_Y_N": float(r[1]),
        "equilibrium_err_naive_pct": float(eq_naive),
        "equilibrium_err_fixed_node_load_excluded_pct": float(eq_corr),
        "gravity_share_on_fixed_nodes_N": float(F_grav_root),
        "uy_spar_axis_mm": float(uy_axis),
        "uy_regression_intercept_mm": float(intercept),
        "twist_slope_mm_per_mm": float(slope),
        "reference_ribs_as_weights_mm": float(v_ribw),
        "pct_vs_reference": float(pct_ref),
        "reference_spar_alone_mm": float(v_spar),
        "stage1_fea_tip_mm": STAGE1_TIP_DEFLECTION_MM,
        "peak_vm_MPa": float(vm[order[0]]),
    }
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
