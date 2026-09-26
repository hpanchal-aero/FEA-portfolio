"""
Project 03, Stage 1 -- extract CalculiX results for the isolated spar
and check against a tapered Euler-Bernoulli beam reference (numerical
integration, since no closed form exists for this taper + combined
load), plus a quick Tier-1 prismatic sanity bound.

Usage: /usr/bin/python3 scripts/extract_stage1_spar.py <job_path_without_ext>
"""
import json
import os
import sys
import numpy as np

SPAN = 1000.0
H_ROOT, H_TIP = 28.8, 14.4
WEB_T, CAP_W, CAP_T = 2.0, 8.0, 3.0
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
TOTAL_LIFT_N = 24.5
TOL = 1e-3


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


def read_inp_nodes(path):
    xyz, mode = [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s: continue
            if s.startswith("*"):
                u = s.upper()
                mode = "node" if (u.startswith("*NODE") and not u.startswith("*NODE FILE")) else None
                continue
            if mode == "node":
                p = s.split(",")
                xyz.append([float(v) for v in p[1:4]])
    return np.array(xyz)


def I_of_z(z):
    H = H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    return (1/12.0) * (CAP_W * H**3 - (CAP_W - WEB_T) * (H - 2*CAP_T)**3)


def A_of_z(z):
    H = H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    return (H - 2*CAP_T) * WEB_T + 2 * CAP_W * CAP_T


def tapered_beam_tip_deflection():
    """Numerically integrate a cantilever (fixed at z=0) under elliptical
    lift + self-weight, tapered I(z), via double integration of M(z)/EI(z)
    (Macaulay/direct numerical double-integration on a fine grid)."""
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]

    q0 = TOTAL_LIFT_N / (np.pi * SPAN / 4.0)
    q_lift = q0 * np.sqrt(np.clip(1.0 - (z / SPAN)**2, 0, None))          # N/mm, up
    q_weight = -RHO * A_of_z(z) * G                                       # N/mm, down (self-weight)
    q_total = q_lift + q_weight

    # shear V(z) = integral of q from z to SPAN (free end V=0 at tip)
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q_total[i] + q_total[i+1]) * dz
    # moment M(z) = integral of V from z to SPAN (free end M=0 at tip)
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dz

    curvature = M / (E_MOD * I_of_z(z))
    # theta(z): integrate curvature from root (theta=0 at root, fixed)
    theta = np.zeros(n)
    for i in range(1, n):
        theta[i] = theta[i-1] + 0.5 * (curvature[i-1] + curvature[i]) * dz
    # deflection: integrate theta from root (v=0 at root)
    v = np.zeros(n)
    for i in range(1, n):
        v[i] = v[i-1] + 0.5 * (theta[i-1] + theta[i]) * dz

    return z, v, M, V, q_total


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = (s[:, i] for i in range(6))
    return np.sqrt(0.5*((sxx-syy)**2+(syy-szz)**2+(szz-sxx)**2) + 3*(sxy**2+syz**2+szx**2))


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_stage1_spar.py <job>"); sys.exit(1)
    job = sys.argv[1]
    inp, frd = job + ".inp", job + ".frd"
    xyz = read_inp_nodes(inp)
    n = len(xyz)
    coords = np.vstack([np.zeros((1,3)), xyz])
    print(f"Job: {job}\nNodes: {n}")

    U = read_frd_block(frd, "DISP", 3, n)
    RF = read_frd_block(frd, "FORC", 3, n)
    S = read_frd_block(frd, "STRESS", 6, n)
    vm = von_mises(S)
    z = coords[:, 2]

    root = np.nonzero(np.abs(z[1:] - 0.0) < TOL)[0] + 1
    tip_band = np.nonzero(np.abs(z[1:] - SPAN) < 2.0)[0] + 1

    print("\n--- Equilibrium ---")
    total_weight = RHO * G * np.trapz([A_of_z(zz) for zz in np.linspace(0, SPAN, 500)],
                                        np.linspace(0, SPAN, 500))
    net_applied_Y = TOTAL_LIFT_N - total_weight   # net applied force on the structure, +Y positive
    expected_reaction_Y = -net_applied_Y           # root reaction must be equal and opposite
    r = RF[root].sum(axis=0)
    print(f"Root nodes: {len(root)}")
    print(f"Sum RF root: X={r[0]:.6f} Y={r[1]:.6f} Z={r[2]:.6f} N")
    print(f"Applied lift: {TOTAL_LIFT_N:.4f} N (+Y);  self-weight: {total_weight:.4f} N (-Y)")
    print(f"Net applied force (+Y): {net_applied_Y:.4f} N;  expected root reaction: {expected_reaction_Y:.4f} N")
    eq_err = 100 * (r[1] - expected_reaction_Y) / expected_reaction_Y
    print(f"Equilibrium error: {eq_err:+.4f}%")

    print("\n--- Tip deflection (Y) ---")
    uy_tip = U[tip_band, 1]
    print(f"Tip-band nodes: {len(tip_band)}")
    print(f"FEA Uy: mean={uy_tip.mean():.6f} min={uy_tip.min():.6f} max={uy_tip.max():.6f} mm")

    z_ref, v_ref, M_ref, V_ref, q_ref = tapered_beam_tip_deflection()
    print(f"\nTier 2 (tapered numerical beam theory) tip deflection: {v_ref[-1]:.6f} mm")
    pct = 100 * (uy_tip.mean() - v_ref[-1]) / v_ref[-1]
    print(f"FEA vs Tier-2 theory: {pct:+.4f}%")

    # Tier 1 quick sanity bound (prismatic root section, UDL approx, no self-weight)
    I_root = I_of_z(0.0)
    w_udl = TOTAL_LIFT_N / SPAN
    v_tier1 = w_udl * SPAN**4 / (8 * E_MOD * I_root)
    print(f"\nTier 1 (prismatic root, UDL approx, sanity bound only): {v_tier1:.6f} mm")
    print(f"(Tier 1 is expected to UNDERESTIMATE vs FEA/Tier-2, since it uses the "
          f"stiffer root section throughout and omits self-weight)")

    print("\n--- Peak von Mises ---")
    order = np.argsort(vm[1:])[::-1][:5] + 1
    for nd in order:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.4f} MPa at ({c[0]:7.3f},{c[1]:7.3f},{c[2]:7.3f})")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_results.json"
    res = {
        "job": job, "nodes": n,
        "root_reaction_Y_N": float(r[1]), "equilibrium_err_pct": float(eq_err),
        "uy_tip_mean_mm": float(uy_tip.mean()),
        "tier2_theory_tip_mm": float(v_ref[-1]), "tier2_pct_diff": float(pct),
        "tier1_sanity_bound_mm": float(v_tier1),
        "peak_vm_MPa": float(vm[order[0]]),
    }
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
