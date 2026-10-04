#!/usr/bin/env python3
"""
Project 03, Stage 3 baseline (no skin), spar v2 ALONE -- extract CalculiX results
and verify against a tapered Euler-Bernoulli beam reference using the ACTUAL
spar v2 section (web + caps trimmed by the airfoil outline), plus:
  * equilibrium with the CalculiX constrained-node reaction convention
    (Addendum A): external load applied directly on fixed nodes is NOT in the
    reported reaction, so the corrected expectation removes CLOAD on root nodes
    and the consistent gravity share on root nodes;
  * flange-stress check sigma_zz vs -M(z)*y/I(z) (Addendum B).

Usage (system Python):
  /usr/bin/python3 scripts/extract_stage3_solid.py <job_path_without_ext>
"""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_stage3_solid as b   # h_env, spar_area, upper_xy, interp_upper

SPAN = 1000.0
WEB_T, CAP_W, CAP_T = b.WEB_T, b.CAP_W, b.CAP_T
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
TOTAL_LIFT_N = 24.5
TOL = 1e-3


def trapz(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


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


def read_inp(path):
    """nodes (xyz), C3D10 elements (list of 10 node ids), CLOAD list, root set."""
    xyz, elems, cloads, root = [], [], [], []
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
                    mode = "root"
                elif u.startswith("*CLOAD"):
                    mode = "cload"
                else:
                    mode = None
                continue
            p = [t.strip() for t in s.split(",") if t.strip()]
            if mode == "node":
                xyz.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                elems.append([int(v) for v in p[1:11]])
            elif mode == "root":
                root += [int(v) for v in p]
            elif mode == "cload":
                cloads.append((int(p[0]), int(p[1]), float(p[2])))
    return np.array(xyz), np.array(elems), cloads, root


# ---- spar v2 section properties (tables on a 5 mm grid, linear interp) ----
ZTAB = np.linspace(0.0, SPAN, 201)


def section_props(z):
    xs, ys = b.upper_xy(z)
    h = b.interp_upper(xs, ys, 0.0)
    y_in = h - CAP_T
    xg = np.linspace(-CAP_W / 2.0, CAP_W / 2.0, 801)
    top = np.array([b.interp_upper(xs, ys, x) for x in xg])
    A = WEB_T * 2.0 * y_in + 2.0 * trapz(top - y_in, xg)
    I = WEB_T * 2.0 * y_in**3 / 3.0 + 2.0 * trapz((top**3 - y_in**3) / 3.0, xg)
    return A, I


_tab = np.array([section_props(z) for z in ZTAB])
A_TAB, I_TAB = _tab[:, 0], _tab[:, 1]


def A_of_z(z):
    return np.interp(z, ZTAB, A_TAB)


def I_of_z(z):
    return np.interp(z, ZTAB, I_TAB)


def tapered_beam_reference():
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = TOTAL_LIFT_N / (np.pi * SPAN / 4.0)
    q_total = q0 * np.sqrt(np.clip(1.0 - (z / SPAN)**2, 0, None)) - RHO * A_of_z(z) * G
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q_total[i] + q_total[i+1]) * dz
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dz
    curv = M / (E_MOD * I_of_z(z))
    theta = np.zeros(n); v = np.zeros(n)
    for i in range(1, n):
        theta[i] = theta[i-1] + 0.5 * (curv[i-1] + curv[i]) * dz
    for i in range(1, n):
        v[i] = v[i-1] + 0.5 * (theta[i-1] + theta[i]) * dz
    return z, v, M


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = (s[:, i] for i in range(6))
    return np.sqrt(0.5*((sxx-syy)**2+(syy-szz)**2+(szz-sxx)**2) + 3*(sxy**2+syz**2+szx**2))


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_stage3_solid.py <job>"); sys.exit(1)
    job = sys.argv[1]
    inp, frd = job + ".inp", job + ".frd"
    xyz, elems, cloads, root_set = read_inp(inp)
    n = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    print(f"Job: {job}\nNodes: {n}, C3D10 elements: {len(elems)}")
    print(f"Section: A(0)={A_of_z(0):.3f} mm^2, I(0)={I_of_z(0):.2f} mm^4; "
          f"A(1000)={A_of_z(SPAN):.3f}, I(1000)={I_of_z(SPAN):.2f}")

    U = read_frd_block(frd, "DISP", 3, n)
    RF = read_frd_block(frd, "FORC", 3, n)
    S = read_frd_block(frd, "STRESS", 6, n)
    vm = von_mises(S)
    z = coords[:, 2]

    root = np.nonzero(np.abs(z[1:]) < TOL)[0] + 1
    tip_band = np.nonzero(np.abs(z[1:] - SPAN) < 2.0)[0] + 1
    if set(root.tolist()) != set(root_set):
        print("WARNING: root nodes from geometry differ from NROOT in .inp")

    # ---- equilibrium ----
    zz = np.linspace(0, SPAN, 2001)
    total_weight = RHO * G * trapz(np.array([A_of_z(x) for x in zz]), zz)
    net_applied_Y = TOTAL_LIFT_N - total_weight
    r = RF[root].sum(axis=0)
    raw_expected = -net_applied_Y
    raw_err = 100 * (r[1] - raw_expected) / raw_expected

    rootset = set(root.tolist())
    cload_root = sum(f for nd, dof, f in cloads if dof == 2 and nd in rootset)
    P = coords[elems[:, :4]]
    vol_e = np.einsum("ij,ij->i", np.cross(P[:, 1]-P[:, 0], P[:, 2]-P[:, 0]), P[:, 3]-P[:, 0]) / 6.0
    grav_root = 0.0
    for e in range(len(elems)):
        for k, nd in enumerate(elems[e]):
            if nd in rootset:
                share = -1.0/20.0 if k < 4 else 1.0/5.0
                grav_root += -RHO * G * vol_e[e] * share      # Y force, N (negative = down)
    applied_on_root = cload_root + grav_root
    corr_expected = -(net_applied_Y - applied_on_root)
    corr_err = 100 * (r[1] - corr_expected) / corr_expected
    print("\n--- Equilibrium (Addendum A convention) ---")
    print(f"Root nodes: {len(root)}; sum RF root Y = {r[1]:.6f} N (X={r[0]:.6f}, Z={r[2]:.6f})")
    print(f"Lift {TOTAL_LIFT_N:.4f} N up; self-weight {total_weight:.4f} N down; net {net_applied_Y:.4f} N")
    print(f"Raw expected reaction {raw_expected:.6f} N -> error {raw_err:+.4f}%")
    print(f"External load applied directly on root nodes: CLOAD {cload_root:+.6f} N, "
          f"gravity share {grav_root:+.6f} N")
    print(f"Corrected expected reaction {corr_expected:.6f} N -> error {corr_err:+.4f}%")

    # ---- tip deflection ----
    z_ref, v_ref, M_ref = tapered_beam_reference()
    uy = U[tip_band, 1]
    pct = 100 * (uy.mean() - v_ref[-1]) / v_ref[-1]
    print("\n--- Tip deflection (Y) ---")
    print(f"Tip-band nodes: {len(tip_band)}; FEA Uy mean={uy.mean():.6f} "
          f"min={uy.min():.6f} max={uy.max():.6f} mm")
    print(f"Tapered numerical beam reference (spar v2 section): {v_ref[-1]:.6f} mm")
    print(f"FEA vs reference: {pct:+.4f}%")

    # ---- flange stress ----
    print("\n--- Flange sigma_zz vs -M(z)*y/I(z) (outer faces, |x|<3, band |z-z0|<1.5) ---")
    stress_diffs = {}
    for z0 in (100.0, 300.0, 500.0):
        for side, sgn in (("top", 1.0), ("bottom", -1.0)):
            ids = []
            for nd in range(1, n + 1):
                x_, y_, z_ = coords[nd]
                if abs(z_ - z0) < 1.5 and abs(x_) < 3.0 and abs(y_ - sgn * b.h_env(z_)) < 0.15:
                    ids.append(nd)
            if not ids:
                print(f"  z0={z0:.0f} {side}: no nodes found"); continue
            ids = np.array(ids)
            fea = S[ids, 2].mean()
            th = np.mean([-np.interp(coords[i, 2], z_ref, M_ref) * coords[i, 1]
                          / I_of_z(coords[i, 2]) for i in ids])
            d = 100 * (fea - th) / th
            stress_diffs[f"{side}_{int(z0)}"] = d
            print(f"  z0={z0:.0f} {side:6s}: {len(ids):4d} nodes, FEA {fea:9.4f} MPa, "
                  f"beam {th:9.4f} MPa, diff {d:+.3f}%")

    print("\n--- Peak von Mises (local, unconverged-type value; not a design number) ---")
    order = np.argsort(vm[1:])[::-1][:5] + 1
    for nd in order:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.4f} MPa at ({c[0]:7.3f},{c[1]:7.3f},{c[2]:7.3f})")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_results.json"
    with open(out, "w") as f:
        json.dump({
            "job": job, "nodes": n,
            "root_reaction_Y_N": float(r[1]),
            "equilibrium_err_raw_pct": float(raw_err),
            "equilibrium_err_corrected_pct": float(corr_err),
            "uy_tip_mean_mm": float(uy.mean()),
            "beam_ref_tip_mm": float(v_ref[-1]), "tip_pct_diff": float(pct),
            "stress_pct_diffs": stress_diffs,
            "peak_vm_MPa": float(vm[order[0]]),
        }, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
