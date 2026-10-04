#!/usr/bin/env python3
"""
Project 03, Stage 3 baseline (no skin) -- spar v2 + 6 ribs. Extract CalculiX
results and verify against a tapered Euler-Bernoulli beam reference using the
spar v2 section, with the ribs as POINT WEIGHTS (zero bending stiffness), as
established in Stage 2.

Rib added weight at each station = (solid volume in the 1.5mm rib slab, i.e.
envelope area integrated over the slab) - (spar volume inside that slab that
the beam's continuous A(z) already counts). Self-check: spar volume + sum of
additions must reproduce the mesh volume.

Tip metric: spar-axis deflection (|x|<4.5mm at the tip band). Equilibrium uses
the CalculiX constrained-node reaction convention (Addendum A).

Usage (system Python):
  /usr/bin/python3 scripts/extract_stage3_ribs.py <job_path_without_ext>
"""
import json
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_stage3_solid as b
from extract_stage3_solid import read_frd_block, von_mises, A_of_z, I_of_z
from extract_stage2_ribs import parse_inp, COEF

SPAN = 1000.0
RIB_T = b.RIB_T
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
TOTAL_LIFT_N = 24.5
TOL = 1e-3
ASSEMBLY_VOLUME_MM3 = 128582.1216      # build_stage3_solid.py ribs

# spar v2 ALONE, FEA tip deflection (extract_stage3_solid.py runs, same mesh sizes)
SPAR_V2_ALONE_TIP_MM = {"h2p0": 2.158439, "h2p5": 2.158980}


def trapz(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def rib_point_weights():
    out = []
    for zs in b.RIB_STATIONS:
        a, c = zs - RIB_T / 2.0, zs + RIB_T / 2.0
        v_env = b.simpson(b.env_area, a, c, 2)
        lo, hi = max(a, 0.0), min(c, SPAN)
        v_spar = b.simpson(b.spar_area, lo, hi, 20) if hi > lo else 0.0
        out.append((zs, v_env, v_spar, (v_env - v_spar) * RHO * G))
    return out


RIBW = rib_point_weights()


def beam(ribs_on):
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = TOTAL_LIFT_N / (np.pi * SPAN / 4.0)
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN)**2, 0, None)) - RHO * A_of_z(z) * G
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q[i] + q[i+1]) * dz
    if ribs_on:
        for zr, _, _, P in RIBW:
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
    return z, v, M


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_stage3_ribs.py <job>"); sys.exit(1)
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
    print(f"Fixed set (NROOT): {len(root)} nodes, z range [{zr.min():.6f}, {zr.max():.6f}] "
          f"(mesh z_min={z_min:.6f}, expected {-RIB_T/2:.4f})")

    P = coords[elems[:, :4]]
    vol = np.abs(np.einsum("ij,ij->i",
                           np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                           P[:, 3] - P[:, 0])) / 6.0
    total_vol = vol.sum()
    weight = RHO * G * total_vol
    print(f"Mesh volume {total_vol:.4f} mm^3 (geometry {ASSEMBLY_VOLUME_MM3:.4f}, "
          f"diff {100*(total_vol - ASSEMBLY_VOLUME_MM3)/ASSEMBLY_VOLUME_MM3:+.4f}%)")

    # reference self-check: spar(0..SPAN) + rib additions == mesh volume
    spar_vol = b.simpson(b.spar_area, 0.0, SPAN, 400)
    add_vol = sum(w[1] - w[2] for w in RIBW)
    ref_vol = spar_vol + add_vol
    print("\n--- Rib point weights for the beam reference ---")
    for zs, v_env, v_spar, Pw in RIBW:
        print(f"  z={zs:6.0f}: slab solid {v_env:10.4f} mm^3, spar in slab {v_spar:9.4f}, "
              f"added weight {Pw:.6f} N")
    print(f"Reference mass check: spar {spar_vol:.4f} + rib additions {add_vol:.4f} = "
          f"{ref_vol:.4f} mm^3 vs mesh {total_vol:.4f} "
          f"({100*(ref_vol - total_vol)/total_vol:+.4f}%)")

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
    print(f"Load applied directly on fixed nodes: gravity share {F_grav_root:+.6f} N, "
          f"CLOAD {F_lift_root:+.6f} N (not in reported RF)")
    eq_naive = 100 * (r[1] - expected_naive) / abs(expected_naive)
    eq_corr = 100 * (r[1] - expected_corr) / abs(expected_corr)
    print(f"Expected reaction, naive:      {expected_naive:+.6f} N -> error {eq_naive:+.4f}%")
    print(f"Expected reaction, corrected:  {expected_corr:+.6f} N -> error {eq_corr:+.4f}%  (PRIMARY)")

    print("\n--- Tip deflection (Y) ---")
    tip_band = np.nonzero(np.abs(z[1:] - SPAN) < 2.0)[0] + 1
    axis = tip_band[np.abs(coords[tip_band, 0]) < 4.5]
    slope, intercept = np.polyfit(coords[tip_band, 0], U[tip_band, 1], 1)
    uy_axis = U[axis, 1].mean()
    print(f"Tip-band nodes: {len(tip_band)}; spar-axis nodes (|x|<4.5): {len(axis)}")
    print(f"Spar-axis FEA Uy (PRIMARY): {uy_axis:.6f} mm  "
          f"(min {U[axis,1].min():.6f}, max {U[axis,1].max():.6f})")
    print(f"Regression Uy = a + b*x: a={intercept:.6f} mm, b={slope:+.6f} mm/mm (twist-like)")
    print(f"(Tip-band mean over all nodes, NOT used: {U[tip_band,1].mean():.6f} mm)")

    z_ref, v_ribs, M_ribs = beam(True)
    _, v_spar_only, _ = beam(False)
    v_ribw = v_ribs[-1]
    pct = 100 * (uy_axis - v_ribw) / v_ribw
    print(f"\nReference (spar v2 beam, ribs as point weights): {v_ribw:.6f} mm")
    print(f"FEA spar-axis vs reference: {pct:+.4f}%")
    print(f"(Reference spar v2 alone: {v_spar_only[-1]:.6f} mm; rib weight alone changes it "
          f"by {100*(v_ribw - v_spar_only[-1])/v_spar_only[-1]:+.2f}%)")
    m = re.search(r"_(h\dp\d)$", job)
    key = m.group(1) if m else None
    if key in SPAR_V2_ALONE_TIP_MM:
        alone = SPAR_V2_ALONE_TIP_MM[key]
        d = uy_axis - alone
        print(f"Change vs spar v2 ALONE FEA at same mesh ({alone:.6f} mm): {d:+.6f} mm "
              f"({100*d/alone:+.3f}%)")

    print("\n--- Flange sigma_zz vs -M(z)*y/I(z) (outer faces, |x|<3, band |z-z0|<1.5) ---")
    stress_diffs = {}
    for z0 in (100.0, 300.0, 500.0):
        for side, sgn in (("top", 1.0), ("bottom", -1.0)):
            ids = [nd for nd in range(1, n + 1)
                   if abs(z[nd] - z0) < 1.5 and abs(coords[nd, 0]) < 3.0
                   and abs(coords[nd, 1] - sgn * b.h_env(z[nd])) < 0.15]
            if not ids:
                print(f"  z0={z0:.0f} {side}: no nodes found"); continue
            ids = np.array(ids)
            fea = S[ids, 2].mean()
            th = np.mean([-np.interp(z[i], z_ref, M_ribs) * coords[i, 1] / I_of_z(z[i]) for i in ids])
            dd = 100 * (fea - th) / th
            stress_diffs[f"{side}_{int(z0)}"] = dd
            print(f"  z0={z0:.0f} {side:6s}: {len(ids):4d} nodes, FEA {fea:9.4f} MPa, "
                  f"beam {th:9.4f} MPa, diff {dd:+.3f}%")

    print("\n--- Peak von Mises (local corner peak; convergence NOT established) ---")
    order = np.argsort(vm[1:])[::-1][:5] + 1
    for nd in order:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.4f} MPa at ({c[0]:7.3f},{c[1]:7.3f},{c[2]:7.3f})")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_results.json"
    with open(out, "w") as f:
        json.dump({
            "job": job, "nodes": n, "elements": int(len(elems)),
            "root_reaction_Y_N": float(r[1]),
            "equilibrium_err_naive_pct": float(eq_naive),
            "equilibrium_err_corrected_pct": float(eq_corr),
            "uy_spar_axis_mm": float(uy_axis),
            "twist_slope_mm_per_mm": float(slope),
            "reference_ribs_as_weights_mm": float(v_ribw),
            "pct_vs_reference": float(pct),
            "reference_spar_alone_mm": float(v_spar_only[-1]),
            "stress_pct_diffs": stress_diffs,
            "peak_vm_MPa": float(vm[order[0]]),
        }, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
