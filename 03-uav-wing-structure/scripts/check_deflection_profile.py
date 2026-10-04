#!/usr/bin/env python3
"""
Deflection-profile check for the skin solve: FEA spar-axis Uy(z) (solid nodes, |x|<4.5, band |z-z0|<2)
vs the tapered beam v(z) with I_tot = spar v2 + skin ring, same loads. Also finite-difference
curvature (step 100 mm) from both profiles.

Usage (system Python):
  /usr/bin/python3 scripts/check_deflection_profile.py <job> [lift_only]
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import read_loadcase, parse_inp, read_frd, E_MOD, SPAN
from verify_skin_stress import moment_and_inertia


def main():
    job = sys.argv[1]
    lift, ng = read_loadcase(job + ".inp")
    if len(sys.argv) > 2 and sys.argv[2] == "lift_only":
        ng = 0.0
    xyz, solid, shell, root, cload, t_skin = parse_inp(job + ".inp")
    coords = np.vstack([np.zeros((1, 3)), xyz])
    U = read_frd(job + ".frd", "DISP", 3)
    print(f"Load case: lift {lift:.3f} N, n={ng}, skin t={t_skin}")

    zb, M, Itot = moment_and_inertia(t_skin, lift, ng)
    curv = M / (E_MOD * Itot)
    dz = zb[1] - zb[0]
    th = np.concatenate([[0.0], np.cumsum(0.5 * (curv[1:] + curv[:-1]) * dz)])
    v = np.concatenate([[0.0], np.cumsum(0.5 * (th[1:] + th[:-1]) * dz)])

    sid = np.unique(solid)
    zs = np.arange(0.0, 1000.1, 100.0)
    fea = []
    for z0 in zs:
        m = (np.abs(coords[sid, 2] - z0) < 2.0) & (np.abs(coords[sid, 0]) < 4.5)
        fea.append(U[sid[m], 1].mean())
    fea = np.array(fea)
    bv = np.interp(zs, zb, v)

    def fd(a):
        return (a[2:] - 2 * a[1:-1] + a[:-2]) / 100.0 ** 2

    print(f"\n{'z':>6s} {'FEA Uy':>10s} {'beam v':>10s} {'FEA/beam':>9s} {'FEA-beam':>10s}")
    for z0, f, b_ in zip(zs, fea, bv):
        r = f / b_ if abs(b_) > 1e-9 else float("nan")
        print(f"{z0:6.0f} {f:10.5f} {b_:10.5f} {r:9.4f} {f - b_:+10.5f}")
    kf, kb = fd(fea), fd(bv)
    print(f"\nFinite-difference curvature (step 100 mm), FEA vs beam:")
    print(f"{'z':>6s} {'kappa FEA':>12s} {'kappa beam':>12s} {'ratio':>8s}")
    for z0, a, b_ in zip(zs[1:-1], kf, kb):
        print(f"{z0:6.0f} {a:12.4e} {b_:12.4e} {a / b_:8.3f}")


if __name__ == "__main__":
    main()
