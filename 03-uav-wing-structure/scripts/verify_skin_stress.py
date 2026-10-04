#!/usr/bin/env python3
"""
Verify the skin-stress interpretation and characterize hot spots for the limit-load skin solve.

(1) sigma_zz at CalculiX expanded nodes (IDs > number of *NODE entries) in mid-bay bands vs the
    beam stress -M(z)*y/I_tot(z), with I_tot = spar v2 + skin ring (mid-surface offset t/2),
    M(z) from the tapered beam with the same loads (lift, skin/solid/rib weights x n).
(2) classify the top-200 skin-node von Mises values by distance to the trailing edge and to the
    nearest rib station; statistics excluding the TE strip (< 15 mm from TE).

Usage (system Python):
  /usr/bin/python3 scripts/verify_skin_stress.py simulation/stage3_skin_s2p5_k8_t1_n3p8
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import (read_loadcase, parse_inp, read_frd, read_frd_nodes, ring_props,
                                 von_mises, A_of_z, I_of_z, ZTAB, RIBW, RHO, G, SPAN)

RIBS = np.array([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])


def moment_and_inertia(t, lift, ng):
    ptab = np.array([ring_props(z, t) for z in ZTAB])
    perim_tab, Iring_tab = ptab[:, 0], ptab[:, 1]
    n = 4000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = lift / (np.pi * SPAN / 4.0)
    A_skin = np.interp(z, ZTAB, perim_tab) * t
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN) ** 2, 0, None)) - RHO * G * ng * (A_of_z(z) + A_skin)
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i + 1] + 0.5 * (q[i] + q[i + 1]) * dz
    for zr, _, _, P in RIBW:
        V[z < zr - 1e-9] -= P * ng
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i + 1] + 0.5 * (V[i] + V[i + 1]) * dz
    return z, M, I_of_z(z) + np.interp(z, ZTAB, Iring_tab)


def main():
    job = sys.argv[1]
    lift, ng = read_loadcase(job + ".inp")
    if len(sys.argv) > 2 and sys.argv[2] == "lift_only":
        ng = 0.0          # job has no gravity (variant *_nograv): weights off in the beam reference
    xyz, solid, shell, root, cload, t_skin = parse_inp(job + ".inp")
    N0 = len(xyz)
    S = read_frd(job + ".frd", "STRESS", 6)
    fid, fxyz = read_frd_nodes(job + ".frd")
    vm = von_mises(S)
    zb, M, Itot = moment_and_inertia(t_skin, lift, ng)
    print(f"Load case: lift {lift:.3f} N, n={ng}; skin t={t_skin} mm")

    exp = fid > N0
    eid, exyz = fid[exp], fxyz[exp]

    print("\n(1) Expanded-node sigma_zz vs beam -M*y/I_tot, |x|<40 mm, |y|>3 mm, band |z-z0|<2:")
    for z0 in (100.0, 300.0, 500.0, 700.0):
        m = (np.abs(exyz[:, 2] - z0) < 2.0) & (np.abs(exyz[:, 0]) < 40.0) & (np.abs(exyz[:, 1]) > 3.0)
        ids, p = eid[m], exyz[m]
        sb = -np.interp(p[:, 2], zb, M) * p[:, 1] / np.interp(p[:, 2], zb, Itot)
        fea = S[ids, 2]
        r = fea / sb
        up = p[:, 1] > 0
        print(f"  z0={z0:5.0f}: n={len(ids):4d}  FEA/beam: median {np.median(r):.3f}  mean {r.mean():.3f}  "
              f"std {r.std():.3f} | upper median {np.median(r[up]):.3f} lower median {np.median(r[~up]):.3f} "
              f"| beam sigma at y=+max: {sb[up].min():.2f} MPa")
        ax = np.abs(p[:, 0])
        bins = []
        for lo, hi in ((0, 4), (4, 12), (12, 25), (25, 40)):
            mm = (ax >= lo) & (ax < hi)
            bins.append(f"{lo}-{hi}mm {np.median(r[mm]):.3f} (n={int(mm.sum())})" if mm.any()
                        else f"{lo}-{hi}mm n/a")
        print("            median ratio by |x|: " + ", ".join(bins))

    print("\n(2) Hot-spot classification (skin expanded nodes):")
    vm_e = vm[eid]
    chord = 300.0 - 0.15 * exyz[:, 2]
    dTE = 0.7 * chord - exyz[:, 0]
    dRib = np.abs(exyz[:, 2][:, None] - RIBS[None, :]).min(axis=1)
    order = np.argsort(vm_e)[::-1]
    for k in (50, 200):
        top = order[:k]
        print(f"  top {k:3d}: within 15 mm of TE: {int((dTE[top] < 15).sum()):3d}; within 3 mm of a rib: "
              f"{int((dRib[top] < 3).sum()):3d}; both: {int(((dTE[top] < 15) & (dRib[top] < 3)).sum()):3d}")
    print("  top 12 nodes:")
    for k in order[:12]:
        print(f"    VM {vm_e[k]:8.3f} MPa  at ({exyz[k,0]:8.3f},{exyz[k,1]:7.3f},{exyz[k,2]:8.3f})  "
              f"dTE {dTE[k]:6.2f} mm  nearest-rib dz {dRib[k]:6.2f} mm  "
              f"(rib z={RIBS[np.abs(RIBS - exyz[k,2]).argmin()]:.0f})")
    for label, m in (("all skin nodes", np.ones(len(eid), bool)),
                     ("excluding TE strip (dTE>=15)", dTE >= 15.0),
                     ("excluding TE strip and z<10", (dTE >= 15.0) & (exyz[:, 2] > 10.0))):
        v = vm_e[m]
        print(f"  {label:30s}: n={int(m.sum()):6d} max {v.max():8.3f}  p99.9 {np.percentile(v, 99.9):8.3f}  "
              f"p99 {np.percentile(v, 99):8.3f}  median {np.median(v):7.3f} MPa")
    zhs = exyz[order[:200], 2]
    print("  z of top-200 nodes by nearest rib station:",
          {int(r): int((np.abs(zhs - r) < 25).sum()) for r in RIBS})


if __name__ == "__main__":
    main()
