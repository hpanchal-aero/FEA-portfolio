#!/usr/bin/env python3
"""
Trailing-edge hot-spot report for a skin static solve: max von Mises in a small window at each
rib's TE (|z - z_rib| < 2 mm, within 6 mm of the TE), for skin (CalculiX expanded) nodes and
solid-only nodes, plus distance rings from the TE at rib z=800 (skin). Window node counts show
how fine the local mesh is. Peaks in this window are suspected to be mesh-dependent (sharp
rib-tip / skin-fold junction); this script only reports, it does not judge.

Usage (system Python):  /usr/bin/python3 scripts/te_hotspot_report.py <job>
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import parse_inp, read_frd, read_frd_nodes, von_mises

job = sys.argv[1]
xyz, solid, shell, root, cload, t = parse_inp(job + ".inp")
N0 = len(xyz)
coords = np.vstack([np.zeros((1, 3)), xyz])
S = read_frd(job + ".frd", "STRESS", 6)
vm = von_mises(S)
fid, fxyz = read_frd_nodes(job + ".frd")
exp = fid > N0
eid, exyz = fid[exp], fxyz[exp]
solid_ids = np.unique(solid)
solid_only = np.setdiff1d(solid_ids, np.unique(shell))
sxyz = coords[solid_only]


def dte(z, x):
    return 0.7 * (300.0 - 0.15 * z) - x


print(f"Job: {job} (skin t={t} mm)")
print(f"{'rib z':>6s} | {'skin: n':>7s} {'max VM':>8s} {'at (x,y,z)':>26s} | {'solid-only: n':>13s} {'max VM':>8s} {'at (x,y,z)':>26s}")
for zr in (200.0, 400.0, 600.0, 800.0, 1000.0):
    ms = (np.abs(exyz[:, 2] - zr) < 2.0) & (dte(exyz[:, 2], exyz[:, 0]) < 6.0)
    mo = (np.abs(sxyz[:, 2] - zr) < 2.0) & (dte(sxyz[:, 2], sxyz[:, 0]) < 6.0)
    row = f"{zr:6.0f} |"
    for m, ids, pts in ((ms, eid, exyz), (mo, solid_only, sxyz)):
        if m.any():
            v = vm[ids[m]]
            k = int(np.argmax(v))
            p = pts[m][k]
            row += f" {int(m.sum()):7d} {v.max():8.2f} ({p[0]:7.2f},{p[1]:6.2f},{p[2]:7.2f}) |"
        else:
            row += f" {0:7d} {'-':>8s} {'':>26s} |"
    print(row)

print("\nSkin, rib z=800 window (|z-800|<2): max VM by distance from the TE")
w = (np.abs(exyz[:, 2] - 800.0) < 2.0)
d = dte(exyz[:, 2], exyz[:, 0])
for lo, hi in ((-1.0, 1.0), (1.0, 3.0), (3.0, 6.0), (6.0, 15.0), (15.0, 40.0)):
    m = w & (d >= lo) & (d < hi)
    print(f"  {lo:5.1f}..{hi:5.1f} mm from TE: n={int(m.sum()):5d}  max VM "
          f"{vm[eid[m]].max() if m.any() else float('nan'):8.2f} MPa")
