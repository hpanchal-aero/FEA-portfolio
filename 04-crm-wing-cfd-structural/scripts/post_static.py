#!/usr/bin/env python3
"""Post-process a CalculiX static run of a uCRM deck: displacement summary and comparison with a reference tip deflection.

Usage: post_static.py /path/to/jobname [reference_tip_uz_in_m]      (without extension)

Force balance is NOT evaluated here. The printed RF at the print set NSPC is K*u, i.e. support reaction plus the external load
at those nodes (CLOAD, gravity of the shell elements, point-mass weight transferred through *EQUATION). A balance that ignores
gravity therefore does not close. Use check_equilibrium.py for the force balance (uCRM-9 coarse: closes to 2.5 N of 6.3e5 N).
"""
import sys
from pathlib import Path

import numpy as np

base = Path(sys.argv[1]).resolve()
inp, dat = base.with_suffix(".inp"), base.with_suffix(".dat")

coords, kw = {}, None
for ln in inp.read_text().splitlines():
    s = ln.strip()
    if not s or s.startswith("**"):
        continue
    if s.startswith("*"):
        kw = s.split(",")[0].upper()
        continue
    if kw == "*NODE":
        t = [x.strip() for x in s.split(",") if x.strip()]
        coords[int(t[0])] = np.array([float(t[1]), float(t[2]), float(t[3])])

disp, mode = {}, None
for ln in dat.read_text().splitlines():
    s = ln.strip()
    if s and s[0].isalpha():
        mode = "U" if s.startswith("displacements") else None
        continue
    tok = ln.split()
    if mode == "U" and len(tok) == 4:
        try:
            disp[int(tok[0])] = np.array([float(x) for x in tok[1:]])
        except ValueError:
            pass

if not disp:
    sys.exit("ABORT: no displacements parsed from the .dat file")

print(f"deck nodes {len(coords)}   displacement nodes parsed {len(disp)}")
ids = np.array(sorted(disp))
U = np.array([disp[i] for i in ids])
X = np.array([coords[i] for i in ids])
mag = np.linalg.norm(U, axis=1)
k = int(np.argmax(mag))
print("\nDISPLACEMENTS [m]")
print(f"  max |u| = {mag[k]:.4f} at node {ids[k]}  x=({X[k][0]:.3f},{X[k][1]:.3f},{X[k][2]:.3f})")
kz = int(np.argmax(U[:, 2]))
print(f"  max uz  = {U[kz, 2]:.4f} at node {ids[kz]}  y={X[kz][1]:.3f}")
ymax = X[:, 1].max()
for tol in (1e-3, 0.05):
    m = X[:, 1] >= ymax - tol
    print(f"  tip section (y >= ymax-{tol}): {m.sum():4d} nodes   uz mean {U[m, 2].mean():.4f}   "
          f"min {U[m, 2].min():.4f}   max {U[m, 2].max():.4f}")
print(f"  ymax = {ymax:.3f} m")
print("\n  NASTRAN reference (ALL loads: aero + gravity + engine + LE/TE + fuel): coarse 2.5239 m "
      "(2.5244 / 2.5247 for medium / fine: column order assumed from the xlsx, verify)")
if len(sys.argv) > 2:
    ref = float(sys.argv[2])
    mt = X[:, 1] >= ymax - 1e-3
    for lab, val in (("tip-section mean uz", U[mt, 2].mean()), ("tip-section max uz", U[mt, 2].max()),
                     ("global max uz", U[:, 2].max())):
        print(f"  {lab:20s} {val:.4f} m   vs {ref:.4f} m : {100 * (val - ref) / ref:+.2f} %")
print("  (a deck without gravity and masses, e.g. v1, is NOT comparable)")
