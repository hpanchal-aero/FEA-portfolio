#!/usr/bin/env python3
"""Post-process a CalculiX static run of a uCRM deck: equilibrium and displacement summary.

Reaction rule: the printed RF is K*u at every node of the print set. At a loaded free DOF it equals the
applied load; at a constrained DOF it equals (support reaction + applied load). The support reaction is
therefore: printed total - applied load on free DOFs of NSPC nodes - applied load on constrained DOFs.

Usage: post_static.py /path/to/jobname [reference_tip_uz_in_m]      (without extension)

Note on reactions: CalculiX prints, at an UNCONSTRAINED but loaded DOF of a node in the print set,
an "RF" equal to the applied load. The deck's print set NSPC contains nodes that are constrained in
only some DOFs (rib at y = 3.1 m: z only) and also loaded in x and y, so the printed total is
corrected below by the applied load on the free DOFs of NSPC nodes. (Verified on uCRM-9 coarse v1:
the correction closes the x and y residuals to the printed precision.)
"""
import sys
from pathlib import Path

import numpy as np

base = Path(sys.argv[1]).resolve()
inp, dat = base.with_suffix(".inp"), base.with_suffix(".dat")

coords, cloads, nspc, cons = {}, [], set(), {}
kw, cur_set = None, None
for ln in inp.read_text().splitlines():
    s = ln.strip()
    if s.startswith("*"):
        kw = s.split(",")[0].upper()
        cur_set = s.upper().split("NSET=")[1].split(",")[0].strip() if (kw == "*NSET" and "NSET=" in s.upper()) else None
        continue
    if not s:
        continue
    t = [x.strip() for x in s.split(",") if x.strip()]
    if kw == "*NODE":
        coords[int(t[0])] = np.array([float(t[1]), float(t[2]), float(t[3])])
    elif kw == "*NSET" and cur_set == "NSPC":
        nspc.update(int(x) for x in t)
    elif kw == "*CLOAD":
        cloads.append((int(t[0]), int(t[1]), float(t[2])))
    elif kw == "*BOUNDARY":
        d1 = int(t[1])
        d2 = int(t[2]) if len(t) > 2 else d1
        cons.setdefault(int(t[0]), set()).update(range(d1, d2 + 1))

load = np.zeros(3)
load_on_cons = np.zeros(3)
load_free_nspc = np.zeros(3)
for n, d, v in cloads:
    load[d - 1] += v
    if d in cons.get(n, set()):
        load_on_cons[d - 1] += v
    elif n in nspc:
        load_free_nspc[d - 1] += v

disp, rf_tot, mode = {}, None, None
for ln in dat.read_text().splitlines():
    s = ln.strip()
    if s.startswith("displacements"):
        mode = "U"
        continue
    if s.startswith("total force"):
        mode = "RFT"
        continue
    tok = ln.split()
    if mode == "U" and len(tok) == 4:
        try:
            disp[int(tok[0])] = np.array([float(x) for x in tok[1:]])
        except ValueError:
            pass
    elif mode == "RFT" and len(tok) == 3:
        rf_tot = np.array([float(x) for x in tok])
        mode = None

if not disp or rf_tot is None:
    sys.exit(f"ABORT: parsed {len(disp)} displacement nodes, reaction total {rf_tot}")

print(f"deck nodes {len(coords)}   displacement nodes parsed {len(disp)}   NSPC {len(nspc)}")
print("\nEQUILIBRIUM")
rf_corr = rf_tot - load_free_nspc - load_on_cons
print(f"  applied load                  : Fx={load[0]:+.4e}  Fy={load[1]:+.4e}  Fz={load[2]:+.4e} N")
print(f"  reaction total (printed)      : Fx={rf_tot[0]:+.4e}  Fy={rf_tot[1]:+.4e}  Fz={rf_tot[2]:+.4e} N")
print(f"  load on free DOFs of NSPC nodes: Fx={load_free_nspc[0]:+.4e}  Fy={load_free_nspc[1]:+.4e}  Fz={load_free_nspc[2]:+.4e} N")
print(f"  reaction total (corrected)    : Fx={rf_corr[0]:+.4e}  Fy={rf_corr[1]:+.4e}  Fz={rf_corr[2]:+.4e} N")
res = load + rf_corr
print(f"  load + corrected reaction     : Fx={res[0]:+.3e}  Fy={res[1]:+.3e}  Fz={res[2]:+.3e} N   "
      f"(max |residual| = {np.abs(res).max():.3e} N = {100 * np.abs(res).max() / abs(load[2]):.2e} % of |Fz|)")
print(f"  load applied on constrained DOFs: Fx={load_on_cons[0]:+.4e}  Fy={load_on_cons[1]:+.4e}  "
      f"Fz={load_on_cons[2]:+.4e} N  (subtracted above)")

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
