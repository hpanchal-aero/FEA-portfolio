#!/usr/bin/env python3
"""Read-only: plate-bending stiffness of NASTRAN PSHELL properties relative to the homogeneous plate that
the CalculiX deck uses, split by element class and span band.

Per property (NASTRAN PSHELL convention, as in section_stiffness.py):
  A = t G(MID1), B = t^2 G(MID4), D = (12I/T^3) t^3/12 G(MID2)
CalculiX deck plate: D_hom,ii = A_ii t^2/12 (shell stiffness from the membrane matrix only).
Ratios in the property's material axes (1, 2).
Usage: /usr/bin/python3 scripts/audit_rib_bending.py [coarse|medium|fine]
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

root = HERE.parent
level = sys.argv[1] if len(sys.argv) > 1 else "coarse"
fem = root / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"


def fnum(s, default=0.0):
    return num(s) if s.strip() else default


grids, quads, pshell, mat2 = {}, [], {}, {}
for name, f in read_cards(fem / f"uCRM-9_wingbox_{level}.bdf"):
    if name == "GRID":
        grids[ival(f[0])] = np.array([num(f[2]), num(f[3]), num(f[4])])
    elif name == "CQUADR":
        quads.append((ival(f[0]), ival(f[1]), [ival(f[i]) for i in range(2, 6)]))
    elif name == "PSHELL":
        pshell[ival(f[0])] = dict(mid1=ival(f[1]), t=num(f[2]), mid2=ival(f[3]),
                                  r=fnum(f[4], 1.0), mid4=ival(f[10]))
    elif name == "MAT2":
        mat2[ival(f[0])] = np.array([[fnum(f[1]), fnum(f[2]), fnum(f[3])],
                                     [fnum(f[2]), fnum(f[4]), fnum(f[5])],
                                     [fnum(f[3]), fnum(f[5]), fnum(f[6])]])

prop = {}
for pid, p in pshell.items():
    t = p["t"]
    A = t * mat2[p["mid1"]]
    B = t ** 2 * mat2[p["mid4"]]
    D = p["r"] * t ** 3 / 12.0 * mat2[p["mid2"]]
    Dh = np.array([A[0, 0], A[1, 1]]) * t * t / 12.0
    prop[pid] = dict(t=t, r=p["r"], d1=D[0, 0] / Dh[0], d2=D[1, 1] / Dh[1],
                     c1=B[0, 0] ** 2 / (A[0, 0] * D[0, 0]) if D[0, 0] > 0 else np.nan,
                     c2=B[1, 1] ** 2 / (A[1, 1] * D[1, 1]) if D[1, 1] > 0 else np.nan)

groups = defaultdict(list)          # name -> list of pids (one per element)
for eid, pid, ns in quads:
    P = np.array([grids[n] for n in ns])
    n = np.cross(P[2] - P[0], P[3] - P[1])
    nn = np.linalg.norm(n)
    if nn < 1e-12:
        groups["degenerate"].append(pid)
        continue
    n /= nn
    yc = P[:, 1].mean()
    if abs(n[1]) >= 0.5:
        groups["ribs: all"].append(pid)
        if 0.70 <= yc <= 0.85:
            groups["ribs: y 0.70-0.85"].append(pid)
        elif 1.45 <= yc <= 1.65:
            groups["ribs: y 1.45-1.65"].append(pid)
    else:
        groups["non-rib (skins, spars)"].append(pid)


def stat(pids, key):
    v = np.array([prop[p][key] for p in pids], float)
    v = v[np.isfinite(v)]
    return f"{np.median(v):9.3g} [{v.min():9.3g}, {v.max():9.3g}]"


print(f"level {level}: {len(quads)} elements, {len(pshell)} properties\n")
print("median [min, max] over elements")
print(f"{'group':26s} {'n':>6s} {'D1/Dhom1':>30s} {'D2/Dhom2':>30s} {'12I/T3':>26s}")
for g in ("non-rib (skins, spars)", "ribs: all", "ribs: y 0.70-0.85", "ribs: y 1.45-1.65"):
    pids = groups.get(g, [])
    if not pids:
        print(f"{g:26s} {0:6d}  (none)")
        continue
    print(f"{g:26s} {len(pids):6d} {stat(pids, 'd1'):>30s} {stat(pids, 'd2'):>30s} {stat(pids, 'r'):>26s}")
print(f"{'degenerate':26s} {len(groups.get('degenerate', [])):6d}")

for g in ("ribs: y 0.70-0.85", "ribs: y 1.45-1.65"):
    cnt = Counter(groups.get(g, []))
    print(f"\nproperties in {g} (pid, elements, t [mm], 12I/T3, D1/Dhom1, D2/Dhom2, B1^2/(A1 D1))")
    for pid, k in cnt.most_common(8):
        q = prop[pid]
        print(f"  {pid:5d} {k:5d} {q['t'] * 1e3:8.3f} {q['r']:8.3f} {q['d1']:10.3g} {q['d2']:10.3g} {q['c1']:8.3f}")
