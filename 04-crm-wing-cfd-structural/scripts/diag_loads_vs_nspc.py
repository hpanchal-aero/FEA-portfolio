#!/usr/bin/env python3
"""Diagnostic: relate loaded nodes to the constrained node set NSPC and to the equilibrium residual.

Usage: diag_loads_vs_nspc.py /path/to/deck.inp
"""
import sys
from pathlib import Path

import numpy as np

inp = Path(sys.argv[1]).resolve()
coords, cloads, nspc, cons, kw, cur_set = {}, [], set(), {}, None, None
for ln in inp.read_text().splitlines():
    s = ln.strip()
    if s.startswith("*"):
        kw = s.split(",")[0].upper()
        cur_set = None
        if kw == "*NSET":
            cur_set = s.upper().split("NSET=")[1].split(",")[0].strip()
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

loaded = sorted({n for n, _, _ in cloads})
print(f"nodes {len(coords)}   NSPC {len(nspc)}   loaded nodes {len(loaded)}")

in_set = [n for n in loaded if n in nspc]
print(f"loaded nodes that belong to NSPC: {len(in_set)}")
if in_set:
    tot_all = np.zeros(3)
    tot_free = np.zeros(3)
    for n, d, v in cloads:
        if n in nspc:
            tot_all[d - 1] += v
            if d not in cons.get(n, set()):
                tot_free[d - 1] += v
    print(f"  load sum on those nodes, all DOFs : Fx={tot_all[0]:+.4e}  Fy={tot_all[1]:+.4e}  Fz={tot_all[2]:+.4e} N")
    print(f"  ... on DOFs that are NOT constrained: Fx={tot_free[0]:+.4e}  Fy={tot_free[1]:+.4e}  Fz={tot_free[2]:+.4e} N")
    print("  (observed residual: Fx=+3.742e+02  Fy=-1.800e+03  Fz=+1.8e-02 N)")
    ys = np.array([coords[n][1] for n in in_set])
    print(f"  y range of those nodes: {ys.min():.3f}..{ys.max():.3f} m")

print("\nnodes 340-344 (largest early y-only forces):")
for n in range(340, 345):
    ld = [(d, v) for m, d, v in cloads if m == n]
    print(f"  node {n}: xyz=({coords[n][0]:.3f},{coords[n][1]:.3f},{coords[n][2]:.3f})  in NSPC: {n in nspc}   loads {ld}")

yl = np.array([coords[n][1] for n in loaded])
print(f"\nloaded nodes: y range {yl.min():.3f}..{yl.max():.3f} m   (y<0.01: {(yl < 0.01).sum()},  "
      f"|y-3.1|<0.01: {(abs(yl - 3.1) < 0.01).sum()})")
