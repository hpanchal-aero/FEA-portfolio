#!/usr/bin/env python3
"""Inventory of the material orientation systems (CORD2R, one per PSHELL) of a uCRM-9 wingbox."""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

level = sys.argv[1] if len(sys.argv) > 1 else "coarse"
fem = (Path(__file__).resolve().parent.parent / "data" / "ucrm9" / "uCRM9" / "Mesh"
       / level.capitalize() / "FEM")
cards = read_cards(fem / f"uCRM-9_wingbox_{level}.bdf")

grids, quads, cord = {}, [], {}
for name, f in cards:
    if name == "GRID":
        grids[ival(f[0])] = np.array([num(f[2]), num(f[3]), num(f[4])])
    elif name == "CQUADR":
        quads.append((ival(f[0]), ival(f[1]), [ival(f[i]) for i in range(2, 6)], ival(f[6])))
    elif name == "CORD2R":
        cord[ival(f[0])] = [np.array([num(f[i]) for i in range(j, j + 3)]) for j in (2, 5, 8)]


def axes(cid):
    A, B, C = cord[cid]
    z = (B - A) / np.linalg.norm(B - A)
    x = (C - A) - np.dot(C - A, z) * z
    x /= np.linalg.norm(x)
    return x, np.cross(z, x), z


groups = defaultdict(lambda: dict(pids=set(), n=0, s=[], nz=[]))
for eid, pid, ns, mcid in quads:
    x, y, z = axes(mcid)
    p1, p2, p3, p4 = (grids[n] for n in ns)
    nv = np.cross(p3 - p1, p4 - p2)
    nv /= np.linalg.norm(nv)
    key = tuple(np.round(np.concatenate([x, y]), 5))
    g = groups[key]
    g["pids"].add(pid)
    g["n"] += 1
    g["s"].append(np.sqrt(max(0.0, 1.0 - np.dot(x, nv) ** 2)))   # |projection of x on the element plane|
    g["nz"].append(abs(np.dot(nv, np.array([0.0, 0.0, 1.0]))))     # |n_z|: ~1 skin panels, ~0 webs/ribs

print(f"level {level}: {len(cord)} CORD2R systems, {len(groups)} distinct (x,y) orientations\n")
print("  x-axis                 y-axis               #props  #elems   min|proj x|  mean|n.z|")
for key, g in sorted(groups.items(), key=lambda kv: -kv[1]["n"]):
    x, y = key[:3], key[3:]
    print(f"  ({x[0]:6.3f},{x[1]:6.3f},{x[2]:6.3f})  ({y[0]:6.3f},{y[1]:6.3f},{y[2]:6.3f})  "
          f"{len(g['pids']):5d}  {g['n']:6d}   {min(g['s']):9.4f}   {np.mean(g['nz']):8.3f}")

all_s = np.concatenate([g["s"] for g in groups.values()])
print(f"\nelements with |projection of x| < 0.1 : {(all_s < 0.1).sum()} of {len(all_s)}")
print(f"elements with |projection of x| < 0.5 : {(all_s < 0.5).sum()} of {len(all_s)}")
print(f"overall minimum |projection of x|      : {all_s.min():.4f}")
