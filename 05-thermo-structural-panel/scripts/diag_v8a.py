#!/usr/bin/python3
"""V8a diagnostic: top and bottom temperature along the line y = L/2 near the plate centre, final output time.

Usage: /usr/bin/python3 scripts/diag_v8a.py JOB MESHTAG
  e.g. /usr/bin/python3 scripts/diag_v8a.py v8a_n20_dt0p5 v6c_n20_z2
Reads simulation/JOB/JOB.dat (NALL temperatures) and the mesh. C = element-corner column, M = mid-side column
(in-plane node index even = corner). Prints x, kind, T_top, T_bottom, T_top - T_bottom for the 9 nodes ending at the centre.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    job, mesh = sys.argv[1], sys.argv[2]
    root = os.path.dirname(HERE)
    mp = os.path.join(root, "mesh", mesh + ".inp")
    par = ccxio.read_params(mp)
    L, h = par["L"], par["h"]
    nodes = ccxio.read_nodes(mp)
    dat = os.path.join(root, "simulation", job, job + ".dat")
    blocks = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("temperatures") and b["set"] == "NALL"]
    b = max(blocks, key=lambda blk: blk["time"])
    T = dict((int(r[0]), r[1]) for r in b["rows"])
    xs = sorted(set(round(p[0], 9) for p in nodes.values()))
    dx = xs[1] - xs[0]
    yc = L / 2.0
    top, bot = {}, {}
    for n, p in nodes.items():
        if abs(p[1] - yc) > 1e-9:
            continue
        i = int(round(p[0] / dx))
        if abs(p[2] - h) < 1e-9:
            top[i] = T[n]
        elif abs(p[2]) < 1e-9:
            bot[i] = T[n]
    ic = int(round(yc / dx))
    print("%s, t = %.2f s, in-plane node spacing %.5f m, line y = %.4f m (node row index %d)" % (job, b["time"], dx, yc, ic))
    print("  index        x [m]  kind   T_top [C]   T_bottom [C]   top - bottom [K]")
    for i in range(ic - 8, ic + 1):
        kind = "C" if i % 2 == 0 else "M"
        print("  %5d %12.5f   %s %11.4f %14.4f %16.4f" % (i, i * dx, kind, top[i], bot[i], top[i] - bot[i]))


if __name__ == "__main__":
    main()
