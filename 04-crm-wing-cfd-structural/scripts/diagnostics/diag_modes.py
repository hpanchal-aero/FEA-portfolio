#!/usr/bin/env python3
"""Read-only: how localized are the eigenmodes in a CalculiX modal .frd?

usage: diag_modes.py <job.frd>

For each mode: max |u|, number of nodes above 10% and 50% of the maximum, the bounding box of the
nodes above 10%, and the coordinates of the node with the largest displacement.
Note: .frd node ids are those of the expanded model, so ids are not compared with the input deck.
"""
import math
import sys
from pathlib import Path

frd = Path(sys.argv[1])
coords, modes = {}, []
blk, cur = None, None
with frd.open() as fh:
    for line in fh:
        if line.startswith("    2C"):
            blk = "coord"
            continue
        if line.startswith(" -4  DISP"):
            blk = "disp"
            cur = {}
            modes.append(cur)
            continue
        if line.startswith(" -3"):
            blk = None
            continue
        if line.startswith(" -1") and blk:
            try:
                nid = int(line[3:13])
                vals = tuple(float(line[13 + 12 * i: 25 + 12 * i]) for i in range(3))
            except ValueError:
                continue
            if blk == "coord":
                coords[nid] = vals
            else:
                cur[nid] = vals

if not modes:
    sys.exit("ABORT: no DISP blocks found")
print(f"DISP blocks: {len(modes)}   nodes with coordinates: {len(coords)}   nodes per block: {len(modes[0])}")
print(f"{'mode':>4} {'max|u|':>10} {'n>10%':>7} {'n>50%':>7}   bounding box of nodes >10% (x, y, z)")
for k, u in enumerate(modes, 1):
    mag = {n: math.sqrt(sum(c * c for c in v)) for n, v in u.items()}
    top = max(mag, key=mag.get)
    umax = mag[top]
    big = [n for n in mag if mag[n] > 0.1 * umax and n in coords]
    half = sum(1 for n in mag if mag[n] > 0.5 * umax)
    if big:
        xs = [coords[n] for n in big]
        bb = "  ".join(f"[{min(c[i] for c in xs):7.2f},{max(c[i] for c in xs):7.2f}]" for i in range(3))
    else:
        bb = "(no coordinates)"
    print(f"{k:>4} {umax:>10.3e} {len(big):>7} {half:>7}   {bb}")
    x = coords.get(top)
    if x:
        print(f"{'':>4}  largest at ({x[0]:.3f}, {x[1]:.3f}, {x[2]:.3f})")
