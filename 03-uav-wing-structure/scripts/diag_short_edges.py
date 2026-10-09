#!/usr/bin/env python3
"""
Shortest edges of an imprint geometry (mesh/stage3_skin/stage3_skin_geom[_rs<N>].step, RIB_SPACING env), with centre,
length and the number of faces using each edge. Uses the distinct-edge key of diag_netgen_faces.py topo.
Usage: RIB_SPACING=50 /usr/bin/python3 -u scripts/diag_short_edges.py [<n_rows>]
"""
import sys
from collections import Counter
import numpy as np
from diag_netgen_faces import load, xyz


def main():
    nrow = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    shape, faces, jc, ja, jsh = load()
    cnt, length = Counter(), {}
    for f in faces:
        for e in f.edges:
            c = xyz(e.center)
            key = (round(c[0], 3), round(c[1], 3), round(c[2], 3), round(float(e.mass), 3))
            cnt[key] += 1
            length[key] = float(e.mass)
    keys = sorted(length, key=length.get)
    L = np.array([length[k] for k in keys])
    print(f"distinct edges {len(keys)}; length < 0.01 mm: {int((L < 0.01).sum())}, < 0.1: {int((L < 0.1).sum())}, "
          f"< 0.5: {int((L < 0.5).sum())}, < 1: {int((L < 1).sum())}, < 2: {int((L < 2).sum())}")
    print("   length [mm]   centre (x, y, z)                  faces")
    for k in keys[:nrow]:
        print(f"   {length[k]:10.5f}   ({k[0]:9.3f}, {k[1]:9.3f}, {k[2]:9.3f})   {cnt[k]}")


if __name__ == "__main__":
    main()
