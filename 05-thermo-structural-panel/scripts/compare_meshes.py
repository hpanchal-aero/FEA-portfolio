#!/usr/bin/python3
"""Compare two panel meshes: identical node ids and coordinates, identical element ids and connectivity.
The element type label may differ (C3D20 versus C3D20R).

Usage: /usr/bin/python3 scripts/compare_meshes.py TAG_A TAG_B
  e.g. /usr/bin/python3 scripts/compare_meshes.py v6c_n20_z2 v8a_n20_z2
Exits with status 1 if the meshes differ.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio


def element_labels(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.upper().startswith("*ELEMENT")]


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    root = os.path.dirname(HERE)
    pa, pb = [os.path.join(root, "mesh", t + ".inp") for t in sys.argv[1:3]]
    na, nb = ccxio.read_nodes(pa), ccxio.read_nodes(pb)
    ea, ca = ccxio.read_elements(pa)
    eb, cb = ccxio.read_elements(pb)
    print("%s: %d nodes, %d elements; %s" % (sys.argv[1], len(na), len(ea), element_labels(pa)))
    print("%s: %d nodes, %d elements; %s" % (sys.argv[2], len(nb), len(eb), element_labels(pb)))
    ok = True
    if set(na) != set(nb):
        print("node id sets differ")
        ok = False
    else:
        dmax = max(float(np.abs(np.array(na[k], dtype=float) - np.array(nb[k], dtype=float)).max()) for k in na)
        print("max node coordinate difference %.3e m" % dmax)
        ok = ok and dmax < 1e-12
    same_el = list(map(int, ea)) == list(map(int, eb)) and all([int(x) for x in a] == [int(x) for x in b] for a, b in zip(ca, cb))
    print("element ids and connectivity identical: %s" % same_el)
    ok = ok and same_el
    print("RESULT: %s" % ("IDENTICAL" if ok else "DIFFERENT"))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
