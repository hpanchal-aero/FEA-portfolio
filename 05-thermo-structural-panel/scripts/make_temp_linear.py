#!/usr/bin/python3
"""Write a CalculiX *TEMPERATURE include file with a linear through-thickness temperature profile.

Usage: /usr/bin/python3 scripts/make_temp_linear.py MESH.inp T_BOTTOM T_TOP OUT.inp

T(z) = T_bottom + (T_top - T_bottom) z / h for every node (uniform in x and y), z = 0 at the bottom face.
The output starts with the *TEMPERATURE keyword and must be included inside a *STEP.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    mesh, tb, tt, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    h = ccxio.read_params(mesh)["h"]
    nodes = ccxio.read_nodes(mesh)
    with open(out, "w") as f:
        f.write("** linear through-thickness temperature, T_bottom = %g, T_top = %g degC (scripts/make_temp_linear.py)\n" % (tb, tt))
        f.write("*TEMPERATURE\n")
        for i in sorted(nodes):
            z = nodes[i][2]
            f.write("%d, %.9f\n" % (i, tb + (tt - tb) * z / h))
    print("written %s (%d nodes)" % (out, len(nodes)))


if __name__ == "__main__":
    main()
