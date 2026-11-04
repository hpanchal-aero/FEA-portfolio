#!/usr/bin/python3
"""Write a CalculiX *DFLUX include file for the top face (z = h) of the panel mesh.

Usage: /usr/bin/python3 scripts/make_flux_top.py MESH.inp Q0 OUT.inp AMPLITUDE_NAME [--gauss SIGMA XC YC]

Flux q = Q0 (uniform), or q = Q0 exp(-((x-xc)^2 + (y-yc)^2) / (2 sigma^2)) evaluated at the centre of
each top element face (piecewise constant per face). The face of each top-layer element is found from the
geometry (the hexahedron face whose four corner nodes all lie at z = h), not assumed from the numbering.
The file starts with *DFLUX, AMPLITUDE=NAME and must be included inside a *STEP; the amplitude is defined
in the deck. Prints the number of faces, the total power sum(q A) and, for the uniform case, Q0 L W.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio

# corner-node indices (0-based) of the six faces of a hexahedron, CalculiX/Abaqus numbering S1..S6
FACES = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]
TOL = 1e-9


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh")
    ap.add_argument("q0", type=float)
    ap.add_argument("out")
    ap.add_argument("amp")
    ap.add_argument("--gauss", type=float, nargs=3, metavar=("SIGMA", "XC", "YC"))
    a = ap.parse_args()
    par = ccxio.read_params(a.mesh)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(a.mesh)
    eids, conn = ccxio.read_elements(a.mesh)
    lines, power, nfaces = [], 0.0, 0
    for e, c in zip(eids, conn):
        P = np.array([nodes[int(n)] for n in c[:8]])
        found = [k for k, f in enumerate(FACES) if np.all(np.abs(P[list(f), 2] - h) < TOL)]
        if not found:
            continue
        if len(found) != 1:
            raise SystemExit("element %d has %d faces at z = h" % (e, len(found)))
        k = found[0]
        q = P[list(FACES[k])]
        cen = q.mean(axis=0)
        area = 0.5 * np.linalg.norm(np.cross(q[2] - q[0], q[3] - q[1]))
        val = a.q0
        if a.gauss:
            s, xc, yc = a.gauss
            val = a.q0 * np.exp(-((cen[0] - xc) ** 2 + (cen[1] - yc) ** 2) / (2 * s * s))
        lines.append("%d, S%d, %.9e\n" % (e, k + 1, val))
        power += val * area
        nfaces += 1
    with open(a.out, "w") as f:
        f.write("** top-face heat flux, Q0 = %g W/m2%s (scripts/make_flux_top.py)\n"
                % (a.q0, "" if not a.gauss else ", Gaussian sigma %g xc %g yc %g" % tuple(a.gauss)))
        f.write("*DFLUX, AMPLITUDE=%s\n" % a.amp)
        f.writelines(lines)
    print("written %s: %d top faces, total power %.6f W%s" % (a.out, nfaces, power, "" if a.gauss else " (Q0 L W = %.6f W)" % (a.q0 * L * W)))


if __name__ == "__main__":
    main()
