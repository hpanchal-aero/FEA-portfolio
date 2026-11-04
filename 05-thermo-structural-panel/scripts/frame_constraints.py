#!/usr/bin/python3
"""Distributed membrane restraint for the frame support of the panel (Project 05).

Usage: /usr/bin/python3 scripts/frame_constraints.py TAG [--faces x|xy]

For every line of nodes across the thickness on an edge face, one *EQUATION forces the thickness-mean
of the in-plane displacement normal to that face to zero:
    sum_i w_i u_i = 0,   sum_i w_i = 1,
with Simpson weights for lines that have mid-side nodes (2 nz + 1 nodes) and trapezoid weights for lines
with corner nodes only (nz + 1 nodes). The face can rotate and warp freely (the edge rotation is not
restrained), and the reaction is spread over the thickness, in contrast with pinning a single mid-thickness
line, which acts as a point load whose effect grows with mesh refinement.
  x = 0 and x = L faces: u_x (dof 1), one equation per y-line      (always)
  y = 0 and y = W faces: u_y (dof 2), one equation per x-line      (--faces xy)
The dependent term of each equation is the node closest to mid-thickness (needs an even nz).
Reads mesh/TAG.inp (node sets FACE_X0, FACE_X1, FACE_Y0, FACE_Y1) and writes mesh/TAG_frame.inp.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio

TOL = 1e-9


def line_weights(z, h, nz):
    """Weights (sum 1) for the nodes of one through-thickness line, z sorted ascending."""
    n = len(z)
    dz = h / nz
    w = np.zeros(n)
    if n == 2 * nz + 1:
        expect = np.arange(n) * h / (2 * nz)
        for j in range(nz):
            w[2 * j] += dz / 6.0
            w[2 * j + 1] += 4.0 * dz / 6.0
            w[2 * j + 2] += dz / 6.0
    elif n == nz + 1:
        expect = np.arange(n) * h / nz
        for j in range(nz):
            w[j] += dz / 2.0
            w[j + 1] += dz / 2.0
    else:
        raise SystemExit("a line has %d nodes, expected %d or %d" % (n, 2 * nz + 1, nz + 1))
    if np.abs(z - expect).max() > 1e-8:
        raise SystemExit("node heights on a through-thickness line do not match the expected levels")
    return w / h


def face_equations(nodes, ids, key_axis, dof, h, nz):
    groups = {}
    for i in ids:
        c = nodes[i]
        groups.setdefault(int(round(c[key_axis] / TOL)), []).append(i)
    eqs = []
    for key in sorted(groups):
        g = sorted(groups[key], key=lambda i: nodes[i][2])
        z = np.array([nodes[i][2] for i in g])
        w = line_weights(z, h, nz)
        dep = int(np.argmin(np.abs(z - h / 2)))
        order = [dep] + [k for k in range(len(g)) if k != dep]
        if abs(w[dep]) < 1e-6:
            raise SystemExit("dependent node has a zero weight")
        eqs.append([(g[k], dof, float(w[k])) for k in order])
    return eqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tag")
    ap.add_argument("--faces", choices=["x", "xy"], default="xy")
    a = ap.parse_args()
    root = os.path.dirname(HERE)
    mesh = os.path.join(root, "mesh", a.tag + ".inp")
    par = ccxio.read_params(mesh)
    h, nz = par["h"], par["nz"]
    if nz % 2:
        raise SystemExit("nz must be even so that a node exists at mid-thickness")
    nodes = ccxio.read_nodes(mesh)
    jobs = [("FACE_X0", 1, 1), ("FACE_X1", 1, 1)]
    if a.faces == "xy":
        jobs += [("FACE_Y0", 0, 2), ("FACE_Y1", 0, 2)]
    out = os.path.join(root, "mesh", a.tag + "_frame.inp")
    count = 0
    with open(out, "w") as f:
        f.write("** distributed membrane restraint (scripts/frame_constraints.py), mesh %s, faces %s\n" % (a.tag, a.faces))
        f.write("** one equation per through-thickness line: thickness-mean of the in-plane displacement normal to the face = 0\n")
        for name, key_axis, dof in jobs:
            ids = ccxio.read_nset(mesh, name)
            eqs = face_equations(nodes, ids, key_axis, dof, h, nz)
            for terms in eqs:
                f.write("*EQUATION\n%d\n" % len(terms))
                for i in range(0, len(terms), 4):
                    f.write(", ".join("%d, %d, %.12e" % t for t in terms[i:i + 4]) + "\n")
                count += 1
            print("%-8s %4d equations (%d nodes)" % (name, len(eqs), len(ids)))
    print("written %s: %d equations" % (out, count))


if __name__ == "__main__":
    main()
