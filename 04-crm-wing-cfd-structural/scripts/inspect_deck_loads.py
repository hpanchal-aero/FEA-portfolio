#!/usr/bin/env python3
"""Read-only inspection of the loads and boundary conditions of a CalculiX deck.
usage: inspect_deck_loads.py <deck.inp>
Reports: node blocks, *CLOAD totals (force, moment about the origin, centre of pressure),
spanwise distribution of the *CLOAD forces, and the *BOUNDARY constraint classes with coordinates.
"""
import sys
import collections
import numpy as np


def parse(path):
    nodes = {}
    block_of = {}
    cload = collections.defaultdict(lambda: [0.0] * 6)
    cload_names = collections.Counter()
    boundary = []
    section = None
    nset = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                parts = [p.strip() for p in line.split(",")]
                section = parts[0].upper()
                nset = None
                for p in parts[1:]:
                    if p.upper().startswith("NSET="):
                        nset = p.split("=", 1)[1]
                continue
            tok = [t.strip() for t in line.split(",")]
            if section == "*NODE":
                nid = int(tok[0])
                nodes[nid] = (float(tok[1]), float(tok[2]), float(tok[3]))
                block_of[nid] = nset
            elif section == "*CLOAD":
                try:
                    nid = int(tok[0])
                except ValueError:
                    cload_names[tok[0]] += 1
                    continue
                dof = int(tok[1])
                if 1 <= dof <= 6:
                    cload[nid][dof - 1] += float(tok[2])
            elif section == "*BOUNDARY":
                d1 = int(tok[1])
                d2 = int(tok[2]) if len(tok) > 2 and tok[2] != "" else d1
                try:
                    key = int(tok[0])
                except ValueError:
                    key = tok[0]
                boundary.append((key, d1, d2))
    return nodes, block_of, cload, cload_names, boundary


def main(path):
    nodes, block_of, cload, cload_names, boundary = parse(path)
    print("deck:", path)
    cnt = collections.Counter(block_of.values())
    print("node blocks (NSET: count):", dict(cnt), " total", len(nodes))

    ids = sorted(cload)
    ok = [i for i in ids if i in nodes]
    miss = [i for i in ids if i not in nodes]
    print("\n=== *CLOAD ===")
    print("loaded nodes:", len(ids), " with coordinates:", len(ok), " without:", len(miss), " named sets:", dict(cload_names))
    F = np.array([cload[i][:3] for i in ok])
    X = np.array([nodes[i] for i in ok])
    tot = F.sum(axis=0)
    M = np.cross(X, F).sum(axis=0)
    print("total force  Fx %.1f  Fy %.1f  Fz %.1f N   (|F| = %.1f N)" % (tot[0], tot[1], tot[2], np.linalg.norm(tot)))
    print("moment about origin  Mx %.4e  My %.4e  Mz %.4e N m" % tuple(M))
    if abs(tot[2]) > 0:
        print("centre of pressure (Fz-weighted): x %.3f  y %.3f m" % ((X[:, 0] * F[:, 2]).sum() / tot[2], (X[:, 1] * F[:, 2]).sum() / tot[2]))
    print("any non-zero moment DOFs (4-6):", any(abs(cload[i][3]) + abs(cload[i][4]) + abs(cload[i][5]) > 0 for i in ids))
    print("loaded-node bbox  x %.3f..%.3f  y %.3f..%.3f  z %.3f..%.3f" % (X[:, 0].min(), X[:, 0].max(), X[:, 1].min(), X[:, 1].max(), X[:, 2].min(), X[:, 2].max()))
    nall = np.array([v for i, v in nodes.items() if block_of[i] == "NALL"])
    print("NALL bbox         x %.3f..%.3f  y %.3f..%.3f  z %.3f..%.3f  (%d nodes)" % (nall[:, 0].min(), nall[:, 0].max(), nall[:, 1].min(), nall[:, 1].max(), nall[:, 2].min(), nall[:, 2].max(), len(nall)))

    ymax = X[:, 1].max()
    edges = np.arange(0.0, np.ceil(ymax / 3.0) * 3.0 + 3.0, 3.0)
    print("\nspanwise distribution of the *CLOAD forces (bins in y)")
    print("   y0    y1   nodes     Fx[kN]    Fy[kN]    Fz[kN]")
    for a, b in zip(edges[:-1], edges[1:]):
        m = (X[:, 1] >= a) & (X[:, 1] < b)
        if m.any():
            s = F[m].sum(axis=0) / 1000.0
            print("%5.1f %5.1f %7d %10.2f %10.2f %10.2f" % (a, b, m.sum(), s[0], s[1], s[2]))

    print("\n=== *BOUNDARY ===")
    print("boundary lines:", len(boundary))
    groups = collections.defaultdict(list)
    for key, d1, d2 in boundary:
        groups[(d1, d2)].append(key)
    for (d1, d2), keys in sorted(groups.items()):
        num = [k for k in keys if isinstance(k, int)]
        names = [k for k in keys if not isinstance(k, int)]
        line = "DOF %d-%d: %d lines (%d unique numeric nodes, %d named)" % (d1, d2, len(keys), len(set(num)), len(names))
        pts = np.array([nodes[k] for k in set(num) if k in nodes])
        if len(pts):
            line += "  x %.2f..%.2f  y %.3f..%.3f  z %.2f..%.2f" % (pts[:, 0].min(), pts[:, 0].max(), pts[:, 1].min(), pts[:, 1].max(), pts[:, 2].min(), pts[:, 2].max())
        print(line)
        if names:
            print("   named:", sorted(set(names))[:10])
    allnum = set(k for k, _, _ in boundary if isinstance(k, int))
    print("unique constrained nodes overall:", len(allnum))
    dofcount = collections.Counter()
    for k, d1, d2 in boundary:
        for d in range(d1, d2 + 1):
            dofcount[d] += 1
    print("constraints per DOF (line counts):", dict(sorted(dofcount.items())))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
