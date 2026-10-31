#!/usr/bin/python3
"""Diagnostic: which frd nodes are NOT assigned to a deck shell node (distance > RMAX), and does it matter?

usage: check_frd_mapping.py <job prefix>      (needs <prefix>.inp and <prefix>.frd)

Same reading and nearest-node assignment as contour_displacement.py / make_bending_video.py.  Reports for the unassigned
frd nodes: how many coincide with deck nodes that are not shell nodes (point-mass truss nodes), the distance statistics of
the rest, their bounding box, and how much their displacement differs from the mean displacement of the nearest shell node.
"""
import sys

import numpy as np
from scipy.spatial import cKDTree

RMAX = 0.05


def read_deck(path):
    nodes, elems = {}, []
    section, is_shell = None, False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                parts = [p.strip() for p in line.split(",")]
                section = parts[0].upper()
                is_shell = False
                if section == "*ELEMENT":
                    typ = ""
                    for p in parts[1:]:
                        if p.upper().startswith("TYPE="):
                            typ = p.split("=", 1)[1].upper()
                    is_shell = typ.startswith("S4")
                continue
            tok = [t.strip() for t in line.split(",") if t.strip() != ""]
            if section == "*NODE" and len(tok) >= 4:
                nodes[int(tok[0])] = (float(tok[1]), float(tok[2]), float(tok[3]))
            elif is_shell and len(tok) >= 5:
                elems.append([int(t) for t in tok[1:5]])
    return nodes, np.array(elems, dtype=np.int64)


def read_frd(path):
    coords, disp, block = {}, {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                block = "coord"
            elif line.startswith(" -4  DISP"):
                block = "disp"
            elif line.startswith(" -4  "):
                block = None
            elif line.startswith(" -3"):
                block = None
            elif line.startswith(" -1") and block:
                node = int(line[3:13])
                vals = [float(line[13 + 12 * i:25 + 12 * i]) for i in range(3)]
                (coords if block == "coord" else disp)[node] = vals
    return coords, disp


prefix = sys.argv[1]
nodes, E = read_deck(prefix + ".inp")
coords, disp = read_frd(prefix + ".frd")
used = np.unique(E)
P0 = np.array([nodes[int(n)] for n in used], dtype=float)
allid = np.array(sorted(nodes))
PA = np.array([nodes[int(n)] for n in allid], dtype=float)
fid = np.array(sorted(coords), dtype=np.int64)
FX = np.array([coords[int(i)] for i in fid], dtype=float)
FU = np.array([disp[int(i)] for i in fid], dtype=float)

tree = cKDTree(P0)
d, k = tree.query(FX)
near = d <= RMAX
cnt = np.bincount(k[near], minlength=len(P0))
U = np.column_stack([np.bincount(k[near], weights=FU[near, c], minlength=len(P0)) / np.maximum(cnt, 1) for c in range(3)])
print(f"{prefix}")
print(f"frd nodes {len(fid)}; deck nodes {len(allid)} = shell {len(used)} + other {len(allid) - len(used)}")
print(f"assigned (<= {RMAX} m) {int(near.sum())}; unassigned {int((~near).sum())}")
far = ~near
if far.any():
    d2 = cKDTree(PA).query(FX[far])[0]
    coin = d2 < 2e-3
    print(f"distance of the unassigned frd nodes to the nearest deck node: min {d2.min():.2e}  median {np.median(d2):.2e}  max {d2.max():.2e} m;  below 1e-3 m: {int((d2 < 1e-3).sum())}")
    print(f"unassigned frd nodes coinciding with a deck node (distance < 2e-3 m; the frd coordinate resolution is about 1e-4 m): {int(coin.sum())}")
    print(f"unassigned frd nodes not coinciding with any deck node: {int((~coin).sum())}")
    for name, m in (("coinciding", coin), ("not coinciding", ~coin)):
        if not m.any():
            continue
        idx = np.where(far)[0][m]
        dd = d[idx]
        dev = np.linalg.norm(FU[idx] - U[k[idx]], axis=1)
        box = f"x {FX[idx, 0].min():.2f}..{FX[idx, 0].max():.2f}  y {FX[idx, 1].min():.2f}..{FX[idx, 1].max():.2f}  z {FX[idx, 2].min():.2f}..{FX[idx, 2].max():.2f}"
        print(f"  {name}: {len(idx)} nodes, distance to nearest shell node min {dd.min():.3f}  median {np.median(dd):.3f}  max {dd.max():.3f} m")
        print(f"     bounding box {box}")
        print(f"     |u - mean u of the nearest shell node|: median {np.median(dev):.2e}  max {dev.max():.2e} m")
    idx = np.where(far)[0][~coin]
    if len(idx):
        h, _ = np.histogram(d[idx], bins=[RMAX, 0.1, 0.2, 0.5, 1.0, 1e9])
        print("  non-coinciding distance bins [0.05-0.1, 0.1-0.2, 0.2-0.5, 0.5-1, >1 m]:", h.tolist())
