#!/usr/bin/python3
"""Check of the frd -> deck displacement mapping used by contour_displacement.py and the video scripts.

usage: check_mapping_vs_dat.py <job prefix>      (needs <prefix>.inp, <prefix>.frd and <prefix>.dat)

For every deck shell node the displacement from the mapping (mean over the frd nodes within RMAX of it) is compared with the
displacement CalculiX prints in <prefix>.dat for the same deck node ID.  This does not depend on the frd nodes that are not
assigned to any shell node.
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


def read_dat_disp(path):
    disp, on = {}, False
    for line in open(path):
        if "displacements (vx,vy,vz)" in line:
            on = True
            continue
        if not on:
            continue
        t = line.split()
        if not t:
            if disp:
                break
            continue
        try:
            disp[int(t[0])] = [float(v) for v in t[1:4]]
        except ValueError:
            break
    return disp


prefix = sys.argv[1]
nodes, E = read_deck(prefix + ".inp")
coords, disp = read_frd(prefix + ".frd")
dat = read_dat_disp(prefix + ".dat")
used = np.unique(E)
P0 = np.array([nodes[int(n)] for n in used], dtype=float)
fid = np.array(sorted(coords), dtype=np.int64)
FX = np.array([coords[int(i)] for i in fid], dtype=float)
FU = np.array([disp[int(i)] for i in fid], dtype=float)
d, k = cKDTree(P0).query(FX)
near = d <= RMAX
cnt = np.bincount(k[near], minlength=len(P0))
if (cnt == 0).any():
    sys.exit("ABORT: deck shell nodes without an assigned frd node")
U = np.column_stack([np.bincount(k[near], weights=FU[near, c], minlength=len(P0)) / cnt for c in range(3)])

common = np.array([int(n) for n in used if int(n) in dat])
print(f"{prefix}")
print(f"deck shell nodes {len(used)}; nodes with displacement in the .dat {len(dat)}; common {len(common)}")
if len(common) == 0:
    sys.exit("ABORT: no common nodes (check the .dat format)")
idx = np.searchsorted(used, common)
D = np.array([dat[int(n)] for n in common])
diff = U[idx] - D
mag = np.linalg.norm(diff, axis=1)
j = int(np.argmax(mag))
umax = float(np.linalg.norm(D, axis=1).max())
print(f"max |u| in the .dat {umax:.4f} m")
print(f"|mapped - .dat|: max {mag.max():.2e} m ({100.0 * mag.max() / umax:.3f} % of max |u|)  p99 {np.percentile(mag, 99):.2e}  "
      f"RMS {np.sqrt(np.mean(mag ** 2)):.2e} m")
print(f"max component differences ux {np.abs(diff[:, 0]).max():.2e}  uy {np.abs(diff[:, 1]).max():.2e}  uz {np.abs(diff[:, 2]).max():.2e} m")
print(f"worst node: deck ID {int(common[j])} at x = ({P0[idx[j], 0]:.3f}, {P0[idx[j], 1]:.3f}, {P0[idx[j], 2]:.3f}) m, "
      f"mapped uz {U[idx[j], 2]:.5f}, .dat uz {D[j, 2]:.5f}")
print(f"nodes with difference above 1e-3 m: {int((mag > 1e-3).sum())}; above 1e-2 m: {int((mag > 1e-2).sum())}")
