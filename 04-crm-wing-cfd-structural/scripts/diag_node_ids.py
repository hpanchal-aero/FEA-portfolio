#!/usr/bin/python3
"""Diagnostic: how do the deck node IDs relate to the node IDs in the CalculiX .frd?

usage: diag_node_ids.py <job prefix>      (needs <prefix>.inp and <prefix>.frd)
"""
import sys
from collections import Counter

import numpy as np

prefix = sys.argv[1]

nodes, elems = {}, []
sec_lines = Counter()
section, is_shell = None, False
with open(prefix + ".inp") as f:
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
                section = "*ELEMENT " + typ
            continue
        sec_lines[section] += 1
        tok = [t.strip() for t in line.split(",") if t.strip() != ""]
        if section == "*NODE" and len(tok) >= 4:
            nodes[int(tok[0])] = (float(tok[1]), float(tok[2]), float(tok[3]))
        elif is_shell and len(tok) >= 5:
            elems.append([int(t) for t in tok[1:5]])
print("deck data lines per keyword section (top 12):")
for k, v in sec_lines.most_common(12):
    print(f"  {str(k):28s} {v}")
E = np.array(elems, dtype=np.int64)
used = np.unique(E) if E.size else np.array([], dtype=np.int64)
print(f"deck nodes parsed {len(nodes)}  (ID range {min(nodes) if nodes else None}..{max(nodes) if nodes else None})")
print(f"S4 elements {len(E)}  unique shell node IDs {used.size}  (ID range {used.min() if used.size else None}..{used.max() if used.size else None})")

coords, disp, block = {}, {}, None
with open(prefix + ".frd") as f:
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
fid = np.array(sorted(coords), dtype=np.int64)
print(f"frd coordinate nodes {len(coords)}  (ID range {fid.min()}..{fid.max()})  with DISP {len(disp)}")
print("first 8 frd IDs:", fid[:8].tolist())
dset, cset, uset = set(nodes), set(coords), set(int(n) for n in used)
print(f"shell node IDs found in deck nodes: {len(uset & dset)}  in frd coords: {len(uset & cset)}  in frd disp: {len(uset & set(disp))}")
above = fid[fid > (max(nodes) if nodes else 0)]
print(f"frd IDs above the largest deck ID: {above.size}" + (f"  (min {above.min()}, max {above.max()})" if above.size else ""))

if nodes and len(coords):
    FX = np.array([coords[int(i)] for i in fid])
    for n in used[:5]:
        p = np.array(nodes[int(n)])
        d = np.linalg.norm(FX - p, axis=1)
        near = np.where(d < 0.08)[0]
        info = [(int(fid[k]), round(float(d[k]), 4), round(float(FX[k, 2] - p[2]), 4)) for k in near][:8]
        print(f"deck node {int(n)} at {tuple(round(v, 3) for v in p)}: frd nodes within 0.08 m (id, dist, dz): {info}")
else:
    print("deck or frd empty: nothing to compare")
