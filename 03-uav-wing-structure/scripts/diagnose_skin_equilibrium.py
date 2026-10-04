#!/usr/bin/env python3
"""
Where is the reaction force?  Sums FORC (RF) over groups of .frd nodes:
  ORIG (id <= number of *NODE entries in the .inp), EXPANDED (id > that), ALL;
and over the root plane (z = z_min) vs elsewhere, to test whether the missing ~4 N sits on
CalculiX's expanded shell nodes.

Usage (system Python):
  /usr/bin/python3 scripts/diagnose_skin_equilibrium.py <job_path_without_ext>
"""
import re
import sys
import numpy as np

job = sys.argv[1]
inp, frd = job + ".inp", job + ".frd"

# number of original nodes and applied loads from the .inp
n0, mode, lift, ids_root = 0, None, 0.0, []
with open(inp) as f:
    for line in f:
        s = line.strip()
        if not s:
            continue
        if s.startswith("*"):
            u = s.upper()
            mode = "node" if (u.startswith("*NODE") and not u.startswith("*NODE FILE")) else \
                   ("cload" if u.startswith("*CLOAD") else None)
            continue
        if mode == "node":
            n0 += 1
        elif mode == "cload":
            p = [t.strip() for t in s.split(",") if t.strip()]
            if int(p[1]) == 2:
                lift += float(p[2])
print(f"Original nodes: {n0}; applied lift (CLOAD Y): {lift:.6f} N")

# frd node coordinates
ids, xyz = [], []
in_nodes = False
with open(frd) as f:
    for line in f:
        if line.startswith("    2C"):
            in_nodes = True
            continue
        if in_nodes:
            if line.startswith(" -3"):
                break
            if line.startswith(" -1"):
                ids.append(int(line[3:13]))
                xyz.append([float(line[13 + 12 * k:25 + 12 * k]) for k in range(3)])
ids = np.array(ids)
xyz = np.array(xyz)
print(f"frd nodes: {len(ids)} (id range {ids.min()}..{ids.max()}); "
      f"original {int((ids <= n0).sum())}, expanded {int((ids > n0).sum())}")

# RF block
rid, rf = [], []
in_blk = False
with open(frd) as f:
    for line in f:
        if line.startswith(" -4"):
            in_blk = (line[5:13].strip() == "FORC")
            continue
        if in_blk:
            if line.startswith(" -3"):
                in_blk = False
            elif line.startswith(" -1"):
                rid.append(int(line[3:13]))
                rf.append([float(line[13 + 12 * k:25 + 12 * k]) for k in range(3)])
rid = np.array(rid)
rf = np.array(rf)
order = np.argsort(rid)
rid, rf = rid[order], rf[order]
o2 = np.argsort(ids)
ids, xyz = ids[o2], xyz[o2]
if not np.array_equal(ids, rid):
    print("ERROR: node lists differ between coordinate and FORC blocks"); sys.exit(1)

zmin = xyz[:, 2].min()
plane = np.abs(xyz[:, 2] - zmin) < 1e-3
orig = ids <= n0
print(f"z_min of all frd nodes: {zmin:.6f}")


def row(label, m):
    s = rf[m].sum(axis=0)
    nz = int((np.abs(rf[m]).max(axis=1) > 1e-9).sum())
    print(f"  {label:34s}: n={int(m.sum()):7d} nonzero-RF={nz:6d}  "
          f"sumRF X={s[0]:+10.5f} Y={s[1]:+10.5f} Z={s[2]:+10.5f} N")


print("\nSum of FORC by group:")
row("ALL nodes", np.ones(len(ids), bool))
row("original IDs", orig)
row("expanded IDs", ~orig)
row("root plane, original", plane & orig)
row("root plane, expanded", plane & ~orig)
row("root plane, all", plane)
row("NOT on root plane, original", ~plane & orig)
row("NOT on root plane, expanded", ~plane & ~orig)

tot = rf.sum(axis=0)
print(f"\nExpected total reaction Y (naive, from extractor): -8.458595 N; "
      f"sum over ALL nodes: {tot[1]:+.6f} N  -> diff {tot[1] + 8.458595:+.6f} N")
print("(Addendum-A-type gravity share on fixed nodes, ~0.05 N, is NOT in RF and is expected to "
      "show up as a diff of that order.)")
big = np.argsort(np.abs(rf[~plane].sum(axis=1)))[::-1][:5]
print("\nLargest |RF| away from the root plane (should be ~0 unless a constraint exists there):")
sub_ids, sub_xyz, sub_rf = ids[~plane], xyz[~plane], rf[~plane]
for k in big:
    print(f"  id {sub_ids[k]:7d} at ({sub_xyz[k,0]:8.3f},{sub_xyz[k,1]:8.3f},{sub_xyz[k,2]:9.3f}) "
          f"RF=({sub_rf[k,0]:+.3e},{sub_rf[k,1]:+.3e},{sub_rf[k,2]:+.3e})")
