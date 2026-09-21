"""
Diagnostic only (NOT the locked Config A boundary condition).

Copies a Config A .inp and replaces the fixed-base BC (NBOTTOM 1-3 = 0) with a
roller base: Uz = 0 on all bottom nodes, plus the minimum extra constraints to
remove rigid-body motion:
    node 1 (-50,-50,0): Ux, Uy = 0
    node 5 ( 50,-50,0): Uy = 0
Purpose: test whether the ~1.7% axial-stiffness offset and the rail-region
stress deficit come from lateral restraint at the fixed base.

Usage: /usr/bin/python3 scripts/make_roller_bc_case.py <src.inp> <dst.inp>
"""
import sys

OLD = "*BOUNDARY\nNBOTTOM, 1, 3, 0.0\n"
NEW = ("*BOUNDARY\nNBOTTOM, 3, 3, 0.0\n"
       "1, 1, 2, 0.0\n"
       "5, 2, 2, 0.0\n")

src, dst = sys.argv[1], sys.argv[2]
text = open(src).read()
if text.count(OLD) != 1:
    print(f"ERROR: expected exactly 1 occurrence of the fixed-base BC, found {text.count(OLD)}")
    sys.exit(1)

# confirm nodes 1 and 5 are the expected corner nodes before constraining them
want = {1: (-50.0, -50.0, 0.0), 5: (50.0, -50.0, 0.0)}
in_nodes = False
for line in text.splitlines():
    s = line.strip()
    if s.startswith("*"):
        in_nodes = s.upper().startswith("*NODE") and not s.upper().startswith("*NODE FILE")
        continue
    if in_nodes and s:
        p = s.split(",")
        nid = int(p[0])
        if nid in want:
            c = tuple(round(float(v), 6) for v in p[1:4])
            if c != want[nid]:
                print(f"ERROR: node {nid} is at {c}, expected {want[nid]}")
                sys.exit(1)
            print(f"node {nid} confirmed at {c}")
        if nid > 5:
            break

open(dst, "w").write(text.replace(OLD, NEW))
print(f"Wrote {dst} (roller base, rigid-body modes removed)")
