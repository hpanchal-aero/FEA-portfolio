#!/usr/bin/env python3
"""Deck-to-deck: replace each point mass (MASS element) by an equivalent T3D2 'truss pair'.

Why: in this CalculiX build a node carrying a MASS element fails as the dependent (first) term of a
multi-term *EQUATION (add_sm_st: coefficient should be 0). A regular element on the node avoids it.

Mass node m keeps its equations. A partner node m' = m + (1,0,0) is tied rigidly to m by three 2-term
equations (dependent = m'). A T3D2 element m--m' (length 1 m) carries the mass rho*A*L = M_k.
Gravity lines on MS<k> are renamed to TM<k>.

Findings (modal runs): with E = 1 Pa the element produces soft spurious modes (eigenvalues ~ E/(rho*A));
with E = 1e9 they vanish and the wing modes are unchanged to 7 digits.

usage: variant_truss_mass.py <src.inp> <dst.inp> [--E <Pa>] [--unit-area]

  default      one material TMAT (rho = 1 kg/m3), section area A_k = M_k m2
  --unit-area  per-mass material TMAT<k> (rho = M_k kg/m3), section area 1 m2
  --E <Pa>     truss Young's modulus (default 1e9)
"""
import re
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
E = 1.0e9
unit_area = False
args = sys.argv[3:]
i = 0
while i < len(args):
    if args[i] == "--E" and i + 1 < len(args):
        E = float(args[i + 1])
        i += 2
    elif args[i] == "--unit-area":
        unit_area = True
        i += 1
    else:
        sys.exit(f"ABORT: bad argument {args[i]}")

lines = src.read_text().splitlines()
blocks, cur = [], [None, []]
for L in lines:
    s = L.strip()
    if s.startswith("*") and not s.startswith("**"):
        blocks.append(cur)
        cur = [L, []]
    else:
        cur[1].append(L)
blocks.append(cur)

mre = re.compile(r"ELSET\s*=\s*MS(\d+)\b", re.I)


def head(kwl):
    return kwl.upper().split(",")[0].strip()


node_of, mass_of, xyz = {}, {}, {}
max_node = max_el = 0
for kwl, data in blocks:
    if kwl is None:
        continue
    h = head(kwl)
    ku = kwl.upper().replace(" ", "")
    rows = [d for d in data if d.strip() and not d.strip().startswith("**")]
    if h == "*NODE":
        for d in rows:
            t = [x.strip() for x in d.split(",")]
            nid = int(t[0])
            max_node = max(max_node, nid)
            xyz[nid] = (float(t[1]), float(t[2]), float(t[3]))
    elif h == "*ELEMENT":
        for d in rows:
            max_el = max(max_el, int(d.split(",")[0]))
        m = mre.search(kwl)
        if m and "TYPE=MASS" in ku:
            t = [x.strip() for x in rows[0].split(",") if x.strip()]
            node_of[int(m.group(1))] = int(t[1])
    elif h == "*MASS":
        m = mre.search(kwl)
        if m:
            mass_of[int(m.group(1))] = float(rows[0].split(",")[0])

ks = sorted(node_of)
if not ks or set(ks) != set(mass_of):
    sys.exit("ABORT: MASS elements and *MASS cards do not match")
if any(node_of[k] not in xyz for k in ks):
    sys.exit("ABORT: a mass node has no coordinates in the deck")

part = {k: max_node + 1 + i for i, k in enumerate(ks)}
eid = {k: max_el + 1 + i for i, k in enumerate(ks)}

new = [("*NODE, NSET=NTRUSS",
        [f"{part[k]}, {xyz[node_of[k]][0] + 1.0:.9e}, {xyz[node_of[k]][1]:.9e}, {xyz[node_of[k]][2]:.9e}"
         for k in ks])]
if unit_area:
    for k in ks:
        new += [(f"*MATERIAL, NAME=TMAT{k}", []),
                ("*ELASTIC", [f"{E:.6e}, 0.0"]),
                ("*DENSITY", [f"{mass_of[k]:.9e}"])]
else:
    new += [("*MATERIAL, NAME=TMAT", []),
            ("*ELASTIC", [f"{E:.6e}, 0.0"]),
            ("*DENSITY", ["1.0"])]
for k in ks:
    mat = f"TMAT{k}" if unit_area else "TMAT"
    area = 1.0 if unit_area else mass_of[k]
    new.append((f"*ELEMENT, TYPE=T3D2, ELSET=TM{k}", [f"{eid[k]}, {node_of[k]}, {part[k]}"]))
    new.append((f"*SOLID SECTION, ELSET=TM{k}, MATERIAL={mat}", [f"{area:.9e}"]))
eq = []
for k in ks:
    for d in (1, 2, 3):
        eq += ["2", f"{part[k]}, {d}, 1.0, {node_of[k]}, {d}, -1.0"]
new.append(("*EQUATION", eq))

out, inserted, n_dl = [], False, 0
for kwl, data in blocks:
    if kwl is None:
        out.append((kwl, data))
        continue
    h = head(kwl)
    if h in ("*ELEMENT", "*MASS") and mre.search(kwl):
        if not inserted:
            out += new
            inserted = True
        continue
    if h == "*DLOAD":
        fixed = []
        for d in data:
            d2 = re.sub(r"^(\s*)MS(\d+)(\s*,)", r"\1TM\2\3", d)
            n_dl += d2 != d
            fixed.append(d2)
        data = fixed
    out.append((kwl, data))
if not inserted:
    sys.exit("ABORT: nothing inserted")

res = []
for kwl, data in out:
    if kwl is not None:
        res.append(kwl)
    res += data
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(res) + "\n")
print(f"masses converted {len(ks)}; total mass {sum(mass_of.values()):.1f} kg; "
      f"E = {E:.3e} Pa; mode: {'unit-area (A = 1 m2, rho = M)' if unit_area else 'A = M m2, rho = 1'}")
print(f"partner nodes {min(part.values())}..{max(part.values())}; truss elements {min(eid.values())}..{max(eid.values())}")
print(f"gravity lines renamed MS->TM: {n_dl}; tie equations added: {3 * len(ks)}")
print(f"wrote {dst}")
