#!/usr/bin/env python3
"""Read-only audit of the *EQUATION and MASS-element blocks of a CalculiX deck.

usage: check_equations.py <deck.inp>
"""
import re
import sys
from collections import Counter
from pathlib import Path

lines = Path(sys.argv[1]).read_text().splitlines()

nodes = set()
mass_elems = {}      # element id -> node id
eqs = []             # list of [(node, dof, coef), ...]
kw, etype = "", ""
i = 0
while i < len(lines):
    L = lines[i].strip()
    if not L or L.startswith("**"):
        i += 1
        continue
    if L.startswith("*"):
        kw = L.split(",")[0].upper()
        m = re.search(r"TYPE\s*=\s*(\w+)", L, re.I)
        etype = m.group(1).upper() if m else ""
        i += 1
        continue
    if kw == "*NODE":
        nodes.add(int(L.split(",")[0]))
    elif kw == "*ELEMENT" and etype == "MASS":
        t = [x for x in L.split(",") if x.strip()]
        mass_elems[int(t[0])] = int(t[1])
    elif kw == "*EQUATION":
        n = int(L.split(",")[0])
        vals = []
        while len(vals) < 3 * n:
            i += 1
            vals += [x.strip() for x in lines[i].split(",") if x.strip()]
        if len(vals) != 3 * n:
            sys.exit(f"ABORT: equation near line {i + 1}: expected {3 * n} values, got {len(vals)}")
        eqs.append([(int(vals[3 * k]), int(vals[3 * k + 1]), float(vals[3 * k + 2])) for k in range(n)])
    i += 1

print(f"nodes defined          : {len(nodes)}")
print(f"MASS elements          : {len(mass_elems)}")
print(f"equations              : {len(eqs)}")
print(f"total terms            : {sum(len(e) for e in eqs)}")
tc = Counter(len(e) for e in eqs)
print(f"terms/equation min/max : {min(tc)} / {max(tc)}")
print(f"DOF histogram          : {dict(sorted(Counter(d for e in eqs for _, d, _ in e).items()))}")

dup = sum(1 for e in eqs if len({(n, d) for n, d, _ in e}) != len(e))
print(f"equations with repeated node/DOF pair : {dup}")

missing = {n for e in eqs for n, _, _ in e if n not in nodes}
print(f"equation nodes missing from *NODE      : {len(missing)}")

dep = Counter(e[0][0] for e in eqs)
mnodes = set(mass_elems.values())
print(f"distinct dependent (first-term) nodes  : {len(dep)}   all are mass nodes: {set(dep) <= mnodes}")
print(f"mass nodes without any equation        : {len(mnodes - set(dep))}")
print(f"mass nodes that also appear as non-first term: "
      f"{len({n for e in eqs for n, _, _ in e[1:]} & mnodes)}")
print(f"dependent terms with coef != 1.0       : {sum(1 for e in eqs if e[0][2] != 1.0)}")
print(f"MASS elements pointing to missing node : {sum(1 for n in mass_elems.values() if n not in nodes)}")

allc = [abs(c) for e in eqs for _, _, c in e]
print(f"|coef| min / max       : {min(allc):.3e} / {max(allc):.3e}")
print(f"terms with |coef|<1e-9 : {sum(1 for c in allc if c < 1e-9)}")
print(f"max |sum of coefs| per equation (translation rigid-body check): "
      f"{max(abs(sum(c for _, _, c in e)) for e in eqs):.3e}")
