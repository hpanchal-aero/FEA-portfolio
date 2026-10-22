#!/usr/bin/env python3
"""Deck-to-deck (static analysis only): replace every MASS element by a nodal *CLOAD of -m*g
in z on its node. The *EQUATION block is kept, so the mass nodes stay tied to the wing.

usage: variant_cload.py <src.inp> <dst.inp> [g]
"""
import re
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
g = float(sys.argv[3]) if len(sys.argv) > 3 else 9.81
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


node_of, mass_of = {}, {}
for kwl, data in blocks:
    if kwl is None:
        continue
    m = mre.search(kwl)
    if not m:
        continue
    k = int(m.group(1))
    rows = [d for d in data if d.strip() and not d.strip().startswith("**")]
    if head(kwl) == "*ELEMENT" and "TYPE=MASS" in kwl.upper().replace(" ", ""):
        t = [x.strip() for x in rows[0].split(",") if x.strip()]
        node_of[k] = int(t[1])
    elif head(kwl) == "*MASS":
        mass_of[k] = float(rows[0].split(",")[0])
if not node_of or set(node_of) != set(mass_of):
    sys.exit("ABORT: MASS elements and *MASS cards do not match")

n_dload = sum(1 for kwl, _ in blocks if kwl and head(kwl) == "*DLOAD")
if n_dload != 1:
    sys.exit(f"ABORT: expected exactly one *DLOAD block, found {n_dload}")

out, n_removed = [], 0
for kwl, data in blocks:
    if kwl is None:
        out.append((kwl, data))
        continue
    h = head(kwl)
    if h in ("*ELEMENT", "*MASS") and mre.search(kwl):
        continue
    if h == "*DLOAD":
        new = []
        for d in data:
            if re.match(r"\s*MS\d+\s*,", d):
                n_removed += 1
            else:
                new.append(d)
        data = new
        cl = [f"{node_of[k]}, 3, {-mass_of[k] * g:.9e}" for k in sorted(node_of)]
        out.append(("*CLOAD", cl))
    out.append((kwl, data))

res = []
for kwl, data in out:
    if kwl is not None:
        res.append(kwl)
    res += data
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(res) + "\n")
M = sum(mass_of.values())
print(f"converted {len(node_of)} MASS elements -> CLOAD; removed {n_removed} gravity lines")
print(f"total mass {M:.1f} kg   weight {M * g:.1f} N (g = {g})")
print(f"wrote {dst}")
