#!/usr/bin/env python3
"""Deck-to-deck variant: keep only the first N point masses (mass node, MASS element,
*MASS card, its 3 equations, its GRAV line). N=0 removes all point masses.

usage: variant_masses.py <src.inp> <dst.inp> <N>
"""
import re
import sys
from pathlib import Path

src, dst, N = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
lines = src.read_text().splitlines()

# split into (keyword line, data lines) blocks
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


# pass 1: which mass nodes are kept
kept_nodes = None
for kwl, data in blocks:
    if kwl and head(kwl) == "*NODE" and "NSET=NMASS" in kwl.upper().replace(" ", ""):
        rows = [d for d in data if d.strip() and not d.strip().startswith("**")]
        kept_nodes = {int(d.split(",")[0]) for d in rows[:N]}
        n_mass_nodes = len(rows)
if kept_nodes is None:
    sys.exit("ABORT: no '*NODE, NSET=NMASS' block found")

stats = dict(node=0, element=0, mass=0, dload=0, eq_kept=0, eq_dropped=0)
out = []
for kwl, data in blocks:
    if kwl is None:
        out.append((kwl, data))
        continue
    h = head(kwl)
    m = mre.search(kwl)
    if h in ("*ELEMENT", "*MASS") and m:
        if int(m.group(1)) >= N:
            stats["element" if h == "*ELEMENT" else "mass"] += 1
            continue
    elif h == "*NODE" and "NSET=NMASS" in kwl.upper().replace(" ", ""):
        rows = [d for d in data if d.strip() and not d.strip().startswith("**")]
        stats["node"] = len(rows) - len(rows[:N])
        if N == 0:
            continue
        data = rows[:N]
    elif h == "*DLOAD":
        new = []
        for d in data:
            mm = re.match(r"\s*MS(\d+)\s*,", d)
            if mm and int(mm.group(1)) >= N:
                stats["dload"] += 1
                continue
            new.append(d)
        data = new
    elif h == "*EQUATION":
        new, i = [], 0
        while i < len(data):
            s = data[i].strip()
            if not s or s.startswith("**"):
                i += 1
                continue
            n = int(s.split(",")[0])
            body, vals = [], []
            while len(vals) < 3 * n:
                i += 1
                body.append(data[i])
                vals += [x.strip() for x in data[i].split(",") if x.strip()]
            if len(vals) != 3 * n:
                sys.exit("ABORT: malformed equation block")
            if int(vals[0]) in kept_nodes:
                new += [data[i - len(body)]] + body
                stats["eq_kept"] += 1
            else:
                stats["eq_dropped"] += 1
            i += 1
        if not new:
            continue
        data = new
    out.append((kwl, data))

res = []
for kwl, data in out:
    if kwl is not None:
        res.append(kwl)
    res += data
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(res) + "\n")
print(f"N={N}: mass nodes kept {len(kept_nodes)} of {n_mass_nodes}; removed: {stats}")
print(f"wrote {dst}")
