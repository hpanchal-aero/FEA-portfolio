#!/usr/bin/env python3
"""Deck-to-deck feature switch for bisection.

usage: variant_features.py <src.inp> <dst.inp> [--no-eq] [--no-mass-elem] [--no-ms-dload]

  --no-eq         remove the whole *EQUATION block
  --no-mass-elem  remove *ELEMENT/*MASS blocks of sets MS<k> and their *DLOAD lines
  --no-ms-dload   remove only the *DLOAD lines of sets MS<k>
"""
import re
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
flags = set(sys.argv[3:])
bad = flags - {"--no-eq", "--no-mass-elem", "--no-ms-dload"}
if bad:
    sys.exit(f"ABORT: unknown flags {bad}")
if not flags:
    sys.exit("ABORT: no flag given")

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

mre = re.compile(r"ELSET\s*=\s*MS\d+\b", re.I)
removed = dict(eq_blocks=0, elem_blocks=0, mass_blocks=0, dload_lines=0)
res = []
for kwl, data in blocks:
    if kwl is None:
        res += data
        continue
    h = kwl.upper().split(",")[0].strip()
    if "--no-eq" in flags and h == "*EQUATION":
        removed["eq_blocks"] += 1
        continue
    if "--no-mass-elem" in flags and h in ("*ELEMENT", "*MASS") and mre.search(kwl):
        removed["elem_blocks" if h == "*ELEMENT" else "mass_blocks"] += 1
        continue
    if h == "*DLOAD" and flags & {"--no-mass-elem", "--no-ms-dload"}:
        keep = []
        for d in data:
            if re.match(r"\s*MS\d+\s*,", d):
                removed["dload_lines"] += 1
            else:
                keep.append(d)
        data = keep
    res.append(kwl)
    res += data

dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(res) + "\n")
print(f"flags {sorted(flags)}  removed {removed}")
print(f"wrote {dst}")
