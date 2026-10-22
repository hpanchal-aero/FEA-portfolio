#!/usr/bin/env python3
"""Deck-to-deck: insert three 2-term translation ties (*EQUATION) before the *STEP.

usage: variant_ties.py <src.inp> <dst.inp> <shellshell|massshell>
  shellshell : tie the two highest-y wing nodes to each other (no point mass involved)
  massshell  : tie the first mass node (NSET=NMASS) to the highest-y wing node
"""
import sys
from pathlib import Path

src, dst, mode = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
if mode not in ("shellshell", "massshell"):
    sys.exit("ABORT: mode must be shellshell or massshell")
lines = src.read_text().splitlines()

wing, mass = {}, []
kw = ""
for L in lines:
    s = L.strip()
    if not s or s.startswith("**"):
        continue
    if s.startswith("*"):
        kw = s.upper().replace(" ", "")
        continue
    if kw.startswith("*NODE") and not kw.startswith(("*NODEFILE", "*NODEPRINT", "*NODEOUTPUT")):
        t = [x.strip() for x in s.split(",")]
        nid = int(t[0])
        if "NSET=NMASS" in kw:
            mass.append(nid)
        else:
            wing[nid] = (float(t[1]), float(t[2]), float(t[3]))

if any(L.strip().upper().startswith("*EQUATION") for L in lines):
    sys.exit("ABORT: source deck already has an *EQUATION block")
idx = [i for i, L in enumerate(lines) if L.strip().upper().startswith("*STEP")]
if len(idx) != 1:
    sys.exit(f"ABORT: expected exactly one *STEP, found {len(idx)}")

by_y = sorted(wing, key=lambda n: -wing[n][1])
if mode == "shellshell":
    n1, n2 = by_y[0], by_y[1]
else:
    if not mass:
        sys.exit("ABORT: no NSET=NMASS node block in the deck")
    n1, n2 = mass[0], by_y[0]

block = ["*EQUATION"]
for d in (1, 2, 3):
    block += ["2", f"{n1}, {d}, 1.0, {n2}, {d}, -1.0"]
out = lines[:idx[0]] + block + lines[idx[0]:]
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(out) + "\n")
print(f"mode {mode}: tie node {n1} <-> node {n2}")
for n in (n1, n2):
    if n in wing:
        print(f"  wing node {n} at {wing[n]}")
    else:
        print(f"  mass node {n}")
print(f"wrote {dst}")
