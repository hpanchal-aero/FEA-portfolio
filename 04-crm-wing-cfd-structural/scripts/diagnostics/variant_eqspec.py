#!/usr/bin/env python3
"""Deck-to-deck: insert user-specified *EQUATION entries before the *STEP.

usage: variant_eqspec.py <src.inp> <dst.inp> --eq "node,dof,coef;node,dof,coef;..." [--eq ...]

node tokens: an integer id, MASS (first node of NSET=NMASS), or TIP<k> (k-th highest-y wing node).
Informational: reports whether each wing node used is in an SPC set found in *BOUNDARY.
"""
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
args = sys.argv[3:]
specs = []
i = 0
while i < len(args):
    if args[i] == "--eq" and i + 1 < len(args):
        specs.append(args[i + 1])
        i += 2
    else:
        sys.exit(f"ABORT: bad argument {args[i]}")
if not specs:
    sys.exit("ABORT: no --eq given")

lines = src.read_text().splitlines()
if any(L.strip().upper().startswith("*EQUATION") for L in lines):
    sys.exit("ABORT: source deck already has an *EQUATION block")
idx = [k for k, L in enumerate(lines) if L.strip().upper().startswith("*STEP")]
if len(idx) != 1:
    sys.exit(f"ABORT: expected exactly one *STEP, found {len(idx)}")

wing, mass, nsets, bnd = {}, [], {}, []
kw, cur, gen = "", None, False
for L in lines:
    s = L.strip()
    if not s or s.startswith("**"):
        continue
    if s.startswith("*"):
        kw = s.upper().replace(" ", "")
        cur, gen = None, False
        if kw.startswith("*NSET"):
            for p in s.split(",")[1:]:
                p = p.strip()
                if p.upper().startswith("NSET="):
                    cur = p[5:].strip().upper()
                    nsets.setdefault(cur, set())
                if p.upper() == "GENERATE":
                    gen = True
        continue
    if kw.startswith("*NODE") and not kw.startswith(("*NODEFILE", "*NODEPRINT", "*NODEOUTPUT")):
        t = [x.strip() for x in s.split(",")]
        nid = int(t[0])
        if "NSET=NMASS" in kw:
            mass.append(nid)
        else:
            wing[nid] = (float(t[1]), float(t[2]), float(t[3]))
    elif kw.startswith("*NSET") and cur:
        try:
            v = [int(x) for x in s.split(",") if x.strip()]
        except ValueError:
            continue
        if gen:
            nsets[cur].update(range(v[0], v[1] + 1, v[2] if len(v) > 2 else 1))
        else:
            nsets[cur].update(v)
    elif kw.startswith("*BOUNDARY"):
        bnd.append(s.split(",")[0].strip())

spc, unresolved = set(), set()
for t in bnd:
    if t.isdigit():
        spc.add(int(t))
    elif t.upper() in nsets:
        spc |= nsets[t.upper()]
    else:
        unresolved.add(t)

by_y = sorted(wing, key=lambda n: -wing[n][1])


def resolve(tok):
    tok = tok.strip()
    if tok == "MASS":
        if not mass:
            sys.exit("ABORT: no NSET=NMASS node in deck")
        return mass[0]
    if tok.startswith("TIP"):
        return by_y[int(tok[3:]) - 1]
    n = int(tok)
    if n not in wing and n not in mass:
        sys.exit(f"ABORT: node {n} not in deck")
    return n


block = ["*EQUATION"]
for e, spec in enumerate(specs, 1):
    terms = []
    for item in spec.split(";"):
        a, d, c = [x.strip() for x in item.split(",")]
        terms.append((resolve(a), int(d), float(c)))
    print(f"equation {e}: {len(terms)} terms")
    for n, d, c in terms:
        if n in wing:
            flag = "SPC" if n in spc else ("not in SPC sets" + (f" (unresolved BC names: {sorted(unresolved)})" if unresolved else ""))
            print(f"  wing node {n} dof {d} coef {c:+.4f}  y={wing[n][1]:.3f}  [{flag}]")
        else:
            print(f"  mass node {n} dof {d} coef {c:+.4f}")
    block.append(str(len(terms)))
    for j in range(0, len(terms), 4):
        block.append(", ".join(f"{n}, {d}, {c:.12e}" for n, d, c in terms[j:j + 4]))

out = lines[:idx[0]] + block + lines[idx[0]:]
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(out) + "\n")
print(f"wrote {dst}")
