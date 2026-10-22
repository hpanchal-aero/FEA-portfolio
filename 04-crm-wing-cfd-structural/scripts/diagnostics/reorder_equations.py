#!/usr/bin/env python3
"""Deck-to-deck: make a wing DOF the dependent (first) term of every *EQUATION.

Finding (variants r1/r2/r3): in this CalculiX build the first term of an *EQUATION is the dependent
DOF, and a MASS-element node as the dependent DOF of a multi-term equation fails with
'add_sm_st: coefficient should be 0'. With a wing DOF first it runs.
For each equation the wing term with the largest |coefficient| that is (a) not in any SPC set and
(b) not already dependent in another equation becomes the first term. Coefficients are unchanged.

usage: reorder_equations.py <src.inp> <dst.inp>
"""
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
lines = src.read_text().splitlines()

# ---- SPC node set from *NSET and *BOUNDARY ----
nsets, bnd = {}, []
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
    if kw.startswith("*NSET") and cur:
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
if unresolved:
    sys.exit(f"ABORT: unresolved *BOUNDARY names {sorted(unresolved)}")
if not spc:
    sys.exit("ABORT: no SPC nodes found")
print(f"SPC nodes found: {len(spc)}")


def is_kw(s):
    return s.startswith("*") and not s.startswith("**")


# ---- rewrite the equations ----
out, i = [], 0
used = set()
indep = set()
n_eq = 0
chosen = []
while i < len(lines):
    L = lines[i]
    if L.strip().upper().startswith("*EQUATION"):
        out.append(L)
        i += 1
        while i < len(lines) and not is_kw(lines[i].strip()):
            s2 = lines[i].strip()
            if not s2 or s2.startswith("**"):
                i += 1
                continue
            n = int(s2.split(",")[0])
            vals = []
            while len(vals) < 3 * n:
                i += 1
                vals += [x.strip() for x in lines[i].split(",") if x.strip()]
            i += 1
            terms = [(int(vals[3 * j]), int(vals[3 * j + 1]), vals[3 * j + 2]) for j in range(n)]
            cand = [(abs(float(c)), j) for j, (nd, dof, c) in enumerate(terms)
                    if j > 0 and nd not in spc and (nd, dof) not in used]
            if not cand:
                sys.exit(f"ABORT: equation {n_eq + 1} has no usable wing DOF for the dependent term")
            _, jb = max(cand)
            nd, dof, c = terms[jb]
            used.add((nd, dof))
            chosen.append(abs(float(c)))
            new = [terms[jb]] + [t for j, t in enumerate(terms) if j != jb]
            for t in new[1:]:
                indep.add((t[0], t[1]))
            out.append(str(len(new)))
            for j in range(0, len(new), 4):
                out.append(", ".join(f"{a}, {b}, {c2}" for a, b, c2 in new[j:j + 4]))
            n_eq += 1
        continue
    out.append(L)
    i += 1

if n_eq == 0:
    sys.exit("ABORT: no equations found")
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(out) + "\n")
print(f"equations {n_eq}; unique dependent wing DOFs {len(used)}")
print(f"|coef| of chosen dependent term: min {min(chosen):.3e}  max {max(chosen):.3e}")
print(f"dependent DOFs that also appear as independent terms elsewhere (cascade cases): {len(used & indep)}")
print(f"wrote {dst}")
