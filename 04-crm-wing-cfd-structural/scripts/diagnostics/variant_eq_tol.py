#!/usr/bin/env python3
"""Deck-to-deck variant: drop *EQUATION terms with |coef| < tol (never the dependent first term).

usage: variant_eq_tol.py <src.inp> <dst.inp> <tol>
"""
import sys
from pathlib import Path

src, dst, tol = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
lines = src.read_text().splitlines()

out, i = [], 0
in_eq = False
n_all = n_drop = n_eq = 0
worst_sum = 0.0
while i < len(lines):
    L = lines[i]
    s = L.strip()
    if s.startswith("*") and not s.startswith("**"):
        in_eq = s.split(",")[0].upper() == "*EQUATION"
        out.append(L)
        i += 1
        continue
    if in_eq and s and not s.startswith("**"):
        n = int(s.split(",")[0])
        vals = []
        while len(vals) < 3 * n:
            i += 1
            vals += [x.strip() for x in lines[i].split(",") if x.strip()]
        if len(vals) != 3 * n:
            sys.exit(f"ABORT: malformed equation near line {i + 1}")
        terms = [(vals[3 * k], vals[3 * k + 1], vals[3 * k + 2]) for k in range(n)]
        keep = [terms[0]] + [t for t in terms[1:] if abs(float(t[2])) >= tol]
        n_all += n
        n_drop += n - len(keep)
        n_eq += 1
        worst_sum = max(worst_sum, abs(sum(float(t[2]) for t in keep)))
        out.append(str(len(keep)))
        for j in range(0, len(keep), 4):
            out.append(", ".join(f"{a}, {b}, {c}" for a, b, c in keep[j:j + 4]))
        i += 1
        continue
    out.append(L)
    i += 1

if n_eq == 0:
    sys.exit("ABORT: no equations found")
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(out) + "\n")
print(f"equations {n_eq}  terms {n_all} -> {n_all - n_drop}  dropped {n_drop}  "
      f"max |sum of coefs| after drop {worst_sum:.3e}")
print(f"wrote {dst}")
