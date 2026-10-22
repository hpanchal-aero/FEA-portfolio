#!/usr/bin/env python3
"""Deck-to-deck: keep a subset of the *EQUATION entries, optionally dropping cross-direction terms.

usage: variant_eqsubset.py <src.inp> <dst.inp> [--keep 1,2,3] [--same-dof-only]

  --keep LIST        keep only these equations (1-based, in deck order)
  --same-dof-only    keep only terms whose DOF equals the dependent (first) term's DOF
                     (diagnostic only: removes the lever-arm coupling, not physical)
"""
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
args = sys.argv[3:]
keep, same = None, False
i = 0
while i < len(args):
    if args[i] == "--keep":
        keep = {int(x) for x in args[i + 1].split(",")}
        i += 2
    elif args[i] == "--same-dof-only":
        same = True
        i += 1
    else:
        sys.exit(f"ABORT: unknown argument {args[i]}")

lines = src.read_text().splitlines()
out, i, found, n_kept = [], 0, False, 0


def is_kw(s):
    return s.startswith("*") and not s.startswith("**")


while i < len(lines):
    L = lines[i]
    if L.strip().upper().startswith("*EQUATION"):
        found = True
        out.append(L)
        i += 1
        k = 0
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
            k += 1
            terms = [(vals[3 * j], int(vals[3 * j + 1]), vals[3 * j + 2]) for j in range(n)]
            if keep is not None and k not in keep:
                continue
            dep = terms[0]
            new = terms
            if same:
                new = [terms[0]] + [t for t in terms[1:] if t[1] == dep[1]]
            cross = sum(1 for t in terms[1:] if t[1] != dep[1])
            print(f"equation {k}: dependent node {dep[0]} dof {dep[1]}; terms {len(terms)} -> {len(new)}; "
                  f"cross-direction terms in original {cross}; distinct wing nodes {len({t[0] for t in terms[1:]})}")
            out.append(str(len(new)))
            for j in range(0, len(new), 4):
                out.append(", ".join(f"{a}, {b}, {c}" for a, b, c in new[j:j + 4]))
            n_kept += 1
        continue
    out.append(L)
    i += 1

if not found or n_kept == 0:
    sys.exit("ABORT: no equations kept")
dst.parent.mkdir(parents=True, exist_ok=True)
dst.write_text("\n".join(out) + "\n")
print(f"kept {n_kept} equation(s); wrote {dst}")
