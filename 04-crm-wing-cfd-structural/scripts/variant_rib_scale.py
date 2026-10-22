#!/usr/bin/env python3
"""Diagnostic variant: scale all elastic constants of selected shell materials in a CalculiX deck.

Usage: variant_rib_scale.py <src.inp> <dst.inp> --pids 94 95 [--factor 10]
Scales E1, E2, E3, G12, G13, G23 of '*MATERIAL, NAME=M<pid>' (ENGINEERING CONSTANTS); Poisson ratios are kept.
Each material block must occur exactly once and have the expected layout, otherwise the script aborts.
"""
import argparse
import sys

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("dst")
ap.add_argument("--pids", type=int, nargs="+", required=True)
ap.add_argument("--factor", type=float, default=10.0)
a = ap.parse_args()

lines = open(a.src).read().split("\n")


def fmt(vals):
    return ", ".join(f"{v:.9e}" for v in vals)


for pid in a.pids:
    key = f"*MATERIAL, NAME=M{pid}"
    idx = [i for i, l in enumerate(lines) if l.strip().upper() == key.upper()]
    if len(idx) != 1:
        sys.exit(f"ABORT: '{key}' occurs {len(idx)} times (need exactly 1)")
    i = idx[0]
    if not lines[i + 1].upper().startswith("*ELASTIC, TYPE=ENGINEERING CONSTANTS"):
        sys.exit(f"ABORT: unexpected line after {key}: {lines[i + 1]!r}")
    l1 = [float(x) for x in lines[i + 2].split(",")]
    l2 = [float(x) for x in lines[i + 3].split(",")]
    if len(l1) != 8 or len(l2) != 2:
        sys.exit(f"ABORT: unexpected field counts for {key}: {len(l1)}, {len(l2)}")
    print(f"M{pid} before: {lines[i + 2]} | {lines[i + 3]}")
    n1 = list(l1)
    for k in (0, 1, 2, 6, 7):
        n1[k] *= a.factor
    n2 = list(l2)
    n2[0] *= a.factor
    lines[i + 2] = fmt(n1)
    lines[i + 3] = fmt(n2)
    print(f"M{pid} after : {lines[i + 2]} | {lines[i + 3]}")

open(a.dst, "w").write("\n".join(lines))
print(f"wrote {a.dst}")
