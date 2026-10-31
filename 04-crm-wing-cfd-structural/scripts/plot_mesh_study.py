#!/usr/bin/python3
"""Mesh-study figure for Project 04 (first-order upwind CFD, alpha = 0.36 deg, levels L8 / L9 / L10).

usage: plot_mesh_study.py [results/mesh_study_p0p36.txt] [figures/mesh_study_p0p36.png]

Reads the table and the volume-based GCI table written by mesh_study.py (nothing is retyped):
Cl, Cd, Cm at the common iteration and the CalculiX tip uz under each mesh's loads.
Panel: the three levels vs cell count, the Richardson-extrapolated value (observed order, dashed line) and the
GCI21 (safety factor 1.25, observed order) as an error bar on the fine mesh.  The GCI with the formal order p = 1
is printed in each panel as well.
"""
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

src = sys.argv[1] if len(sys.argv) > 1 else "results/mesh_study_p0p36.txt"
out = sys.argv[2] if len(sys.argv) > 2 else "figures/mesh_study_p0p36.png"

lev, gci, sect = {}, {}, None
for line in open(src):
    m = re.match(r"^(L\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*\|", line)
    if m:
        parts = line.split("|")
        a = [float(t) for t in parts[1].split()]
        tip = [float(t) for t in parts[4].split()]
        lev[m.group(1)] = dict(cells=int(m.group(2)), Cl=a[0], Cd=a[1], Cm=a[2], tip=tip[0])
        continue
    if line.startswith("=== GCI, volume"):
        sect = "volume"
        continue
    if line.startswith("=== GCI, surface"):
        sect = "surface"
        continue
    if sect == "volume" and "|" in line:
        parts = line.split("|")
        left = parts[0].split()
        try:
            vals = [float(t) for t in left[-3:]]
            mid = [float(t) for t in parts[1].split()]
            p1 = float(parts[2].split()[0])
        except (ValueError, IndexError):
            continue
        gci[" ".join(left[:-3])] = dict(fine=vals[0], p=mid[0], extrap=mid[1], gci=mid[3], gci_p1=p1)

names = ["L8", "L9", "L10"]
need = [("Cl @1148", "Cl", "Cl", "Cl"),
        ("Cd @1148", "Cd", "Cd", "Cd (shown for completeness, not used: Euler drag is numerical)"),
        ("Cm @1148", "Cm", "Cm", "Cm"),
        ("tip uz", "tip", "tip uz (m)", "Tip uz under each mesh's loads (CalculiX, uCRM-9 coarse deck)")]
for n in names:
    if n not in lev:
        sys.exit(f"ABORT: level {n} not found in {src}")
for gk, *_ in need:
    if gk not in gci:
        sys.exit(f"ABORT: '{gk}' not found in the volume GCI table of {src}")

x = np.array([lev[n]["cells"] for n in names]) / 1e6
fig, axs = plt.subplots(2, 2, figsize=(11, 7.6))
for ax, (gk, vk, yl, title) in zip(axs.ravel(), need):
    y = np.array([lev[n][vk] for n in names])
    g = gci[gk]
    err = g["gci"] / 100.0 * abs(y[-1])
    st = [100.0 * (abs(y[i + 1]) - abs(y[i])) / abs(y[i]) for i in range(2)]
    print(f"{gk:9s} values {y[0]:9.5f} {y[1]:9.5f} {y[2]:9.5f}   |value| steps {st[0]:+.2f}% {st[1]:+.2f}%   "
          f"p_obs {g['p']:.2f}  extrapolated {g['extrap']:.5f}  GCI21 {g['gci']:.2f}% (p_obs) / {g['gci_p1']:.2f}% (p=1)")
    ax.plot(x, y, "o-", color="#1f4e79", lw=1.6, ms=7, zorder=3)
    ax.errorbar(x[-1], y[-1], yerr=err, color="#1f4e79", capsize=5, lw=1.4, zorder=3)
    ax.axhline(g["extrap"], color="0.45", ls="--", lw=1.2, zorder=1)
    ax.text(x[0], g["extrap"], " Richardson extrapolation", fontsize=8, color="0.35",
            va="bottom", ha="left")
    for n, xi, yi in zip(names, x, y):
        ax.annotate(n, (xi, yi), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9)
    lo = min(y.min(), g["extrap"], y[-1] - err)
    hi = max(y.max(), g["extrap"], y[-1] + err)
    pad = 0.35 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(x[0] - 0.12, x[-1] + 1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{n}\n{c:.2f} M" for n, c in zip(names, x)])
    up = y[-1] > y[0]
    ax.text(0.97, 0.05 if up else 0.95,
            f"observed order p = {g['p']:.2f}\nextrapolated {g['extrap']:.4f}\n"
            f"GCI21 {g['gci']:.1f}% (p_obs)\nGCI21 {g['gci_p1']:.1f}% (p = 1)",
            transform=ax.transAxes, ha="right", va="bottom" if up else "top", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.set_ylabel(yl)
    ax.set_xlabel("cells")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
fig.suptitle("Mesh study: first-order upwind CFD, alpha = 0.36 deg (common iteration 1148)", fontsize=12)
fig.text(0.5, 0.008, "Three non-nested meshes, refinement ratio 1.20 / 1.22 (volume). Observed orders above 1 for a first-order scheme "
         "mean the asymptotic range is not reached:\nthe p = 1 GCI is the conservative one. Error bar = GCI21 with the observed order.",
         ha="center", fontsize=8.5)
fig.tight_layout(rect=(0, 0.06, 1, 0.96))
fig.savefig(out, dpi=200)
print("saved", out)
