#!/usr/bin/env python3
"""
Stage 6 comparison figures from results/stage6_scores.csv (132 runs: 69 layouts at 1.5 mm rib,
21 selected layouts at 1.2 / 1.0 / 0.8 mm). Writes figures/stage6_comparison.png.
The applied lift is an ANALYTIC PLACEHOLDER; 20 kg mass and 290 MPa yield are assumptions;
the composite weights (mass 0.45, buckling 0.40, tip 0.15) are judgement, not derived.
Usage: /usr/bin/python3 scripts/plot_stage6_figures.py
"""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("results/stage6_scores.csv")))
for r in rows:
    for k in ("rib_t", "buckling", "tip_mm", "mass_kg", "composite"):
        r[k] = float(r[k])
    r["rank_composite"] = int(r["rank_composite"])

colors = {1.5: "#4c72b0", 1.2: "#55a868", 1.0: "#dd8452", 0.8: "#c44e52"}
fig, axs = plt.subplots(1, 3, figsize=(18, 5.6))

ax = axs[0]
for t in (1.5, 1.2, 1.0, 0.8):
    s = [r for r in rows if r["rib_t"] == t]
    ax.scatter([r["mass_kg"] for r in s], [r["buckling"] for r in s], s=26, c=colors[t],
               label="rib %g mm (%d runs)" % (t, len(s)), alpha=0.85)
for r in rows:
    if r["run"] == "L000@1.5" or r["rank_composite"] <= 3:
        ax.annotate(r["run"], (r["mass_kg"], r["buckling"]), textcoords="offset points", xytext=(5, 4), fontsize=8)
ax.set_xlabel("Wing mass [kg]")
ax.set_ylabel("Critical buckling factor [-]")
ax.set_title("Buckling vs mass, all 132 runs")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

ax = axs[1]
ids = sorted({r["id"] for r in rows if r["rib_t"] == 0.8})
hl = {"L006": "#c44e52", "L008": "#8172b2", "L000": "#000000"}
for i in ids:
    s = sorted([r for r in rows if r["id"] == i and r["rib_t"] in (1.5, 1.2, 1.0, 0.8)], key=lambda r: -r["rib_t"])
    ax.plot([r["rib_t"] for r in s], [r["buckling"] for r in s], "-o", ms=3,
            color=hl.get(i, "#bbbbbb"), lw=2.0 if i in hl else 0.9, zorder=3 if i in hl else 1,
            label=i if i in hl else None)
ax.invert_xaxis()
ax.set_xlabel("Rib thickness [mm]")
ax.set_ylabel("Critical buckling factor [-]")
ax.set_title("Buckling vs rib thickness, 21 selected layouts")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

ax = axs[2]
top = sorted(rows, key=lambda r: r["rank_composite"])[:15][::-1]
ax.barh([r["run"] for r in top], [r["composite"] for r in top],
        color=[colors[r["rib_t"]] for r in top])
ax.set_xlabel("Composite score (mass 0.45, buckling 0.40, tip 0.15)")
ax.set_title("Top 15 runs by composite score")
ax.grid(alpha=0.3, axis="x")

fig.suptitle("Project 03 Stage 6: rib placement and rib thickness (placeholder aerodynamic load)", fontsize=11)
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/stage6_comparison.png", dpi=200)
print("wrote figures/stage6_comparison.png")
