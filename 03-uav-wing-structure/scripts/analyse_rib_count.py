#!/usr/bin/env python3
"""
Stage 6B analysis: buckling / tip / mass versus total rib count, using all completed 1.5 mm-rib layouts
(tier 1, replacements, rib-count study) in results/stage6_results.csv.
Writes results/stage6b_by_count.csv, results/stage6b_front.csv, figures/stage6b_rib_count.png
Run: /usr/bin/python3 scripts/analyse_rib_count.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "results/stage6_results.csv"
STAGE5 = dict(ribs=11, buckling=3.91, mass=1.723)   # uniform 100 mm spacing, from Stage 5 notes
REQ = 1.5

rows = []
for r in csv.DictReader(open(SRC)):
    if r["status"] != "complete":
        continue
    sep = ";" if ";" in r["stations"] else ","
    nrib = len([s for s in r["stations"].split(sep) if s.strip()])
    rows.append(dict(id=r["id"], family=r["family"], ribs=nrib, buck=float(r["buckling_crit"]),
                     tip=float(r["tip_mm"]), mass=float(r["mass_kg"]),
                     b2b1=float(r["b2"]) / float(r["b1"]),
                     p99=float(r["skin_p99"]) if r.get("skin_p99") else float("nan"),
                     stations=r["stations"]))
print("complete layouts used: %d ; rib counts present: %s" % (len(rows), sorted({x["ribs"] for x in rows})))

# ---- per-count summary
counts = sorted({x["ribs"] for x in rows})
by = []
for n in counts:
    g = [x for x in rows if x["ribs"] == n]
    bb = max(g, key=lambda x: x["buck"])
    bt = min(g, key=lambda x: x["tip"])
    bm = min(g, key=lambda x: x["mass"])
    by.append(dict(ribs=n, n_layouts=len(g), best_buck_id=bb["id"], best_buck=bb["buck"], best_buck_mass=bb["mass"],
                   best_buck_tip=bb["tip"], best_tip_id=bt["id"], best_tip=bt["tip"], min_mass=bm["mass"],
                   buck_per_kg=bb["buck"] / bb["mass"]))
print("\n%-5s %-4s %-6s %-8s %-8s %-8s %-6s %-8s %-8s %-8s" % (
    "ribs", "n", "bestID", "buck", "mass", "tip", "tipID", "besttip", "minmass", "buck/kg"))
for b in by:
    print("%-5d %-4d %-6s %-8.4f %-8.4f %-8.4f %-6s %-8.4f %-8.4f %-8.3f" % (
        b["ribs"], b["n_layouts"], b["best_buck_id"], b["best_buck"], b["best_buck_mass"], b["best_buck_tip"],
        b["best_tip_id"], b["best_tip"], b["min_mass"], b["buck_per_kg"]))
with open("results/stage6b_by_count.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(by[0])); w.writeheader(); w.writerows(by)

# ---- non-dominated set (buckling up, tip down, mass down), strict dominance
def dominated(a, others):
    return any(o["buck"] >= a["buck"] and o["tip"] <= a["tip"] and o["mass"] <= a["mass"] and
               (o["buck"] > a["buck"] or o["tip"] < a["tip"] or o["mass"] < a["mass"]) for o in others)
front = [x for x in rows if not dominated(x, rows)]
front.sort(key=lambda x: x["mass"])
print("\nnon-dominated layouts (strict, no noise tolerance): %d of %d" % (len(front), len(rows)))
print("%-5s %-5s %-8s %-8s %-8s %-6s %s" % ("id", "ribs", "buck", "tip", "mass", "b2/b1", "stations"))
for x in front:
    print("%-5s %-5d %-8.4f %-8.4f %-8.4f %-6.3f %s" % (x["id"], x["ribs"], x["buck"], x["tip"], x["mass"], x["b2b1"], x["stations"]))
with open("results/stage6b_front.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(front[0])); w.writeheader(); w.writerows(front)

print("\nminimum buckling over all layouts: %.4f (requirement %.1f) ; max skin p99 %.2f MPa ; min %.2f MPa" % (
    min(x["buck"] for x in rows), REQ, max(x["p99"] for x in rows if x["p99"] == x["p99"]),
    min(x["p99"] for x in rows if x["p99"] == x["p99"])))

# ---- figure
os.makedirs("figures", exist_ok=True)
fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
cm = plt.get_cmap("viridis")
def scat(a, key, ylabel):
    a.scatter([x["ribs"] for x in rows], [x[key] for x in rows], s=18, alpha=0.55, color="#4c72b0", edgecolor="none")
    a.set_xlabel("total ribs"); a.set_ylabel(ylabel); a.grid(alpha=0.3); a.set_xticks(counts)
scat(ax[0, 0], "buck", "critical buckling factor")
ax[0, 0].plot([b["ribs"] for b in by], [b["best_buck"] for b in by], "o-", color="#c44e52", label="best per rib count")
ax[0, 0].plot([STAGE5["ribs"]], [STAGE5["buckling"]], "s", color="k", label="Stage 5 uniform (11)")
ax[0, 0].legend(fontsize=8)
scat(ax[0, 1], "tip", "tip deflection [mm]")
ax[0, 1].plot([b["ribs"] for b in by], [b["best_tip"] for b in by], "o-", color="#c44e52", label="best per rib count")
ax[0, 1].legend(fontsize=8)
scat(ax[1, 0], "mass", "structural mass [kg]")
ax[1, 0].plot([b["ribs"] for b in by], [b["min_mass"] for b in by], "o-", color="#c44e52", label="lightest per rib count")
ax[1, 0].legend(fontsize=8)
sc = ax[1, 1].scatter([x["mass"] for x in rows], [x["buck"] for x in rows], c=[x["ribs"] for x in rows], cmap=cm, s=20)
ax[1, 1].scatter([x["mass"] for x in front], [x["buck"] for x in front], s=70, facecolor="none", edgecolor="k", label="non-dominated (3 obj.)")
ax[1, 1].plot([STAGE5["mass"]], [STAGE5["buckling"]], "s", color="r", label="Stage 5 uniform (11)")
ax[1, 1].set_xlabel("structural mass [kg]"); ax[1, 1].set_ylabel("critical buckling factor"); ax[1, 1].grid(alpha=0.3)
ax[1, 1].legend(fontsize=8); fig.colorbar(sc, ax=ax[1, 1], label="total ribs")
fig.suptitle("Stage 6: rib count and placement, 1.5 mm ribs, 1.0 mm skin, n = 3.8 limit load (Stage 4 placeholder pressure)")
fig.savefig("figures/stage6b_rib_count.png", dpi=170)
print("\nwrote figures/stage6b_rib_count.png, results/stage6b_by_count.csv, results/stage6b_front.csv")
