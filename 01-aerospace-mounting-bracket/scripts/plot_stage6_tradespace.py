"""
Stage 6 -- mass-vs-constraint trade-space plot.

Combines the 14-point fixed-scheme (fine L3) dataset and the 16-point
high-t dataset (2 spot-checks + 9 sub-region, all coarse L1 mesh) on
one figure, with resolution level distinguished by marker shape (NOT
blended into one claimed-equal-resolution dataset) and feasibility by
color. Constraint limits drawn as reference lines. lhs_06 highlighted
as the selected design.
"""

import csv
import matplotlib.pyplot as plt

STRESS_LIMIT = 335.0
DISP_LIMIT = 1.5

FIXED_SCHEME_IDS = ["baseline", "lhs_01", "lhs_02", "lhs_03", "lhs_04", "lhs_05",
                     "lhs_07", "lhs_08", "lhs_09", "lhs_12", "lhs_13", "lhs_14",
                     "lhs_15", "lhs_16", "lhs_17", "lhs_18", "lhs_19", "lhs_20",
                     "lhs_21", "lhs_22", "lhs_23", "lhs_24", "lhs_25"]
COARSE_IDS_SUFFIX = "_COARSE"
SUB_PREFIX = "sub_"

BASELINE_MASS = 26859.5774

geom = {}
with open("results/stage6_parametric/geometry_validation_results.csv", newline="") as f:
    for r in csv.DictReader(f):
        geom[r["design_id"]] = r

solved = []
with open("results/stage6_parametric/solve_results.csv", newline="") as f:
    for r in csv.DictReader(f):
        if r["solver_status"] == "ok" and r["parser_status"] == "ok":
            solved.append(r)

fixed_pts, coarse_pts = [], []
for r in solved:
    did = r["design_id"]
    ux = abs(float(r["max_ux_mm"]))
    vm = float(r["peak_vm_mpa"])

    if did == "baseline":
        mass = BASELINE_MASS
    elif did.endswith(COARSE_IDS_SUFFIX):
        gid = did.replace(COARSE_IDS_SUFFIX, "")
        mass = float(geom[gid]["mass"])
    elif did.startswith(SUB_PREFIX) or did == "lhs_06_L2":
        gid = "lhs_06" if did == "lhs_06_L2" else did
        mass = float(geom[gid]["mass"])
    else:
        mass = float(geom[did]["mass"]) if did in geom else None

    if mass is None:
        continue

    feas = (ux <= DISP_LIMIT) and (vm <= STRESS_LIMIT)
    point = (did, mass, ux, vm, feas)

    if did in FIXED_SCHEME_IDS:
        fixed_pts.append(point)
    elif did.endswith(COARSE_IDS_SUFFIX) or did.startswith(SUB_PREFIX):
        coarse_pts.append(point)
    # lhs_06_L2 excluded from both -- plotted separately as the selected design

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, ykey, ylabel, limit in [
    (axes[0], 2, "Peak von Mises (MPa)", STRESS_LIMIT),
    (axes[1], 3, "|Tip Displacement| (mm)", DISP_LIMIT),
]:
    ykey_idx = ykey
    for did, mass, ux, vm, feas in fixed_pts:
        y = vm if ykey_idx == 2 else ux
        ax.scatter(mass, y, marker="o", s=60,
                   color="tab:green" if feas else "tab:red",
                   edgecolors="k", linewidths=0.5, zorder=3,
                   label="_nolegend_")
    for did, mass, ux, vm, feas in coarse_pts:
        y = vm if ykey_idx == 2 else ux
        ax.scatter(mass, y, marker="^", s=70,
                   color="tab:green" if feas else "tab:red",
                   edgecolors="k", linewidths=0.5, zorder=3,
                   label="_nolegend_")

    # selected design
    ux6, vm6, mass6 = 1.39601, 188.4064, 31137.00
    y6 = vm6 if ykey_idx == 2 else ux6
    ax.scatter(mass6, y6, marker="*", s=400, color="gold", edgecolors="k",
               linewidths=1.2, zorder=5, label="Selected: lhs_06 (L2)")

    ax.axhline(limit, color="black", linestyle="--", linewidth=1, label=f"Limit ({limit})")
    ax.set_xlabel("Mass (mm³)")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)

axes[0].set_title("Stress vs. Mass")
axes[1].set_title("Displacement vs. Mass")

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="tab:green",
           markeredgecolor="k", markersize=9, label="Fixed-scheme (fine), feasible"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="tab:red",
           markeredgecolor="k", markersize=9, label="Fixed-scheme (fine), infeasible"),
    Line2D([0], [0], marker="^", color="w", markerfacecolor="tab:green",
           markeredgecolor="k", markersize=9, label="High-t region (coarse), feasible"),
    Line2D([0], [0], marker="^", color="w", markerfacecolor="tab:red",
           markeredgecolor="k", markersize=9, label="High-t region (coarse), infeasible"),
    Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
           markeredgecolor="k", markersize=16, label="Selected: lhs_06 (L2, converged)"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.08))

plt.suptitle("Stage 6 Trade Space -- circles: fine (L3) mesh; triangles: coarse (L1) mesh\n"
             "(two distinct mesh resolutions, not directly comparable point-for-point)")
plt.tight_layout()
plt.savefig("figures/stage6_tradespace.png", dpi=150, bbox_inches="tight")
print("Wrote: figures/stage6_tradespace.png")
print(f"Fixed-scheme points plotted: {len(fixed_pts)}")
print(f"Coarse high-t points plotted: {len(coarse_pts)}")
