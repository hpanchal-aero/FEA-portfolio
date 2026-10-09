#!/usr/bin/env python3
"""
Stage 5 comparison figure from results/stage5_comparison.csv (200 mm vs 100 mm rib spacing).
Run with the conda env python (numpy + matplotlib):  python3 scripts/plot_stage5_figures.py
Writes figures/stage5_comparison.png. The applied lift is an ANALYTIC PLACEHOLDER; 20 kg mass and 290 MPa yield
are assumptions; 50 mm and 20 mm were not run (11 GB RAM).
"""
import csv
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    path = "results/stage5_comparison.csv"
    if not os.path.isfile(path):
        print(f"ERROR: missing {path}")
        sys.exit(1)
    rows = {r["label"]: r for r in csv.DictReader(open(path))}
    a, b = rows["200 mm"], rows["100 mm"]
    g = lambda r, k: float(r[k])
    la = f"200 mm (6 ribs, {g(a, 'mass_total_kg'):.3f} kg)"
    lb = f"100 mm (11 ribs, {g(b, 'mass_total_kg'):.3f} kg)"
    ca, cb = "#4C72B0", "#DD8452"
    w = 0.38
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.0))

    ax = axs[0]
    k = np.arange(1, 7)
    ya = [g(a, f"buckling_mode{i}") for i in k]
    yb = [g(b, f"buckling_mode{i}") for i in k]
    ax.bar(k - w / 2, ya, w, color=ca, label=la)
    ax.bar(k + w / 2, yb, w, color=cb, label=lb)
    for x, y in zip(k - w / 2, ya):
        ax.text(x, y + 0.08, f"{y:.2f}", ha="center", fontsize=7)
    for x, y in zip(k + w / 2, yb):
        ax.text(x, y + 0.08, f"{y:.2f}", ha="center", fontsize=7)
    ax.axhline(1.5, color="k", ls="--", lw=1, label="requirement 1.5")
    ax.set_xlabel("buckling mode number")
    ax.set_ylabel("buckling factor (multiplies the limit load)")
    ax.set_title("Linear buckling factors", fontsize=10)
    ax.set_ylim(0, 8.2)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3, axis="y")

    ax = axs[1]
    cats = ["lift only", "gravity only", "full"]
    keys = [("tip_lift_only_mm", "beam_lift_with_rib_inertia_mm"),
            ("tip_gravity_only_mm", "beam_grav_with_rib_inertia_mm"),
            ("tip_full_mm", "beam_full_with_rib_inertia_mm")]
    x = np.arange(3)
    fa = [g(a, kk[0]) for kk in keys]
    fb = [g(b, kk[0]) for kk in keys]
    ba = [g(a, kk[1]) for kk in keys]
    bb = [g(b, kk[1]) for kk in keys]
    ax.bar(x - w / 2, fa, w, color=ca, label="FEA, 200 mm")
    ax.bar(x + w / 2, fb, w, color=cb, label="FEA, 100 mm")
    ax.plot(x - w / 2, ba, "kD", ms=5, label="beam reference with rib-section inertia")
    ax.plot(x + w / 2, bb, "kD", ms=5)
    for xx, y in list(zip(x - w / 2, fa)) + list(zip(x + w / 2, fb)):
        ax.text(xx, y + (0.3 if y >= 0 else -0.22), f"{y:.3f}", ha="center",
                va="bottom" if y >= 0 else "top", fontsize=7)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("tip deflection of the spar axis [mm]")
    ax.set_title("Tip deflection at limit load (FEA vs beam)", fontsize=10)
    ax.set_ylim(-1.6, 6.4)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, axis="y")

    ax = axs[2]
    items = [("total mass", "mass_total_kg", False),
             ("tip deflection, full", "tip_full_mm", False),
             ("tip deflection, lift only", "tip_lift_only_mm", False),
             ("critical buckling factor", "buckling_critical", False),
             ("tip twist", "tip_twist_deg", False),
             ("skin p99 von Mises", "skin_p99_MPa", False),
             ("design stress, dTE >= 6 mm variant", "design_max_excl_z_ge_990_dTE_lt_6_posthoc_MPa", False),
             ("design stress, dTE >= 3 mm\n(mesh-sensitive, not a stress increase)", "design_stress_dTE_ge3_MPa", True)]
    vals = [100.0 * (g(b, key) / g(a, key) - 1.0) for _, key, _ in items]
    ypos = np.arange(len(items))[::-1]
    for yy, v, (lab, key, hatch) in zip(ypos, vals, items):
        ax.barh(yy, v, color="0.85" if hatch else "0.45", hatch="//" if hatch else None,
                edgecolor="k" if hatch else None)
        ax.text(v + (0.8 if v >= 0 else -0.8), yy, f"{v:+.1f}%", va="center", ha="left" if v >= 0 else "right",
                fontsize=8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([lab for lab, _, _ in items], fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlim(-12, 32)
    ax.set_xlabel("change from 200 mm to 100 mm spacing [%]")
    ax.set_title("Change from 200 mm to 100 mm spacing", fontsize=10)
    ax.grid(alpha=0.3, axis="x")

    fig.suptitle("Stage 5: rib spacing 200 mm vs 100 mm (skin t = 1.0 mm, Stage 4 ANALYTIC PLACEHOLDER pressure; "
                 "20 kg and 290 MPa assumed; 50 and 20 mm not run)", fontsize=10)
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/stage5_comparison.png", dpi=200)
    print("wrote figures/stage5_comparison.png")


if __name__ == "__main__":
    main()
