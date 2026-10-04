#!/usr/bin/env python3
"""
Stage 3 summary figures from the committed result tables (results/stage3_*.csv):
  figures/stage3_thickness_summary.png  buckling factor vs skin thickness (requirement lines),
                                        tip deflection and skin mass vs thickness
  figures/stage3_te_refinement.png      trailing-edge hot-spot refinement study at rib z=800
Run (conda env, matplotlib):  python3 scripts/plot_stage3_summary.py
"""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)


def rd(path):
    with open(path) as f:
        return list(csv.DictReader(f))


bu = rd("results/stage3_buckling_summary.csv")
st = rd("results/stage3_static_summary.csv")
te = rd("results/stage3_te_refinement.csv")

# ---------------- thickness summary ----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
for mesh, mk, lab in (("s4_k8", "o", "skin mesh 8 mm"), ("s4_k6", "s", "skin mesh 6 mm")):
    rows = [r for r in bu if r["mesh"] == mesh]
    a1.plot([float(r["skin_t_mm"]) for r in rows], [float(r["critical"]) for r in rows],
            marker=mk, ls="-" if mesh == "s4_k8" else "--", label=lab)
a1.axhline(1.5, color="r", ls="--", lw=1.2, label="requirement: factor >= 1.5")
a1.axhline(1.0, color="k", ls=":", lw=1.0, label="buckles at limit load (factor 1)")
a1.set_yscale("log")
a1.set_xlabel("skin thickness t [mm]")
a1.set_ylabel("critical buckling factor (x limit load)")
a1.set_title("Skin buckling at limit load (n = 3.8, 200 mm rib spacing)")
a1.set_xticks([0.5, 1.0, 1.5])
a1.set_yticks([1, 1.5, 2, 3, 5, 8])
a1.set_yticklabels(["1", "1.5", "2", "3", "5", "8"])
a1.minorticks_off()
a1.grid(True, which="both", alpha=0.3)
a1.legend(fontsize=8)

rows = sorted([r for r in st if r["mesh"] == "s2p5_k8" and float(r["skin_t_mm"]) in (0.5, 1.0, 1.5)],
              key=lambda r: float(r["skin_t_mm"]))
t = [float(r["skin_t_mm"]) for r in rows]
tip = [float(r["tip_uy_mm"]) for r in rows]
mass = [float(r["W_skin_N"]) / (float(r["n"]) * 9.81) for r in rows]
a2.plot(t, tip, "o-", color="C0")
a2.set_xlabel("skin thickness t [mm]")
a2.set_ylabel("tip deflection (spar axis) [mm]", color="C0")
a2.set_xticks([0.5, 1.0, 1.5])
a2.grid(True, alpha=0.3)
b2 = a2.twinx()
b2.plot(t, mass, "s--", color="C3")
b2.set_ylabel("skin mass, half-wing [kg]", color="C3")
a2.set_title("Tip deflection and skin mass")
fig.suptitle("Project 03, Stage 3: skin thickness study (aircraft mass 20 kg assumed, yield 290 MPa assumed)", fontsize=9)
fig.tight_layout()
fig.savefig("figures/stage3_thickness_summary.png", dpi=200)
plt.close(fig)

# ---------------- TE refinement ----------------
x = [float(r["local_median_edge_mm"]) for r in te]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
a1.loglog(x, [float(r["skin_peak_rib800"]) for r in te], "o-", label="skin peak (rib 800, TE vertex)")
a1.loglog(x, [float(r["solid_only_peak_rib800"]) for r in te], "s-", label="solid-only peak (rib 800)")
a1.axhline(float(te[0]["skin_peak_rib600_control"]), color="gray", ls=":", label="unrefined rib 600 (control)")
a1.invert_xaxis()
a1.set_xlabel("local median surface edge near the TE vertex [mm]  (finer ->)")
a1.set_ylabel("peak von Mises [MPa]")
a1.set_title("Peak keeps growing: singularity")
a1.grid(True, which="both", alpha=0.3)
a1.legend(fontsize=8)
for key, lab in (("ring_1_3mm", "1-3 mm from TE"), ("ring_3_6mm", "3-6 mm"),
                 ("ring_6_15mm", "6-15 mm"), ("ring_15_40mm", "15-40 mm")):
    a2.semilogx(x, [float(r[key]) for r in te], "o-", label=lab)
a2.invert_xaxis()
a2.set_xlabel("local median surface edge near the TE vertex [mm]  (finer ->)")
a2.set_ylabel("max von Mises in ring [MPa]")
a2.set_title("Away from the vertex the stress converges")
a2.grid(True, which="both", alpha=0.3)
a2.legend(fontsize=8)
fig.suptitle("Trailing-edge hot spot at rib z = 800, t = 1.0 mm, limit load (solid 4 mm, skin 8 mm)", fontsize=9)
fig.tight_layout()
fig.savefig("figures/stage3_te_refinement.png", dpi=200)
plt.close(fig)
print("wrote figures/stage3_thickness_summary.png and figures/stage3_te_refinement.png")
