"""
Project 01 - Aerospace Mounting Bracket
Generate convergence plots from results/convergence_study_hex.csv
using Matplotlib - the authoritative, reproducible plotting pipeline
for this project (per project policy: only script-generated plots,
saved to figures/, are treated as trusted results).
"""

import csv
import matplotlib.pyplot as plt
import os

CSV_PATH = "../results/convergence_study_hex.csv"
FIG_DIR = "../figures"
os.makedirs(FIG_DIR, exist_ok=True)

labels, n_elements, tip_uz, root_vm = [], [], [], []
with open(CSV_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        labels.append(row["label"])
        n_elements.append(int(row["n_elements"]))
        tip_uz.append(float(row["tip_uz_avg_mm"]))
        root_vm.append(float(row["max_vm_root_mpa"]))

TARGET_UZ = -1.108
TARGET_VM = 132.4

# --- Displacement convergence ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(n_elements, tip_uz, "o-", color="#1f77b4", label="FEA tip Uz (avg)")
ax.axhline(TARGET_UZ, color="gray", linestyle="--",
           label=f"Euler-Bernoulli target ({TARGET_UZ} mm)")
ax.set_xscale("log")
ax.set_xlabel("Number of elements (log scale)")
ax.set_ylabel("Tip Uz (mm)")
ax.set_title("Tip Displacement Convergence - Structured Hex Mesh")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
for x, y, lbl in zip(n_elements, tip_uz, labels):
    ax.annotate(lbl, (x, y), textcoords="offset points", xytext=(0, 8),
                fontsize=8, ha="center")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "convergence_displacement.png"), dpi=150)
plt.close(fig)

# --- Stress convergence ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(n_elements, root_vm, "o-", color="#d62728",
        label="FEA root-band max von Mises")
ax.axhline(TARGET_VM, color="gray", linestyle="--",
           label=f"Euler-Bernoulli target ({TARGET_VM} MPa)")
ax.set_xscale("log")
ax.set_xlabel("Number of elements (log scale)")
ax.set_ylabel("Max von Mises in root band [3,9]mm (MPa)")
ax.set_title("Root-Band Stress Convergence - Structured Hex Mesh")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
for x, y, lbl in zip(n_elements, root_vm, labels):
    ax.annotate(lbl, (x, y), textcoords="offset points", xytext=(0, 8),
                fontsize=8, ha="center")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "convergence_stress.png"), dpi=150)
plt.close(fig)

print("Wrote convergence_displacement.png")
print("Wrote convergence_stress.png")
print(f"Saved to {os.path.abspath(FIG_DIR)}")
