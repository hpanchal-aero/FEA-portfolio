"""
Project 01 - Aerospace Mounting Bracket
Primary analytical verification: sigma_xx (longitudinal bending
stress) vs Euler-Bernoulli theory, at the same x-stations used in
the biaxial-stress investigation (check_biaxial_stress.py).

This is the correct headline verification comparison, since
Euler-Bernoulli theory directly predicts sigma_xx = M(x)c/I, not
von Mises stress or a windowed maximum (see chat log discussion).

Also generates the sigma_yy/sigma_xx vs x plot showing the
Poisson-ratio-governed biaxial stress development toward the root.

Uses ONLY the already-converged 0.3mm result - no re-solving.
"""

import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, ".")
import hex_lib as lib

CASE_DIR = "../simulation/convergence_hex/mesh_0.3mm"
FIG_DIR = "../figures"
RESULTS_DIR = "../results"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

L, B, H = lib.L, lib.B, lib.H
F = 235.44
E = 71700.0
NU = 0.33
I_SECTION = B * H**3 / 12.0
C = H / 2.0

def sigma_theory(x):
    M = F * (L - x)
    return M * C / I_SECTION

# Same station set as the biaxial-stress investigation
X_STATIONS = [3, 6, 9, 12, 20, 30, 40]

clean_inp = os.path.join(CASE_DIR, "mesh_clean.inp")
frd_path = os.path.join(CASE_DIR, "analysis.frd")

nodes, elements = lib.parse_clean_mesh(clean_inp)
with open(frd_path) as f:
    frd_lines = f.readlines()
stress = lib.parse_frd_block(frd_lines, "STRESS", 6)

TOL = 0.05

def find_node(xs):
    candidates = [nid for nid, (x, y, z) in nodes.items()
                  if abs(x - xs) < TOL and abs(y - B/2) < TOL and abs(z - H) < TOL]
    if not candidates:
        candidates = [nid for nid, (x, y, z) in nodes.items()
                      if abs(x - xs) < 0.4 and abs(y - B/2) < 0.4 and abs(z - H) < 0.4]
    return candidates[0] if candidates else None

rows = []
for xs in X_STATIONS:
    nid = find_node(xs)
    if nid is None or nid not in stress:
        print(f"WARNING: no valid node/stress data found at x={xs}mm - "
              f"reporting limitation, not inventing a value.")
        continue
    sxx, syy, szz, sxy, syz, szx = stress[nid]
    th = sigma_theory(xs)
    abs_diff = sxx - th
    pct_diff = abs_diff / th * 100
    ratio = syy / sxx if abs(sxx) > 1e-6 else float('nan')
    rows.append({
        "x_mm": xs,
        "sigma_xx_theory_mpa": round(th, 3),
        "sigma_xx_fea_mpa": round(sxx, 3),
        "abs_diff_mpa": round(abs_diff, 3),
        "pct_diff": round(pct_diff, 2),
        "sigma_yy_fea_mpa": round(syy, 3),
        "ratio_yy_xx": round(ratio, 4),
    })

# --- Write CSV table ---
csv_path = os.path.join(RESULTS_DIR, "sigma_xx_verification.csv")
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote: {csv_path}")

# --- Print table ---
print(f"\n{'x (mm)':>8} {'theory (MPa)':>14} {'FEA (MPa)':>12} "
      f"{'abs diff':>10} {'% diff':>8}")
for r in rows:
    print(f"{r['x_mm']:>8} {r['sigma_xx_theory_mpa']:>14.3f} "
          f"{r['sigma_xx_fea_mpa']:>12.3f} {r['abs_diff_mpa']:>10.3f} "
          f"{r['pct_diff']:>7.2f}%")

# --- Plot 1: sigma_xx FEA vs theory ---
x_vals = [r["x_mm"] for r in rows]
theory_vals = [r["sigma_xx_theory_mpa"] for r in rows]
fea_vals = [r["sigma_xx_fea_mpa"] for r in rows]

# Dense theory curve for a smooth reference line
x_dense = np.linspace(0.5, 55, 200)
theory_dense = [sigma_theory(x) for x in x_dense]

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(x_dense, theory_dense, "-", color="gray", linewidth=1.5,
        label="Euler-Bernoulli theory ($\\sigma_{xx} = M(x)c/I$)")
ax.plot(x_vals, fea_vals, "o", color="#1f77b4", markersize=8,
        label="FEA $\\sigma_{xx}$ (top fiber, mid-width, 0.3mm mesh)")
for x, th, fe in zip(x_vals, theory_vals, fea_vals):
    ax.annotate(f"{fe:.1f}", (x, fe), textcoords="offset points",
                xytext=(6, -12), fontsize=8, color="#1f77b4")
ax.set_xlabel("Distance from fixed root, x (mm)")
ax.set_ylabel("$\\sigma_{xx}$ (MPa)")
ax.set_title("Longitudinal Bending Stress: FEA vs Euler-Bernoulli Theory")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "sigma_xx_verification.png"), dpi=150)
plt.close(fig)
print("\nWrote sigma_xx_verification.png")

# --- Plot 2: sigma_yy/sigma_xx ratio vs x, showing Poisson decay ---
ratio_vals = [r["ratio_yy_xx"] for r in rows]

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(x_vals, ratio_vals, "o-", color="#d62728", markersize=8,
        label="FEA $\\sigma_{yy}/\\sigma_{xx}$ ratio")
ax.axhline(NU, color="gray", linestyle="--",
           label=f"Poisson's ratio target ($\\nu$ = {NU})")
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xlabel("Distance from fixed root, x (mm)")
ax.set_ylabel("$\\sigma_{yy} / \\sigma_{xx}$")
ax.set_title("Biaxial Stress Ratio: Anticlastic Curvature Restraint Decay")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "biaxial_ratio_decay.png"), dpi=150)
plt.close(fig)
print("Wrote biaxial_ratio_decay.png")

print(f"\nAll outputs saved to {os.path.abspath(FIG_DIR)} and "
      f"{os.path.abspath(RESULTS_DIR)}")
