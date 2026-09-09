"""
Stage 4b - toe convergence figures. Reuses the parsing logic from
diagnose_toe_convergence.py, run across all three solved levels, to
produce two figures:
  1. Peak toe stress vs. mesh refinement (non-convergence evidence)
  2. Stress vs. distance from the toe, overlaid across levels
     (Saint-Venant-style decay, converged away from the singularity)
"""

import numpy as np
import matplotlib.pyplot as plt

LEVELS = [
    ("Baseline (~0.5mm)", ".", 0.5, 8000),
    ("toe_fine (0.15mm)", "toe_convergence/level_toe_fine", 0.15, 83082),
    ("toe_intermediate (0.10mm)", "toe_convergence/level_toe_intermediate", 0.10, 139252),
]


def parse_frd(path):
    coords, stress = {}, {}
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("    2C"):
            i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                coords[nid] = (float(rec[13:25]), float(rec[25:37]), float(rec[37:49]))
                i += 1
            continue
        if "STRESS" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                vals = [float(rec[13 + 12*k:25 + 12*k]) for k in range(6)]
                stress[nid] = tuple(vals)
                i += 1
            continue
        i += 1
    return coords, stress


peak_values = []
element_counts = []
trend_data = []

for label, level_dir, toe_size, n_elem in LEVELS:
    frd_path = f"{level_dir}/analysis_frame_gusset.frd"
    print(f"Reading {frd_path}")
    coords, stress = parse_frd(frd_path)

    toe_candidates = []
    for nid, (x, y, z) in coords.items():
        if 37.5 <= x <= 38.5 and 1.5 <= z <= 2.5:
            if nid in stress:
                toe_candidates.append(abs(stress[nid][0]))
    peak = max(toe_candidates)
    peak_values.append(peak)
    element_counts.append(n_elem)
    print(f"  {label}: peak |sigma_xx| = {peak:.2f} MPa")

    x_stations = [38, 36, 34, 32, 30, 25, 20]
    trend = []
    for x_target in x_stations:
        matches = []
        for nid, (x, y, z) in coords.items():
            if abs(x - x_target) < 0.5 and abs(z - 2.0) < 0.3 and abs(y - 20.0) < 1.5:
                if nid in stress:
                    matches.append(stress[nid][0])
        trend.append(np.mean(matches) if matches else np.nan)
    trend_data.append((label, x_stations, trend))

# ---- Figure 1: peak toe stress vs. mesh refinement ----
fig1, ax1 = plt.subplots(figsize=(7, 5))
ax1.plot(element_counts, peak_values, 'o-', color='crimson', markersize=8, linewidth=2)
for ec, pv, (label, *_) in zip(element_counts, peak_values, LEVELS):
    ax1.annotate(label.split(' (')[0], (ec, pv), textcoords="offset points",
                 xytext=(8, -4), fontsize=9)
ax1.set_xlabel("Total elements in mesh")
ax1.set_ylabel("Peak |sigma_xx| at gusset toe (MPa)")
ax1.set_title("Stage 4b: Gusset Toe Peak Stress vs. Mesh Refinement\n(monotonic growth = non-convergent singularity)")
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig("01-aerospace-mounting-bracket/figures/stage4b_toe_nonconvergence.png", dpi=150)
print("\nSaved: stage4b_toe_nonconvergence.png")

# ---- Figure 2: stress vs. distance from toe, overlaid ----
fig2, ax2 = plt.subplots(figsize=(8, 5.5))
colors = ['tab:blue', 'tab:orange', 'tab:green']
for (label, x_stations, trend), color in zip(trend_data, colors):
    distance_from_toe = [38 - x for x in x_stations]
    ax2.plot(distance_from_toe, trend, 'o-', color=color, label=label, linewidth=2, markersize=6)

ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5, label='Toe location (x=38)')
ax2.set_xlabel("Distance from toe, along Flange A (mm)")
ax2.set_ylabel("sigma_xx (MPa)")
ax2.set_title("Stage 4b: Stress vs. Distance from Gusset Toe\n(converged everywhere except at the toe itself)")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig("01-aerospace-mounting-bracket/figures/stage4b_toe_decay_profile.png", dpi=150)
print("Saved: stage4b_toe_decay_profile.png")
