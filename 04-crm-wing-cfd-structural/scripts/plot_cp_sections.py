#!/usr/bin/env python3
"""Chordwise Cp on the wing surface at several span stations (Euler runs).

Usage: plot_cp_sections.py [time] [case_dir] [tag]
Reads VTK/wing/wing_<time>.vtk and writes figures/cfd_cp_sections_<tag>.png.
Cp = (p - pInf)/qInf. Prints a table of chord and extreme Cp per station.
"""
import sys
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

time = sys.argv[1] if len(sys.argv) > 1 else "600"
case = sys.argv[2] if len(sys.argv) > 2 else "simulation/cfd_v1/euler"
tag = sys.argv[3] if len(sys.argv) > 3 else f"t{time}"

P_INF, Q_INF = 21662.7, 10955.9
CP_FLOOR = (2166.27 - P_INF) / Q_INF          # limitPressure floor as Cp
stations = [2.0, 6.0, 10.0, 14.0, 18.0, 22.0, 26.3, 29.0]

w = pv.read(f"{case}/VTK/wing/wing_{time}.vtk").cell_data_to_point_data()

fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=True)
print(f"{'y':>6} {'chord':>7} {'x_LE':>7} {'Cp_min':>8} {'x/c@min':>8} "
      f"{'n@floor':>8}")
for ax, y0 in zip(axes.ravel(), stations):
    s = w.slice(normal=(0, 1, 0), origin=(0, y0, 0))
    if s.n_points == 0:
        ax.set_title(f"y = {y0} m (no data)")
        continue
    pts = np.asarray(s.points)
    cp = (np.asarray(s.point_data["p"]) - P_INF) / Q_INF
    x, z = pts[:, 0], pts[:, 2]
    iLE, iTE = int(np.argmin(x)), int(np.argmax(x))
    chord = x[iTE] - x[iLE]
    xc = (x - x[iLE]) / chord
    zline = z[iLE] + (z[iTE] - z[iLE]) * xc
    up = z >= zline
    ax.plot(xc[up], cp[up], ".", ms=2, label="upper")
    ax.plot(xc[~up], cp[~up], ".", ms=2, label="lower")
    ax.axhline(CP_FLOOR, color="k", ls="--", lw=0.8)
    ax.set_title(f"y = {y0} m, c = {chord:.2f} m")
    ax.set_xlabel("x/c")
    ax.grid(alpha=0.3)
    k = int(np.argmin(cp))
    print(f"{y0:6.1f} {chord:7.3f} {x[iLE]:7.2f} {cp[k]:8.3f} {xc[k]:8.3f} "
          f"{int(np.sum(cp <= CP_FLOOR + 0.002)):8d}")
axes[0, 0].invert_yaxis()
axes[0, 0].set_ylabel("Cp")
axes[1, 0].set_ylabel("Cp")
axes[0, 0].legend()
fig.suptitle(f"Cp sections, time {time} ({tag}); dashed = pressure floor")
fig.tight_layout()
out = f"figures/cfd_cp_sections_{tag}.png"
fig.savefig(out, dpi=150)
print("saved", out)
