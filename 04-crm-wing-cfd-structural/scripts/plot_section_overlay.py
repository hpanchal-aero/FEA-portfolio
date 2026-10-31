#!/usr/bin/env python3
"""Wing sections with clipped cells overlaid (Euler runs).

Usage: plot_section_overlay.py <time> <case_dir> <tag> [y1 y2 ...]

Draws x-z outlines of the wing surface at the given span stations and marks
cell centres within +-0.15 m of each station that sit at the pressure floor
(blue), the temperature floor (green) or the temperature ceiling (red).
Writes figures/cfd_section_overlay_<tag>.png and prints the counts.
"""
import glob
import sys
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

time = sys.argv[1]
case = sys.argv[2]
tag = sys.argv[3]
stations = [float(v) for v in sys.argv[4:]] or [0.5, 2.0, 4.0, 6.0, 26.3, 29.0]

P_FLOOR, T_FLOOR, T_CEIL = 2166.27, 100.0, 400.0

mesh = pv.read(sorted(glob.glob(f"{case}/VTK/*_{time}.vtk"))[0])
p = np.asarray(mesh.cell_data["p"])
T = np.asarray(mesh.cell_data["T"])
ctr = mesh.cell_centers().points
wing = pv.read(f"{case}/VTK/wing/wing_{time}.vtk")

fl = p <= P_FLOOR * 1.0005
cold = T <= T_FLOOR + 0.05
hot = T >= T_CEIL - 0.05

ncol = 3
nrow = int(np.ceil(len(stations) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(16, 5 * nrow))
print(f"{'y':>6} {'p floor':>8} {'T floor':>8} {'T ceil':>8}")
for ax, y0 in zip(np.atleast_1d(axes).ravel(), stations):
    s = wing.slice(normal=(0, 1, 0), origin=(0, y0, 0))
    pts = np.asarray(s.points)
    band = np.abs(ctr[:, 1] - y0) < 0.15
    ax.plot(pts[:, 0], pts[:, 2], "k.", ms=1.5)
    for mask, col, lab in ((band & fl, "tab:blue", "p floor"),
                           (band & cold, "tab:green", "T floor"),
                           (band & hot, "tab:red", "T ceiling")):
        ax.plot(ctr[mask, 0], ctr[mask, 2], ".", color=col, ms=3, label=lab)
    ax.set_aspect("equal")
    if pts.size:
        ax.set_xlim(pts[:, 0].min() - 1.0, pts[:, 0].max() + 1.0)
        zc = 0.5 * (pts[:, 2].min() + pts[:, 2].max())
        ax.set_ylim(zc - 2.5, zc + 2.5)
    ax.set_title(f"y = {y0} m")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.grid(alpha=0.3)
    print(f"{y0:6.1f} {int((band & fl).sum()):8d} {int((band & cold).sum()):8d} "
          f"{int((band & hot).sum()):8d}")
np.atleast_1d(axes).ravel()[0].legend()
fig.suptitle(f"Clipped cells near the wing, time {time} ({tag})")
fig.tight_layout()
out = f"figures/cfd_section_overlay_{tag}.png"
fig.savefig(out, dpi=150)
print("saved", out)
