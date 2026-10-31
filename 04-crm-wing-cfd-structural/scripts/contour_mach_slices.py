#!/usr/bin/python3
"""Mach-number contours on spanwise slices (planes y = const) of the CFD internal field.

usage: contour_mach_slices.py <out.png> "label1=internal1.vtk" ["label2=internal2.vtk" ...]
Rows = cases, columns = y stations (6, 14, 22 m).  View along +y, x to the right, z up, parallel
projection, identical scale in every panel (half-height 7.5 m), each panel centred on the wing section.
Colour: Ma (coolwarm, 0.6..1.4, neutral at 1.0).  Thick black line = sonic line (Ma = 1.0), thin black
lines = Ma 1.2 and 1.4, heavy black outline = wing section, cut from the wing-patch VTK of the same time
(VTK/wing/wing_<t>.vtk).
"""
import os
import sys

import numpy as np
import pyvista as pv

STATIONS = (6.0, 14.0, 22.0)
CLIM = (0.6, 1.4)
SCALE = 7.5
pv.OFF_SCREEN = True

out = sys.argv[1]
cases = [a.split("=", 1) for a in sys.argv[2:]]
if not cases or any(len(c) != 2 for c in cases):
    sys.exit('usage: contour_mach_slices.py <out.png> "label=internal.vtk" ...')
nr, nc = len(cases), len(STATIONS)
pl = pv.Plotter(shape=(nr, nc), off_screen=True, window_size=(1000 * nc, 600 * nr), border=False)
pl.set_background("white")

for r, (label, path) in enumerate(cases):
    m = pv.read(path)
    base = os.path.dirname(path)
    t = os.path.basename(path)[:-4].rsplit("_", 1)[1]
    wing = pv.read(os.path.join(base, "wing", f"wing_{t}.vtk"))
    for c, y in enumerate(STATIONS):
        sl = m.slice(normal=(0.0, 1.0, 0.0), origin=(0.0, y, 0.0))
        wl = wing.slice(normal=(0.0, 1.0, 0.0), origin=(0.0, y, 0.0))
        b = wl.bounds
        xc, zc = 0.5 * (b[0] + b[1]), 0.5 * (b[4] + b[5])
        ma = np.asarray(sl.point_data["Ma"])
        print(f"{label} y={y:4.1f}: slice cells {sl.n_cells}  Ma max {ma.max():.3f}  "
              f"section x {b[0]:.2f}..{b[1]:.2f}  z {b[4]:.2f}..{b[5]:.2f}")
        pl.subplot(r, c)
        show_bar = (r == nr - 1 and c == nc - 1)
        sb = dict(title="Mach", vertical=True, height=0.7, position_x=0.90, position_y=0.15, color="black",
                  title_font_size=18, label_font_size=16, n_labels=5) if show_bar else None
        pl.add_mesh(sl, scalars="Ma", cmap="coolwarm", lighting=False, clim=CLIM,
                    show_scalar_bar=show_bar, scalar_bar_args=sb)
        son = sl.contour(isosurfaces=[1.0], scalars="Ma")
        if son.n_points:
            son.points = son.points + np.array([0.0, -0.05, 0.0]); pl.add_mesh(son, color="black", line_width=2.5)
        hi = sl.contour(isosurfaces=[1.2, 1.4], scalars="Ma")
        if hi.n_points:
            hi.points = hi.points + np.array([0.0, -0.05, 0.0]); pl.add_mesh(hi, color="black", line_width=0.8)
        wl.points = wl.points + np.array([0.0, -0.05, 0.0]); pl.add_mesh(wl, color="black", line_width=3.0)
        pl.add_text(f"{label} - y = {y:.0f} m", font_size=11, color="black", position="upper_left")
        pl.camera.position = (xc, y - 100.0, zc)
        pl.camera.focal_point = (xc, y, zc)
        pl.camera.up = (0.0, 0.0, 1.0)
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = SCALE
pl.screenshot(out)
print("saved", out)
