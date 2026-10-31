#!/usr/bin/python3
"""Cp contours on the CFD wing surface, planform views (upper and lower surface).

usage: contour_cp_surface.py <out.png> "label[|note]=wing1.vtk" ["label2[|note2]=wing2.vtk" ...]

Cp = (p - pinf) / q with pinf = 21662.7 Pa, q = 10955.9 Pa.  Red = suction (negative Cp), blue = high pressure.
One column per case; row 1 = view from above (upper surface), row 2 = view from below (lower surface).
Faces are drawn flat with their own Cp (no smoothing).  Isolines (from point-averaged Cp) at Cp = -1.0, 0.0, 0.5.
Parallel projection, common colour scale (-1, 1).  An optional note after '|' is printed in the upper-surface panel.
"""
import sys

import numpy as np
import pyvista as pv

PINF, Q = 21662.7, 10955.9
CLIM = (-1.0, 1.0)
LEVELS = [-1.0, 0.0, 0.5]
pv.OFF_SCREEN = True

out = sys.argv[1]
cases = [a.split("=", 1) for a in sys.argv[2:]]
if not cases or any(len(c) != 2 for c in cases):
    sys.exit('usage: contour_cp_surface.py <out.png> "label[|note]=wing.vtk" ...')
n = len(cases)
pl = pv.Plotter(shape=(2, n), off_screen=True, window_size=(800 * n, 1500), border=False)
pl.set_background("white")

for j, (label, path) in enumerate(cases):
    title, _, note = label.partition("|")
    m = pv.read(path)
    m.cell_data["Cp"] = (np.asarray(m.cell_data["p"]) - PINF) / Q
    mp = m.cell_data_to_point_data()
    b = m.bounds
    xc, yc = 0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])
    print(f"{title}: faces {m.n_cells}  Cp min {m.cell_data['Cp'].min():.3f}  max {m.cell_data['Cp'].max():.3f}  "
          f"bounds x {b[0]:.2f}..{b[1]:.2f} y {b[2]:.2f}..{b[3]:.2f} z {b[4]:.2f}..{b[5]:.2f}")
    lines = mp.contour(isosurfaces=LEVELS, scalars="Cp")
    for i, (name, sgn) in enumerate((("upper surface (view from above)", 1.0), ("lower surface (view from below)", -1.0))):
        pl.subplot(i, j)
        sb = dict(title="Cp", vertical=True, height=0.6, position_x=0.90, position_y=0.2, color="black",
                  title_font_size=18, label_font_size=16, n_labels=5) if (j == n - 1 and i == 0) else None
        pl.add_mesh(m, scalars="Cp", preference="cell", lighting=False, cmap="RdBu", clim=CLIM,
                    show_scalar_bar=sb is not None, scalar_bar_args=sb)
        if lines.n_points:
            ln = lines.copy()
            ln.points = ln.points + np.array([0.0, 0.0, 0.03 * sgn])
            pl.add_mesh(ln, color="black", line_width=1.2)
        pl.add_text(f"{title} - {name}" + (chr(10) + note if (note and i == 0) else ""), font_size=11, color="black", position="upper_left")
        pl.camera.position = (xc, yc, 120.0 * sgn)
        pl.camera.focal_point = (xc, yc, 0.0)
        pl.camera.up = (0.0, 1.0 * sgn, 0.0)
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = 17.5
pl.screenshot(out)
print("saved", out)
