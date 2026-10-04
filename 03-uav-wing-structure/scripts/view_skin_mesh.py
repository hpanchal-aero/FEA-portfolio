#!/usr/bin/env python3
"""
Viewer for the classified skin+solid surface mesh (PyVista).

  grey   = solid boundary that is not skin (rib flat faces, end faces)
  orange = skin SHARED with the solid (rib rims, spar cap tops)
  blue   = free skin (shell only)

Interactive (WSLg window, opacity slider):
  python3 scripts/view_skin_mesh.py mesh/stage3_skin/probe_maxh12_classified.npz
Image mode (off-screen; writes <png> = full wing, and <png stem>_root.png = root close-up):
  python3 scripts/view_skin_mesh.py mesh/stage3_skin/probe_maxh12_classified.npz figures/skin_probe.png

Run inside the conda env (ccx223).
"""
import sys
import numpy as np
import pyvista as pv

d = np.load(sys.argv[1])
P, tris, cls = d["P"], d["tris"], d["tri_cls"]
png = sys.argv[2] if len(sys.argv) > 2 else None


def poly(mask):
    t = tris[mask]
    faces = np.hstack([np.full((len(t), 1), 3), t]).ravel()
    return pv.PolyData(P, faces)


solid, shared, free = poly(cls == 0), poly(cls == 1), poly(cls == 2)
print(f"triangles: solid-nonskin {int((cls == 0).sum())}, shared skin {int((cls == 1).sum())}, "
      f"free skin {int((cls == 2).sum())}")

p = pv.Plotter(off_screen=png is not None, window_size=(1600, 800))
p.set_background("white")
p.add_mesh(solid, color="lightgrey", show_edges=False, label="solid (non-skin faces)")
p.add_mesh(shared, color="darkorange", show_edges=False, label="skin shared with solid")
free_actor = p.add_mesh(free, color="steelblue", opacity=0.35, show_edges=True,
                        edge_color="navy", line_width=0.3, label="free skin")
p.add_legend(bcolor="white", size=(0.25, 0.12))
p.add_axes()
p.camera_position = "iso"


def set_opacity(v):
    free_actor.GetProperty().SetOpacity(float(v))


if png is None:
    p.add_slider_widget(set_opacity, [0.0, 1.0], value=0.35, title="free-skin opacity",
                        pointa=(0.05, 0.08), pointb=(0.30, 0.08), style="modern")
    p.show()
else:
    p.show(screenshot=png, auto_close=False)
    print(f"saved {png}")
    p.camera_position = [(260.0, 170.0, -140.0), (30.0, 0.0, 90.0), (0.0, 1.0, 0.0)]
    p.reset_camera_clipping_range()
    root_png = png[:-4] + "_root.png" if png.endswith(".png") else png + "_root.png"
    p.screenshot(root_png)
    print(f"saved {root_png}")
    p.close()
