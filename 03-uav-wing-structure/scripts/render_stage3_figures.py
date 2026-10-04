#!/usr/bin/env python3
"""
Stage 3 3D figures (PyVista, off-screen). Each image uses a FRESH plotter (the first viewer version
saved two identical images because the second camera was set after the render).
  figures/stage3_mesh_overview.png   skin (blue, translucent), skin shared with the solid (orange),
                                     other solid faces (grey); t=1.0, solid 4 mm / skin 8 mm
  figures/stage3_mesh_root.png       close-up looking along +z at the root rib
  figures/stage3_buckling_mode1.png  critical buckling mode (factor 3.14), upper skin, plan view
                                     (span to the right, leading edge at the bottom), nodes colored
                                     by relative |U|
Needs the local viewer npz and buckle .frd (not in git; regenerate with the build scripts).
Run (conda env, pyvista):  python3 scripts/render_stage3_figures.py
"""
import os
import sys
import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_buckle_modes import read_modes

os.makedirs("figures", exist_ok=True)


def mesh_figures():
    d = np.load("mesh/stage3_skin/stage3_skin_s4_k8_t1_n3p8_viewer.npz")
    P, tris, cls = d["P"], d["tris"], d["tri_cls"]

    def poly(mask):
        t = tris[mask]
        return pv.PolyData(P, np.hstack([np.full((len(t), 1), 3), t]).ravel())

    solid, shared, free = poly(cls == 0), poly(cls == 1), poly(cls == 2)

    def build(free_opacity):
        p = pv.Plotter(off_screen=True, window_size=(1600, 800))
        p.set_background("white")
        p.add_mesh(solid, color="lightgrey", label="solid (rib faces, end faces)")
        p.add_mesh(shared, color="darkorange", label="skin shared with solid (rib rims, spar caps)")
        p.add_mesh(free, color="steelblue", opacity=free_opacity, show_edges=True, edge_color="navy",
                   line_width=0.3, label="free skin")
        p.add_legend(bcolor="white", size=(0.36, 0.12), loc="upper left")
        p.add_axes()
        return p

    p = build(0.35)
    p.camera_position = "iso"
    p.show(screenshot="figures/stage3_mesh_overview.png")
    p = build(0.25)
    p.camera_position = [(60.0, 90.0, -330.0), (60.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    p.show(screenshot="figures/stage3_mesh_root.png")


def mode_figure():
    job = "simulation/stage3_skin_s4_k8_t1_n3p8_buckle"
    ids, xyz, blocks = read_modes(job + ".frd")
    blk = blocks[1]                      # block 2 = buckling mode 1 (block 1 is the static preload)
    bid = np.array([b[0] for b in blk])
    U = np.array([b[1] for b in blk])
    o = np.argsort(ids)
    pos = np.searchsorted(ids[o], bid)
    P = xyz[o][pos]
    rel = np.linalg.norm(U, axis=1)
    rel = rel / rel.max()
    sel = (P[:, 1] > 0) & (P[:, 2] < 450.0)
    cloud = pv.PolyData(P[sel])
    cloud["relative |U|"] = rel[sel]
    p = pv.Plotter(off_screen=True, window_size=(1600, 700))
    p.set_background("white")
    p.add_mesh(cloud, scalars="relative |U|", cmap="viridis", point_size=4,
               render_points_as_spheres=True,
               scalar_bar_args=dict(title="relative displacement, mode 1", color="black"))
    p.add_text("Critical buckling mode (factor 3.14), upper skin, plan view: span to the right, "
               "leading edge at the bottom", color="black", font_size=9, position="upper_left")
    p.camera_position = [(60.0, 800.0, 200.0), (60.0, 0.0, 200.0), (1.0, 0.0, 0.0)]
    p.show(screenshot="figures/stage3_buckling_mode1.png")


if __name__ == "__main__":
    mesh_figures()
    mode_figure()
    for f in ("stage3_mesh_overview", "stage3_mesh_root", "stage3_buckling_mode1"):
        print(f, os.path.getsize(f"figures/{f}.png"), "bytes")
