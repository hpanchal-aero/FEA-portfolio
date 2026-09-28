"""Quick PyVista render of the Stage 2 preview STL -- visualization only."""
import pyvista as pv
import os

STL_PATH = "mesh/stage2_ribs/stage2_ribs_PREVIEW.stl"
OUT_PNG = "figures/stage2_ribs_preview.png"
os.makedirs("figures", exist_ok=True)

mesh = pv.read(STL_PATH)

pl = pv.Plotter(off_screen=True, window_size=(1600, 1000))
pl.add_mesh(mesh, color="lightsteelblue", show_edges=True, edge_color="gray", line_width=0.3)
pl.add_axes()
pl.camera_position = "iso"
pl.camera.azimuth = 30
pl.camera.elevation = 20
pl.reset_camera()
pl.screenshot(OUT_PNG)
print(f"Saved: {OUT_PNG}")

pl2 = pv.Plotter(off_screen=True, window_size=(1600, 800))
pl2.add_mesh(mesh, color="lightsteelblue", show_edges=True, edge_color="gray", line_width=0.3)
pl2.view_xy()
pl2.reset_camera()
pl2.screenshot("figures/stage2_ribs_preview_top.png")
print("Saved: figures/stage2_ribs_preview_top.png")

pl3 = pv.Plotter(off_screen=True, window_size=(1600, 800))
pl3.add_mesh(mesh, color="lightsteelblue", show_edges=True, edge_color="gray", line_width=0.3)
pl3.view_yz()
pl3.reset_camera()
pl3.screenshot("figures/stage2_ribs_preview_end.png")
print("Saved: figures/stage2_ribs_preview_end.png")
