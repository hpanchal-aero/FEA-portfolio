"""
Stage 5 -- PyVista visual verification of the toe-fillet geometry.
Generates a coarse surface mesh (visualization only, not for solving)
from the verified BRep, then renders:
  1. Full solid overview
  2. A y=20 (mid-width) cross-section slice, zoomed to the toe region,
     to visually confirm the fillet is a smooth concave blend and not
     a malformed or misplaced surface.
"""

import gmsh
import pyvista as pv
import numpy as np
import os

INPUT_BREP = "mesh/full_bracket_study/frame_gusset_toefillet_check.brep"
STL_OUT = "mesh/full_bracket_study/toe_fillet_visual_check.stl"
FIG_DIR = "figures"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(os.path.dirname(STL_OUT), exist_ok=True)

# --- Generate a coarse surface mesh purely for visualization ---
gmsh.initialize()
gmsh.model.add("toe_fillet_visual")
gmsh.model.occ.importShapes(INPUT_BREP)
gmsh.model.occ.synchronize()

# Coarse but not so coarse it hides the 1mm fillet -- cap element size
gmsh.option.setNumber("Mesh.MeshSizeMax", 1.5)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
gmsh.model.mesh.generate(2)
gmsh.write(STL_OUT)
gmsh.finalize()

print(f"Wrote visualization surface mesh: {STL_OUT}")

# --- Load into PyVista ---
mesh = pv.read(STL_OUT)
print(f"Loaded mesh: {mesh.n_points} points, {mesh.n_cells} cells")

# --- Figure 1: full solid overview ---
pl = pv.Plotter(off_screen=True, window_size=(1400, 1000))
pl.add_mesh(mesh, color="lightsteelblue", show_edges=True, edge_color="gray", line_width=0.5)
pl.add_axes()
pl.camera_position = [(150, -150, 120), (30, 20, 15), (0, 0, 1)]
pl.add_text("Stage 5 - Frame + Gusset + Toe Fillet (overview)", font_size=12)
out1 = os.path.join(FIG_DIR, "stage5_toe_fillet_overview.png")
pl.screenshot(out1)
pl.close()
print(f"Saved: {out1}")

# --- Figure 2: mid-width slice at y=20, zoomed to toe region ---
slice_plane = mesh.slice(normal=(0, 1, 0), origin=(0, 20, 0))

pl2 = pv.Plotter(off_screen=True, window_size=(1400, 1000))
pl2.add_mesh(slice_plane, color="black", line_width=3, render_lines_as_tubes=True)
pl2.view_xz()
pl2.camera.parallel_projection = True
# Zoom to toe region: x in [33,43], z in [-3,10]
pl2.set_position((38, -100, 3.5))
pl2.set_focus((38, 0, 3.5))
pl2.set_viewup((0, 0, 1))
pl2.camera.parallel_scale = 6  # tight zoom in mm
pl2.add_text("Stage 5 - Toe fillet region, y=20mm slice", font_size=12)
out2 = os.path.join(FIG_DIR, "stage5_toe_fillet_slice_zoom.png")
pl2.screenshot(out2)
pl2.close()
print(f"Saved: {out2}")

# --- Figure 3: same slice, wider view for context (shows gusset + corner) ---
pl3 = pv.Plotter(off_screen=True, window_size=(1400, 1000))
pl3.add_mesh(slice_plane, color="black", line_width=2.5, render_lines_as_tubes=True)
pl3.view_xz()
pl3.camera.parallel_projection = True
pl3.set_position((45, -100, 20))
pl3.set_focus((45, 0, 20))
pl3.set_viewup((0, 0, 1))
pl3.camera.parallel_scale = 35
pl3.add_text("Stage 5 - Full corner/gusset region, y=20mm slice", font_size=12)
out3 = os.path.join(FIG_DIR, "stage5_toe_fillet_slice_context.png")
pl3.screenshot(out3)
pl3.close()
print(f"Saved: {out3}")

print("\nDone. Please inspect all three PNGs:")
print(f"  1. {out1} -- full solid sanity check")
print(f"  2. {out2} -- tight zoom on toe fillet, should show a smooth")
print(f"     concave arc blend at x~38, z~2, NOT a sharp corner or spike")
print(f"  3. {out3} -- wider context: R3 base fillet, corner-bite void,")
print(f"     and toe fillet all visible in one cross-section")
