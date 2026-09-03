"""
Project 01 - Aerospace Mounting Bracket
Generate stress/displacement contour plots from the finest (0.3mm)
converged hex mesh result, using PyVista to read the CalculiX .frd
file directly.

Outputs (saved to figures/):
  von_mises_contour.png   - full-body von Mises stress contour
  von_mises_root_zoom.png - zoomed view of the root region
  displacement_contour.png - deformed shape colored by |U|, exaggerated
"""

import pyvista as pv
import numpy as np
import os

CASE_DIR = "../simulation/convergence_hex/mesh_0.3mm"
FIG_DIR = "../figures"
os.makedirs(FIG_DIR, exist_ok=True)

# --- Parse mesh (nodes + hex elements) using the same validated logic ---
import sys
sys.path.insert(0, ".")
import hex_lib as lib

clean_inp = os.path.join(CASE_DIR, "mesh_clean.inp")
nodes, elements = lib.parse_clean_mesh(clean_inp)

node_ids_sorted = sorted(nodes.keys())
id_to_idx = {nid: i for i, nid in enumerate(node_ids_sorted)}
points = np.array([nodes[nid] for nid in node_ids_sorted])

# Build PyVista unstructured grid of quadratic hexahedra (VTK type 25)
cells = []
for eid, node_list in elements.items():
    idx20 = [id_to_idx[n] for n in node_list]
    cells.append(20)
    cells.extend(idx20)
cells = np.array(cells)
cell_types = np.full(len(elements), 25, dtype=np.uint8)  # VTK_QUADRATIC_HEXAHEDRON

grid = pv.UnstructuredGrid(cells, cell_types, points)

# --- Parse displacement and stress from .frd (same fixed-width parser) ---
frd_path = os.path.join(CASE_DIR, "analysis.frd")
with open(frd_path) as f:
    frd_lines = f.readlines()

disp = lib.parse_frd_block(frd_lines, "DISP", 3)
stress = lib.parse_frd_block(frd_lines, "STRESS", 6)

disp_arr = np.array([disp.get(nid, (0, 0, 0)) for nid in node_ids_sorted])
vm_arr = np.array([
    lib.von_mises(stress[nid]) if nid in stress else 0.0
    for nid in node_ids_sorted
])

grid["displacement"] = disp_arr
grid["von_mises"] = vm_arr
grid["disp_mag"] = np.linalg.norm(disp_arr, axis=1)

# --- Plot 1: full-body von Mises contour ---
pl = pv.Plotter(off_screen=True, window_size=(1600, 900))
pl.add_mesh(grid, scalars="von_mises", cmap="turbo", show_edges=False,
            scalar_bar_args={"title": "von Mises Stress (MPa)"})
pl.camera_position = [(120, -80, 60), (30, 20, 2), (0, 0, 1)]
pl.add_text("Von Mises Stress - Full Plate (0.3mm mesh)", font_size=14)
pl.screenshot(os.path.join(FIG_DIR, "von_mises_contour.png"))
pl.close()
print("Wrote von_mises_contour.png")

# --- Plot 2: zoomed root region (where peak bending stress occurs) ---
root_zone = grid.clip_box(bounds=(0, 15, 0, 40, 0, 4), invert=False)
pl = pv.Plotter(off_screen=True, window_size=(1600, 900))
pl.add_mesh(root_zone, scalars="von_mises", cmap="turbo", show_edges=True,
            scalar_bar_args={"title": "von Mises Stress (MPa)"})
pl.camera_position = [(30, -30, 20), (7, 20, 2), (0, 0, 1)]
pl.add_text("Von Mises Stress - Root Region Detail", font_size=14)
pl.screenshot(os.path.join(FIG_DIR, "von_mises_root_zoom.png"))
pl.close()
print("Wrote von_mises_root_zoom.png")

# --- Plot 3: deformed shape (exaggerated) colored by displacement magnitude ---
warp_factor = 20.0  # exaggerate for visibility (real tip deflection ~1mm)
warped = grid.warp_by_vector("displacement", factor=warp_factor)
pl = pv.Plotter(off_screen=True, window_size=(1600, 900))
pl.add_mesh(grid, color="lightgray", opacity=0.25, show_edges=False)
pl.add_mesh(warped, scalars="disp_mag", cmap="viridis", show_edges=False,
            scalar_bar_args={"title": "Displacement Magnitude (mm)"})
pl.camera_position = [(120, -80, 60), (30, 20, 2), (0, 0, 1)]
pl.add_text(f"Deformed Shape (displacement exaggerated {warp_factor}x)",
            font_size=14)
pl.screenshot(os.path.join(FIG_DIR, "displacement_contour.png"))
pl.close()
print("Wrote displacement_contour.png")

print(f"\nAll figures written to {os.path.abspath(FIG_DIR)}")
