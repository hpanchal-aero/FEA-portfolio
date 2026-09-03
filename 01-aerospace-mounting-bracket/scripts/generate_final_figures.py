"""
Project 01 - Aerospace Mounting Bracket
Final 3D visualization set for the verification plate, from the
converged 0.3mm result. Uses PyVista, reads .frd directly.

Figures (each chosen for genuine diagnostic/physical insight, not
decoration):
  1. von_mises_contour.png     - already generated, reused
  2. von_mises_root_zoom.png   - already generated, reused
  3. displacement_contour.png  - already generated, reused
  4. sigma_xx_contour.png      - NEW: shows the primary verification
                                  quantity spatially, ties directly
                                  to the sigma_xx table/plot
  5. mesh_visualization.png    - NEW: shows the actual structured
                                  hex mesh density/quality used for
                                  this converged result
"""

import os
import sys
import numpy as np
import pyvista as pv

sys.path.insert(0, ".")
import hex_lib as lib

CASE_DIR = "../simulation/convergence_hex/mesh_0.3mm"
FIG_DIR = "../figures"
os.makedirs(FIG_DIR, exist_ok=True)

clean_inp = os.path.join(CASE_DIR, "mesh_clean.inp")
nodes, elements = lib.parse_clean_mesh(clean_inp)

node_ids_sorted = sorted(nodes.keys())
id_to_idx = {nid: i for i, nid in enumerate(node_ids_sorted)}
points = np.array([nodes[nid] for nid in node_ids_sorted])

cells = []
for eid, node_list in elements.items():
    idx20 = [id_to_idx[n] for n in node_list]
    cells.append(20)
    cells.extend(idx20)
cells = np.array(cells)
cell_types = np.full(len(elements), 25, dtype=np.uint8)
grid = pv.UnstructuredGrid(cells, cell_types, points)

frd_path = os.path.join(CASE_DIR, "analysis.frd")
with open(frd_path) as f:
    frd_lines = f.readlines()
stress = lib.parse_frd_block(frd_lines, "STRESS", 6)

sxx_arr = np.array([stress[nid][0] if nid in stress else 0.0
                     for nid in node_ids_sorted])
grid["sigma_xx"] = sxx_arr

# --- Figure: sigma_xx contour (the primary verification quantity) ---
pl = pv.Plotter(off_screen=True, window_size=(1600, 900))
pl.add_mesh(grid, scalars="sigma_xx", cmap="coolwarm", show_edges=False,
            scalar_bar_args={"title": "Sigma_xx (MPa)"})
pl.camera_position = [(120, -80, 60), (30, 20, 2), (0, 0, 1)]
pl.add_text("Longitudinal Bending Stress (sigma_xx) - 0.3mm mesh",
            font_size=14)
pl.screenshot(os.path.join(FIG_DIR, "sigma_xx_contour.png"))
pl.close()
print("Wrote sigma_xx_contour.png")

# --- Figure: mesh visualization (use a coarser level for visibility -
#     0.3mm mesh has 344k elements and would be solid black wireframe;
#     using 1.3mm level here purely for VISUAL clarity of structure,
#     NOT as a result - clearly labeled as such) ---
COARSE_CASE = "../simulation/convergence_hex/mesh_1.3mm"
coarse_clean = os.path.join(COARSE_CASE, "mesh_clean.inp")
c_nodes, c_elements = lib.parse_clean_mesh(coarse_clean)

c_ids_sorted = sorted(c_nodes.keys())
c_id_to_idx = {nid: i for i, nid in enumerate(c_ids_sorted)}
c_points = np.array([c_nodes[nid] for nid in c_ids_sorted])

c_cells = []
for eid, node_list in c_elements.items():
    idx20 = [c_id_to_idx[n] for n in node_list]
    c_cells.append(20)
    c_cells.extend(idx20)
c_cells = np.array(c_cells)
c_cell_types = np.full(len(c_elements), 25, dtype=np.uint8)
c_grid = pv.UnstructuredGrid(c_cells, c_cell_types, c_points)

pl = pv.Plotter(off_screen=True, window_size=(1600, 900))
pl.add_mesh(c_grid, color="lightblue", show_edges=True,
            edge_color="black", line_width=1)
pl.camera_position = [(120, -80, 60), (30, 20, 2), (0, 0, 1)]
pl.add_text("Structured Hex Mesh (1.3mm level shown for visual clarity; "
            "0.3mm mesh used for reported results)", font_size=12)
pl.screenshot(os.path.join(FIG_DIR, "mesh_visualization.png"))
pl.close()
print("Wrote mesh_visualization.png")

print(f"\nAll figures saved to {os.path.abspath(FIG_DIR)}")
print("\nReused from earlier stage (already generated, still valid):")
print("  von_mises_contour.png, von_mises_root_zoom.png, "
      "displacement_contour.png")
