"""
Stage 5 -- PyVista von Mises contour visualization, Level 3 (converged) result.

Builds a VTK QUADRATIC_TETRA (cell type 24) unstructured grid directly from
the C3D10 .inp connectivity -- confirmed the edge_pairs order used in
mesh_netgen_level1_v2.py ((0,1)(1,2)(0,2)(0,3)(1,3)(2,3)) matches VTK's
native quadratic-tet node ordering, so no reordering is needed here.

v2 fix: full-bracket front/back views in v1 used hand-picked camera
coordinates that did not correctly frame the model (blank/cropped output).
Replaced with PyVista's view_vector() + reset_camera(), which auto-fits
the camera to the actual mesh bounds along a given direction -- removes
the need to know exact model extents in advance.

Produces:
  figures/stage5_full_bracket_vm_iso.png     -- isometric
  figures/stage5_full_bracket_vm_top.png     -- top-down (+z), flange A face-on
  figures/stage5_full_bracket_vm_front.png   -- front, viewing along -x
  figures/stage5_full_bracket_vm_back.png    -- back, viewing along +x (gusset side)
  figures/stage5_hole_vm_top.png             -- hole close-up, top-down (286 MPa peak)
  figures/stage5_hole_vm_angle.png           -- hole close-up, angled
  figures/stage5_toe_vm_top.png              -- toe close-up, top-down (ridge, not artifact)
  figures/stage5_toe_vm_angle.png            -- toe close-up, angled
"""

import numpy as np
import pyvista as pv
import os

INP_FILE = "simulation/stage5_netgen_L3_solve.inp"
FRD_FILE = "simulation/stage5_netgen_L3_solve.frd"
FIG_DIR = "figures"

TOE_POINT = np.array([38.0, 20.0, 2.0])
HOLE_POINT = np.array([20.0, 20.0, 2.0])

VTK_QUADRATIC_TETRA = 24

def parse_inp_nodes_elements(inp_path):
    nodes = {}
    elements = []
    mode = None
    with open(inp_path) as f:
        for line in f:
            s = line.strip()
            su = s.upper()
            if su == "*NODE":
                mode = "node"; continue
            if su.startswith("*ELEMENT") and "C3D10" in su:
                mode = "elem"; continue
            if su.startswith("*"):
                mode = None; continue
            if not s:
                continue
            parts = [p.strip() for p in s.split(",")]
            if mode == "node":
                tag = int(parts[0])
                nodes[tag] = [float(parts[1]), float(parts[2]), float(parts[3])]
            elif mode == "elem":
                conn = [int(x) for x in parts[1:]]
                elements.append(conn)
    return nodes, elements

def parse_frd_stress(frd_path):
    stress = {}
    with open(frd_path) as f:
        lines = f.readlines()
    mode = False
    for line in lines:
        if line.startswith(" -4") and "STRESS" in line:
            mode = True; continue
        if mode and line.startswith(" -1"):
            node = int(line[3:13])
            stress[node] = [
                float(line[13:25]), float(line[25:37]), float(line[37:49]),
                float(line[49:61]), float(line[61:73]), float(line[73:85]),
            ]
            continue
        if mode and not line.startswith(" -1") and not line.startswith(" -5"):
            mode = False
    return stress

def von_mises(s):
    s11, s22, s33, s12, s13, s23 = s
    return np.sqrt(0.5 * ((s11-s22)**2 + (s22-s33)**2 + (s33-s11)**2
                          + 6*(s12**2 + s13**2 + s23**2)))

def build_grid():
    print("Reading nodes/elements from .inp ...")
    nodes, elements = parse_inp_nodes_elements(INP_FILE)
    print(f"  {len(nodes)} nodes, {len(elements)} C3D10 elements")

    print("Reading STRESS block from .frd ...")
    stress = parse_frd_stress(FRD_FILE)
    print(f"  {len(stress)} nodes with stress data")

    tag_list = sorted(nodes.keys())
    tag_to_idx = {tag: i for i, tag in enumerate(tag_list)}
    points = np.array([nodes[t] for t in tag_list])

    n_missing = sum(1 for t in tag_list if t not in stress)
    if n_missing:
        print(f"  WARNING: {n_missing} nodes have no stress entry (set to 0)")

    vm = np.zeros(len(tag_list))
    for t in tag_list:
        if t in stress:
            vm[tag_to_idx[t]] = von_mises(stress[t])

    n_cells = len(elements)
    cells = np.zeros((n_cells, 11), dtype=np.int64)
    cells[:, 0] = 10
    for i, conn in enumerate(elements):
        cells[i, 1:11] = [tag_to_idx[n] for n in conn]
    cells = cells.flatten()
    cell_types = np.full(n_cells, VTK_QUADRATIC_TETRA, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells, cell_types, points)
    grid["von_mises"] = vm
    print(f"Peak VM (full model): {vm.max():.2f} MPa")
    print(f"Bounds: {grid.bounds}")
    return grid, points

def render_direction(mesh, out_path, view_vector, viewup=(0, 0, 1), clim=None,
                      title="von Mises (MPa)", window_size=(1400, 1000), zoom=1.0):
    """Auto-fits the camera to mesh bounds along a given view direction."""
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.add_mesh(mesh, scalars="von_mises", cmap="jet", show_edges=False,
                      clim=clim, scalar_bar_args={"title": title})
    plotter.add_axes()
    plotter.set_background("white")
    plotter.view_vector(view_vector, viewup=viewup)
    plotter.reset_camera()
    plotter.camera.zoom(zoom)
    plotter.screenshot(out_path)
    plotter.close()
    print(f"Wrote: {out_path}")

def render_iso(mesh, out_path, clim=None, title="von Mises (MPa)",
               window_size=(1400, 1000), zoom=1.2):
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.add_mesh(mesh, scalars="von_mises", cmap="jet", show_edges=False,
                      clim=clim, scalar_bar_args={"title": title})
    plotter.add_axes()
    plotter.set_background("white")
    plotter.camera_position = "iso"
    plotter.reset_camera()
    plotter.camera.zoom(zoom)
    plotter.screenshot(out_path)
    plotter.close()
    print(f"Wrote: {out_path}")

def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    grid, points = build_grid()
    global_clim = [0, float(grid["von_mises"].max())]

    # ---- Full bracket, multiple angles (auto-fit, bounds-safe) ----
    render_iso(grid, os.path.join(FIG_DIR, "stage5_full_bracket_vm_iso.png"),
               clim=global_clim)

    render_direction(grid, os.path.join(FIG_DIR, "stage5_full_bracket_vm_top.png"),
                      view_vector=(0, 0, 1), viewup=(0, 1, 0), clim=global_clim)

    render_direction(grid, os.path.join(FIG_DIR, "stage5_full_bracket_vm_front.png"),
                      view_vector=(1, 0, 0), viewup=(0, 0, 1), clim=global_clim)

    render_direction(grid, os.path.join(FIG_DIR, "stage5_full_bracket_vm_back.png"),
                      view_vector=(-1, 0, 0), viewup=(0, 0, 1), clim=global_clim)

    # ---- Hole close-up (governing peak, 286 MPa) -- unchanged, already good ----
    dist_hole = np.linalg.norm(points - HOLE_POINT, axis=1)
    hole_region = grid.extract_points(dist_hole < 8.0, adjacent_cells=True)
    hole_clim = [0, float(hole_region["von_mises"].max())]

    render_direction(hole_region, os.path.join(FIG_DIR, "stage5_hole_vm_top.png"),
                      view_vector=(0, 0, 1), viewup=(0, 1, 0), clim=hole_clim)

    render_direction(hole_region, os.path.join(FIG_DIR, "stage5_hole_vm_angle.png"),
                      view_vector=(1, -0.4, 0.6), viewup=(0, 0, 1), clim=hole_clim)

    # ---- Toe close-up (regularized ridge, ~230 MPa) -- unchanged, already good ----
    dist_toe = np.linalg.norm(points - TOE_POINT, axis=1)
    toe_region = grid.extract_points(dist_toe < 8.0, adjacent_cells=True)
    toe_clim = [0, float(toe_region["von_mises"].max())]

    render_direction(toe_region, os.path.join(FIG_DIR, "stage5_toe_vm_top.png"),
                      view_vector=(1, 0, 1), viewup=(0, 1, 0), clim=toe_clim)

    render_direction(toe_region, os.path.join(FIG_DIR, "stage5_toe_vm_angle.png"),
                      view_vector=(0.5, -1, 0.4), viewup=(0, 0, 1), clim=toe_clim)

if __name__ == "__main__":
    main()
