"""
Stage 6 -- PyVista von Mises visualization of the selected design (lhs_06),
at its finest verified mesh level (L2: 198,510 nodes).

Reuses the exact approach validated in Stage 5's visualize_stage5_stress.py
(VTK_QUADRATIC_TETRA=24 cell type, view_vector auto-fit camera framing).
"""

import numpy as np
import pyvista as pv
import os

INP_FILE = "simulation/stage6_parametric/stage6_lhs_06_L2_mesh.inp"
FRD_FILE = "simulation/stage6_parametric/stage6_lhs_06_L2_mesh.frd"
FIG_DIR = "figures"

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


def render_direction(mesh, out_path, view_vector, viewup=(0, 0, 1), clim=None,
                      title="von Mises (MPa)", window_size=(1400, 1000), zoom=1.0):
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

    print("Reading nodes/elements from .inp ...")
    nodes, elements = parse_inp_nodes_elements(INP_FILE)
    print(f"  {len(nodes)} nodes, {len(elements)} C3D10 elements")

    print("Reading STRESS block from .frd ...")
    stress = parse_frd_stress(FRD_FILE)
    print(f"  {len(stress)} nodes with stress data")

    tag_list = sorted(nodes.keys())
    tag_to_idx = {tag: i for i, tag in enumerate(tag_list)}
    points = np.array([nodes[t] for t in tag_list])

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
    print(f"Peak VM: {vm.max():.2f} MPa")

    clim = [0, float(vm.max())]

    render_iso(grid, os.path.join(FIG_DIR, "stage6_lhs06_vm_iso.png"), clim=clim)
    render_direction(grid, os.path.join(FIG_DIR, "stage6_lhs06_vm_top.png"),
                      view_vector=(0, 0, 1), viewup=(0, 1, 0), clim=clim)
    render_direction(grid, os.path.join(FIG_DIR, "stage6_lhs06_vm_front.png"),
                      view_vector=(1, 0, 0), viewup=(0, 0, 1), clim=clim)
    render_direction(grid, os.path.join(FIG_DIR, "stage6_lhs06_vm_back.png"),
                      view_vector=(-1, 0, 0), viewup=(0, 0, 1), clim=clim)


if __name__ == "__main__":
    main()
