"""
Project 02, Config A -- mesh visualization from a CalculiX .inp.

Usage:
    /usr/bin/python3 scripts/view_config_a_mesh.py <path/to/file.inp> [--interactive]

Reads *NODE and *ELEMENT (C3D10) from the .inp, uses the 4 corner nodes of each
tet to build a linear-tet PyVista grid, and:
  - saves the outer surface as STL (viewing only, not an analysis file)
  - saves 3 PNG figures: surface with edges, z=50 cross-section, corner zoom
  - prints a through-thickness node check for a mid-panel location

Default is off-screen (PNG only). --interactive opens a window afterwards.
"""

import os
import sys
import numpy as np
import pyvista as pv


def read_inp(path):
    tags, xyz, conn = [], [], []
    mode = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                if u.startswith("*NODE") and not u.startswith("*NODE FILE"):
                    mode = "node"
                elif u.startswith("*ELEMENT"):
                    mode = "elem"
                else:
                    mode = None
                continue
            if mode == "node":
                p = s.split(",")
                tags.append(int(p[0]))
                xyz.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                p = s.split(",")
                conn.append([int(v) for v in p[1:5]])   # 4 corner nodes only
    tags = np.array(tags)
    xyz = np.array(xyz)
    conn = np.array(conn)
    if not np.array_equal(tags, np.arange(1, len(tags) + 1)):
        print("ERROR: node tags are not contiguous 1..N; script assumes they are.")
        sys.exit(1)
    return xyz, conn - 1


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    interactive = "--interactive" in sys.argv
    if len(args) != 1:
        print("Usage: view_config_a_mesh.py <file.inp> [--interactive]")
        sys.exit(1)
    inp = args[0]
    stem = os.path.splitext(os.path.basename(inp))[0]
    os.makedirs("figures", exist_ok=True)
    stl_out = f"mesh/config_a_study/{stem}_surface.stl"

    xyz, conn = read_inp(inp)
    n_el = len(conn)
    cells = np.hstack([np.full((n_el, 1), 4), conn]).ravel()
    celltypes = np.full(n_el, pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, celltypes, xyz)

    surf = grid.extract_surface().triangulate()
    surf.save(stl_out)
    print(f"Read {inp}: {len(xyz)} nodes (incl. midside), {n_el} tets")
    print(f"Bounds: {grid.bounds}")
    print(f"Surface triangles: {surf.n_cells}")
    print(f"Wrote STL: {stl_out}")

    # ---- through-thickness check in the mid-region of the y=-50..-48 panel ----
    corner_ids = np.unique(conn)
    c = xyz[corner_ids]
    in_panel_span = (np.abs(c[:, 0]) < 20.0) & (c[:, 2] > 10.0) & (c[:, 2] < 90.0)
    on_faces = in_panel_span & ((np.abs(c[:, 1] + 50.0) < 1e-3) | (np.abs(c[:, 1] + 48.0) < 1e-3))
    interior = in_panel_span & (c[:, 1] > -49.99) & (c[:, 1] < -48.01)
    print(f"\nThrough-thickness check (panel y in [-50,-48], |x|<20, 10<z<90, corner nodes):")
    print(f"  nodes on the two panel faces: {int(on_faces.sum())}")
    print(f"  nodes strictly inside the 2 mm wall: {int(interior.sum())}")

    # ---- figure 1: outer surface with edges ----
    p = pv.Plotter(off_screen=True, window_size=(1400, 1000))
    p.add_mesh(surf, color="lightsteelblue", show_edges=True, edge_color="black",
               line_width=0.3)
    p.add_axes()
    p.camera_position = "iso"
    p.add_text(f"Config A outer surface, {n_el} tets", font_size=10)
    f1 = f"figures/{stem}_surface.png"
    p.screenshot(f1)
    p.close()

    # ---- figure 2: z=50 cross-section, full ----
    sl = grid.slice(normal="z", origin=(0, 0, 50.0))
    p = pv.Plotter(off_screen=True, window_size=(1400, 1400))
    p.add_mesh(sl, color="lightsteelblue", show_edges=True, edge_color="black",
               line_width=0.6)
    p.camera_position = "xy"
    p.add_text("z = 50 mm cross-section", font_size=10)
    f2 = f"figures/{stem}_section_z50.png"
    p.screenshot(f2)
    p.close()

    # ---- figure 3: corner zoom of the same section (rail + two panel ends) ----
    zoom = sl.clip_box(bounds=(-50, -25, -50, -25, 0, 100), invert=False)
    p = pv.Plotter(off_screen=True, window_size=(1400, 1400))
    p.add_mesh(zoom, color="lightsteelblue", show_edges=True, edge_color="black",
               line_width=1.0)
    p.camera_position = "xy"
    p.add_text("z = 50 mm section, corner (-50,-50) zoom", font_size=10)
    f3 = f"figures/{stem}_section_corner_zoom.png"
    p.screenshot(f3)
    p.close()

    print(f"\nWrote figures:\n  {f1}\n  {f2}\n  {f3}")

    if interactive:
        p = pv.Plotter()
        p.add_mesh(surf, color="lightsteelblue", show_edges=True, edge_color="black",
                   line_width=0.3)
        p.add_axes()
        p.show()


if __name__ == "__main__":
    main()
