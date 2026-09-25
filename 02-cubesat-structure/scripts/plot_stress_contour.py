"""
Project 02 -- von Mises stress contour from a CalculiX .inp/.frd pair.

Usage: /usr/bin/python3 scripts/plot_stress_contour.py <job_path_without_ext>

Reads *NODE / *ELEMENT (C3D10) from .inp (4 corner nodes per tet, same
simplification used in view_config_a_mesh.py), and the STRESS block from
.frd (fixed-width parser reused from extract_config_a.py), computes von
Mises per node, and renders an off-screen PyVista contour on the
UNDEFORMED shape (true displacements are ~1e-4mm -- not visually
meaningful at true scale; showing them undeformed avoids implying a scale
that was not applied). Color scale is auto-ranged per figure (job-to-job
comparison is NOT implied by this script alone -- ranges vary by orders
of magnitude across this project's cases).

Saves: figures/contours/<job>_vm_contour.png
"""
import os
import sys
import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_config_a import read_frd_block, von_mises  # noqa: E402


def read_inp(path):
    tags, xyz, conn, conn_full, mode = [], [], [], [], None
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
                vals = [int(v) for v in p[1:]]
                conn.append(vals[:4])
                conn_full.append(vals)
    tags = np.array(tags)
    if not np.array_equal(tags, np.arange(1, len(tags) + 1)):
        raise ValueError("node tags in .inp are not contiguous 1..N")
    return np.array(xyz), np.array(conn) - 1, np.array(conn_full) - 1


def main():
    if len(sys.argv) != 2:
        print("Usage: plot_stress_contour.py <job_path_without_ext>")
        sys.exit(1)
    job = sys.argv[1]
    inp, frd = job + ".inp", job + ".frd"
    stem = os.path.basename(job)
    os.makedirs("figures/contours", exist_ok=True)
    out_png = f"figures/contours/{stem}_vm_contour.png"

    xyz, conn, conn_full = read_inp(inp)
    n = len(xyz)
    coords_1based = np.vstack([np.zeros((1, 3)), xyz])

    S = read_frd_block(frd, "STRESS", 6, n)
    vm_1based = von_mises(S)
    vm = vm_1based[1:]  # drop the dummy row-0

    n_el = len(conn)
    cells = np.hstack([np.full((n_el, 1), 4), conn]).ravel()
    celltypes = np.full(n_el, pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, celltypes, xyz)
    grid["von_mises_MPa"] = vm

    surf = grid.extract_surface(progress_bar=False)

    p = pv.Plotter(off_screen=True, window_size=(1400, 1100))
    p.add_mesh(surf, scalars="von_mises_MPa", cmap="turbo", show_edges=False,
               scalar_bar_args={"title": "von Mises (MPa)"})
    p.add_axes()
    p.camera_position = "iso"
    p.add_text(f"{stem}\nvon Mises stress (undeformed shape, auto-scaled)",
               font_size=9)
    p.screenshot(out_png)
    p.close()

    connected = np.zeros(n, dtype=bool)
    connected[np.unique(conn_full)] = True
    vm_connected = vm[connected]
    n_excluded = int((~connected).sum())
    excl_note = f" ({n_excluded} unconnected node(s) excluded, e.g. *EQUATION reference node)" if n_excluded else ""
    print(f"{stem}: nodes={n} elements={n_el} "
          f"VM min={vm_connected.min():.4f} max={vm_connected.max():.4f} MPa{excl_note}  -> {out_png}")


if __name__ == "__main__":
    main()
