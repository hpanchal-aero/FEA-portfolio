"""
*TIE verification: solid block (C3D10) + shell plate (S6) tied via a
surface-based *TIE, independently meshed (different densities), no
coincident nodes assumed. Cantilever, tip load, compared against the
closed-form composite-beam reference (0.041112mm).

Block: 100(L) x 10(W) x 10(H) mm, fixed at z=0 (renaming length axis
to z for consistency with the wing's spanwise convention).
Shell: 100 x 10mm, 1.0mm thick, flush on the block's top face (y=10).
Tip load: 10N in -Y at the free end (z=100).

Usage: /usr/bin/python3 scripts/tie_verify_build.py
"""
import gmsh
import math
import os
import sys
import numpy as np
from netgen.occ import OCCGeometry, Box, Pnt

L, W, H_BLOCK, T_SHELL = 100.0, 10.0, 10.0, 1.0
E_MOD, NU, RHO = 73100.0, 0.33, 2780e-12
TIP_LOAD = 10.0
TOL = 1e-3

OUT_DIR = "mesh/tie_verify"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs("simulation", exist_ok=True)


def build_block_step():
    box = Box(Pnt(0, 0, 0), Pnt(W, H_BLOCK, L))
    step_path = os.path.join(OUT_DIR, "block.step")
    box.WriteStep(step_path)
    return step_path


def naca_dummy():
    pass  # not used here; placeholder removed below


def mesh_block(step_path, maxh):
    geo = OCCGeometry(step_path)
    mesh = geo.GenerateMesh(maxh=maxh)
    points = mesh.Points()
    n_pts = len(points)
    coords = np.zeros((n_pts + 1, 3))
    for i, p in enumerate(points, start=1):
        coords[i] = [p.p[0], p.p[1], p.p[2]]
    tets = []
    for el in mesh.Elements3D():
        v = [pt.nr for pt in el.vertices]
        tets.append(v)
    tets = np.array(tets)
    tets = tets[:, [0, 1, 3, 2]]
    return coords, tets, n_pts


def promote_c3d10(coords, tets, n_pts):
    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    midside_cache, extra_coords, quad_elems = {}, {}, []
    next_tag = n_pts + 1
    for t in tets:
        conn = list(t)
        for a, b in edge_pairs:
            na, nb = t[a], t[b]
            key = (min(na, nb), max(na, nb))
            if key not in midside_cache:
                extra_coords[next_tag] = (coords[na] + coords[nb]) / 2.0
                midside_cache[key] = next_tag
                next_tag += 1
            conn.append(midside_cache[key])
        quad_elems.append(conn)
    total_nodes = n_pts + len(midside_cache)
    all_coords = np.vstack([coords, np.array([extra_coords[k] for k in sorted(extra_coords)])])
    return all_coords, quad_elems, total_nodes


C3D10_FACES = {
    frozenset((0, 1, 2)): 1,
    frozenset((0, 1, 3)): 2,
    frozenset((1, 2, 3)): 3,
    frozenset((0, 2, 3)): 4,
}


def find_top_faces(quad_elems, all_coords, y_target):
    """For each C3D10 element, check each of its 4 possible corner-triples;
    if all 3 corners sit at y=y_target, that's a top face. Returns list of
    (elem_index_1based, face_id)."""
    faces = []
    for ei, conn in enumerate(quad_elems, start=1):
        corners = conn[:4]
        ys = [all_coords[c, 1] for c in corners]
        on_top = [i for i, y in enumerate(ys) if abs(y - y_target) < TOL]
        if len(on_top) == 3:
            face_id = C3D10_FACES[frozenset(on_top)]
            faces.append((ei, face_id))
    return faces


def build_shell(maxh, y0):
    """Flat W x L shell at height y0, independently gridded (structured,
    different density than the block's tet mesh -- deliberately not
    coincident with the block's surface nodes)."""
    nx = max(2, round(W / maxh))
    nz = max(2, round(L / maxh))
    xs = np.linspace(0, W, 2*nx + 1)
    zs = np.linspace(0, L, 2*nz + 1)
    coords = {}
    tag_grid = -np.ones((2*nx + 1, 2*nz + 1), dtype=int)
    n = 0
    for j in range(2*nz + 1):
        for i in range(2*nx + 1):
            if i % 2 == 1 and j % 2 == 1:
                continue
            n += 1
            tag_grid[i, j] = n
            coords[n] = (xs[i], y0, zs[j])
    elems = []
    for ez in range(nz):
        for ex in range(nx):
            i0, j0 = 2*ex, 2*ez
            c1, c2, c3, c4 = (tag_grid[i0, j0], tag_grid[i0+2, j0],
                              tag_grid[i0+2, j0+2], tag_grid[i0, j0+2])
            m12, m23 = tag_grid[i0+1, j0], tag_grid[i0+2, j0+1]
            m34, m41 = tag_grid[i0+1, j0+2], tag_grid[i0, j0+1]
            elems.append([c1, c2, c3, c4, m12, m23, m34, m41])
    return coords, elems, n


def main():
    print("=" * 70)
    print("TIE VERIFICATION: block + shell, cantilever tip load")
    print("=" * 70)

    step_path = build_block_step()
    print(f"Block STEP: {step_path}")

    coords, tets, n_pts = mesh_block(step_path, maxh=3.0)
    print(f"Block mesh: {len(tets)} tets, {n_pts} linear nodes")
    all_coords_b, quad_elems_b, total_nodes_b = promote_c3d10(coords, tets, n_pts)
    print(f"Block C3D10: {total_nodes_b} nodes, {len(quad_elems_b)} elements")

    top_faces = find_top_faces(quad_elems_b, all_coords_b, y_target=H_BLOCK)
    print(f"Block top-face (y={H_BLOCK}) facets found: {len(top_faces)}")
    if not top_faces:
        print("ERROR: no top faces found -- STOP."); sys.exit(1)

    shell_coords_local, shell_elems, n_shell = build_shell(maxh=4.0, y0=H_BLOCK)
    # offset shell node tags to avoid collision with block node tags
    offset = total_nodes_b
    shell_coords = {tg + offset: c for tg, c in shell_coords_local.items()}
    shell_elems_offset = [[v + offset for v in e] for e in shell_elems]
    print(f"Shell mesh: {n_shell} nodes (S6), {len(shell_elems_offset)} elements, "
          f"tag offset {offset}")

    fixed_block = [tg for tg in range(1, total_nodes_b + 1) if abs(all_coords_b[tg, 2] - 0.0) < TOL]
    fixed_shell = [tg for tg, c in shell_coords.items() if abs(c[2] - 0.0) < TOL]
    tip_shell = [tg for tg, c in shell_coords.items() if abs(c[2] - L) < TOL]
    print(f"Fixed block nodes (z=0): {len(fixed_block)}")
    print(f"Fixed shell nodes (z=0): {len(fixed_shell)}")
    print(f"Tip shell nodes (z={L}, load applied): {len(tip_shell)}")
    if not (fixed_block and fixed_shell and tip_shell):
        print("ERROR: an edge node set is empty -- STOP."); sys.exit(1)

    force_per_node = -TIP_LOAD / len(tip_shell)

    inp_out = "simulation/tie_verify.inp"
    with open(inp_out, "w") as f:
        f.write("*HEADING\nTIE verification: block+shell cantilever\n")

        f.write("*NODE\n")
        for tg in range(1, total_nodes_b + 1):
            c = all_coords_b[tg]
            f.write(f"{tg}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")
        for tg, c in shell_coords.items():
            f.write(f"{tg}, {c[0]:.10f}, {c[1]:.10f}, {c[2]:.10f}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=EBLOCK\n")
        for i, conn in enumerate(quad_elems_b, start=1):
            f.write(f"{i}, " + ", ".join(str(v) for v in conn) + "\n")

        f.write("*ELEMENT, TYPE=S8R, ELSET=ESHELL\n")
        shell_elem_start = len(quad_elems_b) + 1
        for i, conn in enumerate(shell_elems_offset, start=shell_elem_start):
            f.write(f"{i}, " + ", ".join(str(v) for v in conn) + "\n")

        f.write("*SURFACE, NAME=STOP, TYPE=ELEMENT\n")
        for ei, face_id in top_faces:
            f.write(f"{ei}, S{face_id}\n")
        f.write("*SURFACE, NAME=SSHELL, TYPE=ELEMENT\n")
        for i in range(shell_elem_start, shell_elem_start + len(shell_elems_offset)):
            f.write(f"{i}, SPOS\n")

        f.write("*NSET, NSET=NFIXB\n")
        for i in range(0, len(fixed_block), 10):
            f.write(", ".join(str(v) for v in fixed_block[i:i+10]) + "\n")
        f.write("*NSET, NSET=NFIXS\n")
        for i in range(0, len(fixed_shell), 10):
            f.write(", ".join(str(v) for v in fixed_shell[i:i+10]) + "\n")
        f.write("*NSET, NSET=NTIP\n")
        for i in range(0, len(tip_shell), 10):
            f.write(", ".join(str(v) for v in tip_shell[i:i+10]) + "\n")

        f.write("*MATERIAL, NAME=AL2024T3\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=EBLOCK, MATERIAL=AL2024T3\n")
        f.write("*SHELL SECTION, ELSET=ESHELL, MATERIAL=AL2024T3\n")
        f.write(f"{T_SHELL}\n")

        f.write("*TIE, NAME=T1, POSITION TOLERANCE=0.5\n")
        f.write("SSHELL, STOP\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\n")
        f.write("NFIXB, 1, 6, 0.0\n")
        f.write("NFIXS, 1, 6, 0.0\n")
        f.write("*CLOAD\n")
        for nd in tip_shell:
            f.write(f"{nd}, 2, {force_per_node:.8f}\n")
        f.write("*NODE FILE\nU\n*END STEP\n")

    print(f"\nWrote: {inp_out}")


if __name__ == "__main__":
    main()
