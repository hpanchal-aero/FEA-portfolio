"""
Stage 3 buckling-verification sub-stage: flat plate, simply supported,
uniaxial compression, S8R shell elements on a directly-generated
structured mesh (no Gmsh/Netgen for this standalone case -- see
development log for why).

Plate: 200 x 30 x 1.0mm (A x B x T), material Al 2024-T3.
Loaded edge: x=A, compression in -X. Reference load = closed-form
P_cr (scripts/buckle_verify_reference.py), so the CalculiX eigenvalue
should come out near 1.0 if the model is correct.

Usage: /usr/bin/python3 scripts/buckle_verify_build_mesh.py <nx> <ny>
  nx, ny: number of elements along A (200mm) and B (30mm)
"""
import os
import sys
import numpy as np

A, B, T = 200.0, 30.0, 1.0
E_MOD, NU, RHO = 73100.0, 0.33, 2780e-12
P_CR = 9017.4053  # N, closed-form reference (buckle_verify_reference.py)
TOL = 1e-6


def main():
    if len(sys.argv) != 3:
        print("Usage: buckle_verify_build_mesh.py <nx> <ny>"); sys.exit(1)
    nx, ny = int(sys.argv[1]), int(sys.argv[2])
    tag = f"buckle_verify_{nx}x{ny}"
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    # --- structured S8R grid: full point grid at half-element spacing,
    # omitting the (odd,odd) "center" points a serendipity element has no node at.
    xs = np.linspace(0.0, A, 2*nx + 1)
    ys = np.linspace(0.0, B, 2*ny + 1)
    tag_grid = -np.ones((2*nx + 1, 2*ny + 1), dtype=int)
    coords = {}
    n = 0
    for j in range(2*ny + 1):
        for i in range(2*nx + 1):
            if i % 2 == 1 and j % 2 == 1:
                continue  # no center node in a serendipity (8-node) element
            n += 1
            tag_grid[i, j] = n
            coords[n] = (xs[i], ys[j], 0.0)
    print(f"Grid: {nx}x{ny} elements, {n} nodes")

    elems = []
    for ey in range(ny):
        for ex in range(nx):
            i0, j0 = 2*ex, 2*ey
            c1, c2, c3, c4 = (tag_grid[i0, j0], tag_grid[i0+2, j0],
                              tag_grid[i0+2, j0+2], tag_grid[i0, j0+2])
            m12, m23 = tag_grid[i0+1, j0], tag_grid[i0+2, j0+1]
            m34, m41 = tag_grid[i0+1, j0+2], tag_grid[i0, j0+1]
            assert min(c1, c2, c3, c4, m12, m23, m34, m41) > 0
            elems.append([c1, c2, c3, c4, m12, m23, m34, m41])
    print(f"Elements: {len(elems)} (S8R)")

    def edge_nodes(i_fixed=None, j_fixed=None):
        out = []
        for tg, (x, y, z) in coords.items():
            if i_fixed is not None and abs(x - i_fixed) < TOL:
                out.append(tg)
            if j_fixed is not None and abs(y - j_fixed) < TOL:
                out.append(tg)
        return sorted(out)

    edge_x0 = edge_nodes(i_fixed=0.0)
    edge_xA = edge_nodes(i_fixed=A)
    edge_y0 = edge_nodes(j_fixed=0.0)
    edge_yB = edge_nodes(j_fixed=B)
    print(f"Edge node counts: x=0:{len(edge_x0)} x=A:{len(edge_xA)} "
          f"y=0:{len(edge_y0)} y=B:{len(edge_yB)}")

    # simply-supported: UZ=0 on all 4 edges (union, dedup)
    uz_nodes = sorted(set(edge_x0) | set(edge_xA) | set(edge_y0) | set(edge_yB))
    corner00 = tag_grid[0, 0]

    # consistent nodal loads on the x=A edge, quadratic-edge Simpson weights
    # (corner:mid:corner = 1:4:1 over each element-edge segment length),
    # representing the uniform line load q = P_CR/B = N_cr [N/mm]
    q = P_CR / B  # N/mm
    node_force = {}
    for ey in range(ny):
        j0 = 2*ey
        c_lo, c_hi, m = tag_grid[2*nx, j0], tag_grid[2*nx, j0+2], tag_grid[2*nx, j0+1]
        seg_len = ys[j0+2] - ys[j0]
        F = q * seg_len  # total force on this segment
        for tg, w in [(c_lo, 1/6), (m, 4/6), (c_hi, 1/6)]:
            node_force[tg] = node_force.get(tg, 0.0) + F * w
    total_F = sum(node_force.values())
    print(f"Consistent nodal load total: {total_F:.4f} N (target {P_CR:.4f})")
    if abs(total_F - P_CR) / P_CR > 1e-4:
        print("ERROR: consistent load does not sum to P_CR -- STOP."); sys.exit(1)

    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Stage 3 buckling verification, {nx}x{ny} S8R, plate {A}x{B}x{T}mm\n")

        f.write("*NODE\n")
        for tg in range(1, n + 1):
            x, y, z = coords[tg]
            f.write(f"{tg}, {x:.10f}, {y:.10f}, {z:.10f}\n")

        f.write("*ELEMENT, TYPE=S8R, ELSET=EALL\n")
        for i, c in enumerate(elems, start=1):
            f.write(f"{i}, " + ", ".join(str(v) for v in c) + "\n")

        f.write("*NSET, NSET=NUZ\n")
        for i in range(0, len(uz_nodes), 10):
            f.write(", ".join(str(v) for v in uz_nodes[i:i+10]) + "\n")
        f.write("*NSET, NSET=NX0\n")
        for i in range(0, len(edge_x0), 10):
            f.write(", ".join(str(v) for v in edge_x0[i:i+10]) + "\n")
        f.write("*NSET, NSET=NY0\n")
        for i in range(0, len(edge_y0), 10):
            f.write(", ".join(str(v) for v in edge_y0[i:i+10]) + "\n")

        f.write("*MATERIAL, NAME=AL2024T3\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SHELL SECTION, ELSET=EALL, MATERIAL=AL2024T3\n")
        f.write(f"{T}\n")

        f.write("*BOUNDARY\n")
        f.write("NUZ, 3, 3, 0.0\n")
        f.write("NX0, 1, 1, 0.0\n")
        f.write("NY0, 2, 2, 0.0\n")
        f.write(f"{corner00}, 6, 6, 0.0\n")

        f.write("*STEP\n*BUCKLE\n4\n")
        f.write("*CLOAD\n")
        for tg, F in node_force.items():
            f.write(f"{tg}, 1, {-F:.8f}\n")
        f.write("*NODE FILE\nU\n*END STEP\n")

    print(f"\nWrote: {inp_out}")


if __name__ == "__main__":
    main()
