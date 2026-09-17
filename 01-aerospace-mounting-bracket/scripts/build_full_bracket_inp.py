"""
Stage 5 -- Build the solver-ready CalculiX .inp from the verified,
quality-checked Level 1 mesh (stage5_level1_v3.msh).

Adds: material (Al 7075-T6), solid section, root encastre BC,
kinematic-coupled tip load (F=235.44N in -x), output requests,
and *STATIC step -- on top of the raw node/element data already
exported by Gmsh.

Node-set identification is coordinate-based (tolerance-matched),
consistent with the identification approach used throughout this
project rather than relying on Gmsh's own (unreliable/renumbered)
entity groupings for BC purposes.
"""

import gmsh
import numpy as np
import os

MSH_IN = "mesh/full_bracket_study/stage5_level1_v3.msh"
INP_OUT = "simulation/stage5_level1_v3_solve.inp"

E_MOD = 71700.0      # MPa
NU = 0.33
RHO = 2810e-9        # t/mm^3 (consistent mm-N-t-s-MPa unit system, per prior stages)

FORCE_X = -235.44     # N, applied in -x at Flange B tip

ROOT_X_TARGET = 0.0
ROOT_TOL = 1e-3

TIP_Z_TARGET = 60.0
TIP_TOL = 1e-3

REF_NODE_COORD = (60.0, 20.0, 63.0)  # off-mesh, ahead of tip face for coupling

def main():
    gmsh.initialize()
    gmsh.open(MSH_IN)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.array(node_tags, dtype=int)
    coords = np.array(node_coords).reshape(-1, 3)
    n_nodes = len(node_tags)
    print(f"Read {n_nodes} nodes from {MSH_IN}")

    # --- Identify root nodes (x approx 0) ---
    root_mask = np.abs(coords[:, 0] - ROOT_X_TARGET) < ROOT_TOL
    root_nodes = node_tags[root_mask]
    print(f"Root nodes (x={ROOT_X_TARGET}): {len(root_nodes)}")
    if len(root_nodes) == 0:
        print("ERROR: no root nodes found -- check coordinate convention.")
        gmsh.finalize()
        return

    # --- Identify tip nodes (z approx 60) ---
    tip_mask = np.abs(coords[:, 2] - TIP_Z_TARGET) < TIP_TOL
    tip_nodes = node_tags[tip_mask]
    print(f"Tip nodes (z={TIP_Z_TARGET}): {len(tip_nodes)}")
    if len(tip_nodes) == 0:
        print("ERROR: no tip nodes found -- check coordinate convention.")
        gmsh.finalize()
        return

    # Sanity: tip nodes should span the Flange B cross-section (x=58..62, y=0..40)
    tip_coords = coords[tip_mask]
    print(f"  tip node x-range: [{tip_coords[:,0].min():.3f}, {tip_coords[:,0].max():.3f}] "
          f"(expect ~58 to 62)")
    print(f"  tip node y-range: [{tip_coords[:,1].min():.3f}, {tip_coords[:,1].max():.3f}] "
          f"(expect ~0 to 40)")

    # --- Element connectivity (C3D10, type 11 = 10-node tet in Gmsh) ---
    elem_types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(3)
    tet_idx = None
    for i, et in enumerate(elem_types):
        if et == 11:  # 10-node tetrahedron
            tet_idx = i
            break
    if tet_idx is None:
        print("ERROR: no 10-node tet elements found.")
        gmsh.finalize()
        return

    etags = np.array(elem_tags_list[tet_idx], dtype=int)
    enodes = np.array(elem_node_tags_list[tet_idx], dtype=int).reshape(len(etags), 10)
    print(f"Elements (C3D10): {len(etags)}")

    # Gmsh node order for 10-node tet vs Abaqus/CalculiX C3D10 node order differ
    # in the mid-side node sequence. Gmsh: [n1,n2,n3,n4, n12,n13,n14,n23,n24,n34]
    # Abaqus/CCX C3D10: [n1,n2,n3,n4, n12,n23,n13,n14,n24,n34]
    # Reorder columns 4-9 (0-indexed) accordingly: gmsh idx -> ccx position
    # gmsh: [0,1,2,3, 4(12),5(13),6(14),7(23),8(24),9(34)]
    # ccx:  [0,1,2,3, 12,23,13,14,24,34] -> from gmsh indices [4,7,5,6,8,9]
    reorder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  # CORRECTED: gmsh and CCX C3D10 mid-edge order already match; no permutation needed
    enodes_ccx = enodes[:, reorder]

    os.makedirs(os.path.dirname(INP_OUT), exist_ok=True)
    with open(INP_OUT, "w") as f:
        f.write("*HEADING\n")
        f.write("Stage 5 -- Full L-bracket integration, Level 1 mesh\n")
        f.write("Al 7075-T6, F=235.44N in -x at Flange B tip, root encastre\n")
        f.write("Verification/consistency study (see README) -- NOT closed-form validation\n")

        f.write("*NODE\n")
        for tag, c in zip(node_tags, coords):
            f.write(f"{tag}, {c[0]!r}, {c[1]!r}, {c[2]!r}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=EBRACKET\n")
        for tag, nn in zip(etags, enodes_ccx):
            f.write(f"{tag}, " + ", ".join(str(n) for n in nn) + "\n")

        f.write("*NSET, NSET=NROOT\n")
        for i in range(0, len(root_nodes), 10):
            f.write(", ".join(str(n) for n in root_nodes[i:i+10]) + "\n")

        f.write("*NSET, NSET=NTIP\n")
        for i in range(0, len(tip_nodes), 10):
            f.write(", ".join(str(n) for n in tip_nodes[i:i+10]) + "\n")

        # Reference node for kinematic coupling -- pick an unused high tag
        ref_tag = int(node_tags.max()) + 1
        f.write("*NODE, NSET=NREF\n")
        f.write(f"{ref_tag}, {REF_NODE_COORD[0]:.6f}, {REF_NODE_COORD[1]:.6f}, "
                f"{REF_NODE_COORD[2]:.6f}\n")

        f.write("*SURFACE, NAME=STIP, TYPE=NODE\n")
        f.write("NTIP\n")

        f.write("*MATERIAL, NAME=AL7075T6\n")
        f.write("*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n")
        f.write("*DENSITY\n")
        f.write(f"{RHO}\n")

        f.write("*SOLID SECTION, ELSET=EBRACKET, MATERIAL=AL7075T6\n")

        f.write("*STEP\n")
        f.write("*STATIC\n")
        f.write("*BOUNDARY\n")
        f.write("NROOT, 1, 3, 0.0\n")

        f.write("*COUPLING, CONSTRAINT NAME=CTIP, REF NODE=" + str(ref_tag) +
                ", SURFACE=STIP\n")
        f.write("*KINEMATIC\n")

        f.write("*CLOAD\n")
        f.write(f"{ref_tag}, 1, {FORCE_X}\n")

        f.write("*NODE FILE\n")
        f.write("U, RF\n")
        f.write("*EL FILE\n")
        f.write("S, E\n")
        f.write("*END STEP\n")

    print(f"\nWrote solver-ready input: {INP_OUT}")
    print(f"Reference node for tip load: {ref_tag} at {REF_NODE_COORD}")
    print(f"Root nodes fixed (1,3): {len(root_nodes)}")
    print(f"Tip nodes coupled to ref node: {len(tip_nodes)}")
    print(f"\nTotal DOF estimate: {3 * n_nodes} (unconstrained), "
          f"minus {3 * len(root_nodes)} constrained at root")

    gmsh.finalize()

if __name__ == "__main__":
    main()
