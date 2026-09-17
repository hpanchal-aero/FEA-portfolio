"""
Stage 5 -- Build solver-ready CalculiX .inp from the TEMPORARY linear
(C3D4) baseline mesh. See mesh_full_bracket_baseline_linear.py header
for why this is linear-only and temporary.

C3D4 has only 4 nodes per element (no midside nodes), so there is no
node-ordering ambiguity between Gmsh and CalculiX conventions -- this
sidesteps that entire class of issue that affected the C3D10 attempts.
"""

import gmsh
import numpy as np
import os

MSH_IN = "mesh/full_bracket_study/stage5_level1_LINEAR_TEMP.msh"
INP_OUT = "simulation/stage5_level1_LINEAR_TEMP_solve.inp"

E_MOD = 71700.0
NU = 0.33
RHO = 2810e-9

FORCE_X = -235.44

ROOT_X_TARGET = 0.0
ROOT_TOL = 1e-3
TIP_Z_TARGET = 60.0
TIP_TOL = 1e-3

REF_NODE_COORD = (60.0, 20.0, 63.0)

def main():
    gmsh.initialize()
    gmsh.open(MSH_IN)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.array(node_tags, dtype=int)
    coords = np.array(node_coords).reshape(-1, 3)
    n_nodes = len(node_tags)
    print(f"Read {n_nodes} nodes from {MSH_IN}")

    root_mask = np.abs(coords[:, 0] - ROOT_X_TARGET) < ROOT_TOL
    root_nodes = node_tags[root_mask]
    print(f"Root nodes (x={ROOT_X_TARGET}): {len(root_nodes)}")

    tip_mask = np.abs(coords[:, 2] - TIP_Z_TARGET) < TIP_TOL
    tip_nodes = node_tags[tip_mask]
    print(f"Tip nodes (z={TIP_Z_TARGET}): {len(tip_nodes)}")
    tip_coords = coords[tip_mask]
    print(f"  tip x-range: [{tip_coords[:,0].min():.3f}, {tip_coords[:,0].max():.3f}]")
    print(f"  tip y-range: [{tip_coords[:,1].min():.3f}, {tip_coords[:,1].max():.3f}]")

    if len(root_nodes) == 0 or len(tip_nodes) == 0:
        print("ERROR: root or tip node set empty.")
        gmsh.finalize()
        return

    elem_types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(3)
    tet_idx = None
    for i, et in enumerate(elem_types):
        if et == 4:  # 4-node tetrahedron
            tet_idx = i
            break
    if tet_idx is None:
        print("ERROR: no 4-node tet elements found.")
        gmsh.finalize()
        return

    etags = np.array(elem_tags_list[tet_idx], dtype=int)
    enodes = np.array(elem_node_tags_list[tet_idx], dtype=int).reshape(len(etags), 4)
    print(f"Elements (C3D4): {len(etags)}")
    # No reordering needed -- 4-node tet has no ambiguous midside sequence.

    os.makedirs(os.path.dirname(INP_OUT), exist_ok=True)
    with open(INP_OUT, "w") as f:
        f.write("*HEADING\n")
        f.write("Stage 5 -- Full L-bracket integration, TEMPORARY LINEAR (C3D4) baseline\n")
        f.write("Al 7075-T6, F=235.44N in -x at Flange B tip, root encastre\n")
        f.write("NOTE: C3D4 elements are known to be overly stiff in bending.\n")
        f.write("This is a temporary sanity-check mesh, NOT the final Stage 5 formulation.\n")

        f.write("*NODE\n")
        for tag, c in zip(node_tags, coords):
            f.write(f"{tag}, {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}\n")

        f.write("*ELEMENT, TYPE=C3D4, ELSET=EBRACKET\n")
        for tag, nn in zip(etags, enodes):
            f.write(f"{tag}, " + ", ".join(str(n) for n in nn) + "\n")

        f.write("*NSET, NSET=NROOT\n")
        for i in range(0, len(root_nodes), 10):
            f.write(", ".join(str(n) for n in root_nodes[i:i+10]) + "\n")

        f.write("*NSET, NSET=NTIP\n")
        for i in range(0, len(tip_nodes), 10):
            f.write(", ".join(str(n) for n in tip_nodes[i:i+10]) + "\n")

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
        f.write("1, 3\n")  # explicit DOF list: solid elements have only 3 translational DOF/node

        f.write("*CLOAD\n")
        f.write(f"{ref_tag}, 1, {FORCE_X}\n")

        f.write("*NODE FILE\n")
        f.write("U, RF\n")
        f.write("*EL FILE\n")
        f.write("S, E\n")
        f.write("*END STEP\n")

    print(f"\nWrote: {INP_OUT}")
    print(f"Reference node: {ref_tag} at {REF_NODE_COORD}")
    print(f"Root nodes fixed: {len(root_nodes)}")
    print(f"Tip nodes coupled: {len(tip_nodes)}")
    print(f"DOF estimate: {3*n_nodes} (unconstrained)")

    gmsh.finalize()

if __name__ == "__main__":
    main()
