"""
Neighbor patch test, v2: properly constrained.

v1 fixed only 3 arbitrary nodes (9 DOF) on an 83-element, 210-node
patch -- grossly under-constrained, allowing the patch to behave as
a near-mechanism and produce spurious large local distortions/Jacobian
failures at essentially random elements unrelated to the real defect.
Confirmed: target element 170300 passed in v1's patch test while a
different, wider set of elements failed -- inconsistent with a valid
per-element geometric check, and traced to inadequate BCs, not a
data-extraction bug (node/element line locations were verified correct
via block-position-aware parsing).

v2 fix: identify every "boundary" node of the patch (any node
belonging to an element OUTSIDE the 83-element patch, i.e. was
supported by the surrounding 157k-element model in the real mesh)
and fix it fully (0 displacement, all 3 DOF) -- approximating "held
by the stiffer surrounding material," which is what those nodes
effectively experience in the real full-model context. Only nodes
fully interior to the patch (touched ONLY by patch elements) are
left free.

Note: since Jacobian evaluation happens during element stiffness
assembly regardless of the eventual displacement solution, a load is
included but its exact magnitude/location is not critical to whether
the Jacobian error reproduces.
"""

import numpy as np

INP_FILE = "simulation/stage5_level1_v3_solve.inp"
OUT_INP = "simulation/neighbor_patch_test_v2.inp"
TARGET_ELEM = 170300

def main():
    nodes = {}
    all_elements = {}

    with open(INP_FILE) as f:
        lines = f.readlines()

    mode = None
    for line in lines:
        s = line.strip()
        su = s.upper()
        if su == "*NODE":
            mode = "node"
            continue
        if su.startswith("*ELEMENT"):
            mode = "elem"
            continue
        if su.startswith("*"):
            mode = None
            continue
        if mode == "node" and s:
            parts = s.split(",")
            nodes[int(parts[0])] = [p.strip() for p in parts[1:4]]
        elif mode == "elem" and s:
            parts = [p.strip() for p in s.split(",")]
            all_elements[int(parts[0])] = [int(p) for p in parts[1:11]]

    print(f"Parsed {len(nodes)} nodes, {len(all_elements)} elements from full model")

    target_nodes = set(all_elements[TARGET_ELEM])
    patch_elems = {TARGET_ELEM: all_elements[TARGET_ELEM]}
    for etag, enodes in all_elements.items():
        if etag == TARGET_ELEM:
            continue
        if target_nodes & set(enodes):
            patch_elems[etag] = enodes

    print(f"Patch size: {len(patch_elems)} elements")

    patch_node_set = set()
    for enodes in patch_elems.values():
        patch_node_set.update(enodes)
    print(f"Patch node count: {len(patch_node_set)}")

    # --- Find which patch nodes ALSO belong to elements OUTSIDE the patch ---
    # These are the "boundary" nodes -- in the real model, they're supported
    # by the surrounding material. We fix them to approximate that support.
    outside_elem_ids = set(all_elements.keys()) - set(patch_elems.keys())
    boundary_nodes = set()
    for etag in outside_elem_ids:
        enodes = all_elements[etag]
        shared = patch_node_set & set(enodes)
        if shared:
            boundary_nodes.update(shared)

    interior_nodes = patch_node_set - boundary_nodes
    print(f"\nPatch boundary nodes (touch elements outside patch, will be FIXED): "
          f"{len(boundary_nodes)}")
    print(f"Patch interior nodes (touched ONLY by patch elements, left FREE): "
          f"{len(interior_nodes)}")

    if len(interior_nodes) == 0:
        print("WARNING: no interior nodes -- entire patch is boundary. "
              "This would over-constrain the patch (no meaningful test possible).")

    patch_node_list = sorted(patch_node_set)

    with open(OUT_INP, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Neighbor patch test v2 -- element {TARGET_ELEM} + neighbors, "
                f"properly constrained boundary\n")
        f.write("*NODE\n")
        for n in patch_node_list:
            c = nodes[n]
            f.write(f"{n}, {c[0]}, {c[1]}, {c[2]}\n")
        f.write("*ELEMENT, TYPE=C3D10, ELSET=EPATCH\n")
        for etag, enodes in patch_elems.items():
            f.write(f"{etag}, " + ", ".join(str(n) for n in enodes) + "\n")
        f.write("*MATERIAL, NAME=AL7075T6\n*ELASTIC\n71700.0, 0.33\n")
        f.write("*SOLID SECTION, ELSET=EPATCH, MATERIAL=AL7075T6\n")

        f.write("*NSET, NSET=NBOUND\n")
        bnd_list = sorted(boundary_nodes)
        for i in range(0, len(bnd_list), 10):
            f.write(", ".join(str(n) for n in bnd_list[i:i+10]) + "\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNBOUND, 1, 3, 0.0\n")

        # Small load on one interior node (or a boundary node if no interior
        # exists) -- exact placement doesn't matter for the Jacobian check itself
        load_node = sorted(interior_nodes)[0] if interior_nodes else patch_node_list[0]
        f.write(f"*CLOAD\n{load_node}, 1, -1.0\n")
        f.write("*NODE FILE\nU\n*END STEP\n")

    print(f"\nWrote: {OUT_INP}")
    print(f"Fixed (boundary) nodes: {len(boundary_nodes)}")
    print(f"Load applied at node: {load_node}")

if __name__ == "__main__":
    main()
