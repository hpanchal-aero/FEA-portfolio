"""
Bisection step: build a standalone .inp containing element 170300
PLUS all elements that share at least one node with it (its full
1-ring neighbor patch, 82 elements found previously). If this patch
reproduces the Jacobian failure, we've found a minimal failing
context far smaller than the full 157k-element model, and can
inspect it directly. If it does NOT fail, the trigger requires an
even larger neighborhood, and we grow the patch further (2-ring).
"""

import numpy as np

INP_FILE = "simulation/stage5_level1_v3_solve.inp"
OUT_INP = "simulation/neighbor_patch_test.inp"
TARGET_ELEM = 170300

def main():
    nodes = {}
    elements = {}

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
            elements[int(parts[0])] = [int(p) for p in parts[1:11]]

    target_nodes = set(elements[TARGET_ELEM])
    patch_elems = {TARGET_ELEM: elements[TARGET_ELEM]}
    for etag, enodes in elements.items():
        if etag == TARGET_ELEM:
            continue
        if target_nodes & set(enodes):
            patch_elems[etag] = enodes

    print(f"Patch size (1-ring neighbors + target): {len(patch_elems)} elements")

    patch_nodes = set()
    for enodes in patch_elems.values():
        patch_nodes.update(enodes)
    print(f"Patch node count: {len(patch_nodes)}")

    # Find the patch's outer boundary nodes (any node NOT shared by
    # every element touching it is fine to leave free; for a minimal
    # rigid-body-safe test, just fix a handful of well-spread nodes)
    patch_node_list = sorted(patch_nodes)

    with open(OUT_INP, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Neighbor patch test -- element {TARGET_ELEM} + {len(patch_elems)-1} neighbors\n")
        f.write("*NODE\n")
        for n in patch_node_list:
            c = nodes[n]
            f.write(f"{n}, {c[0]}, {c[1]}, {c[2]}\n")
        f.write("*ELEMENT, TYPE=C3D10, ELSET=EPATCH\n")
        for etag, enodes in patch_elems.items():
            f.write(f"{etag}, " + ", ".join(str(n) for n in enodes) + "\n")
        f.write("*MATERIAL, NAME=AL7075T6\n*ELASTIC\n71700.0, 0.33\n")
        f.write("*SOLID SECTION, ELSET=EPATCH, MATERIAL=AL7075T6\n")
        f.write("*STEP\n*STATIC\n")
        # Minimal rigid-body constraint: fix first 3 non-collinear nodes fully
        fixed = patch_node_list[:4]
        f.write(f"*BOUNDARY\n{fixed[0]}, 1, 3, 0.0\n{fixed[1]}, 1, 3, 0.0\n"
                f"{fixed[2]}, 1, 3, 0.0\n")
        f.write(f"*CLOAD\n{patch_node_list[-1]}, 1, -1.0\n")
        f.write("*NODE FILE\nU\n*END STEP\n")

    print(f"\nWrote: {OUT_INP}")

if __name__ == "__main__":
    main()
