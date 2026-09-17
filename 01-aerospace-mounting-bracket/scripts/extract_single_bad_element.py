"""
Extract ONE failing element (default: 170300) directly from the exact
.inp file CalculiX rejected, along with its node coordinates, and
build a minimal standalone single-element .inp to test in isolation.
This bypasses all Gmsh-side diagnostics entirely -- we test exactly
the data CalculiX itself consumed.
"""

import sys

# Use the C3D10 curved-mesh solve .inp (the one that actually failed on ccx)
SOURCE_INP = "simulation/stage5_level1_v3_solve.inp"
TARGET_ELEM = 170300
OUT_INP = "simulation/single_element_test.inp"

def main():
    with open(SOURCE_INP) as f:
        lines = f.readlines()

    # --- Find the target element's connectivity line ---
    elem_nodes = None
    mode = None
    for line in lines:
        s = line.strip()
        su = s.upper()
        if su.startswith("*ELEMENT"):
            mode = "elem"
            continue
        if su.startswith("*") and mode == "elem":
            mode = None
        if mode == "elem" and s:
            parts = [p.strip() for p in s.split(",")]
            if int(parts[0]) == TARGET_ELEM:
                elem_nodes = [int(p) for p in parts[1:11]]
                break

    if elem_nodes is None:
        print(f"ERROR: element {TARGET_ELEM} not found in {SOURCE_INP}")
        sys.exit(1)
    print(f"Element {TARGET_ELEM} nodes: {elem_nodes}")

    # --- Find those 10 nodes' coordinates ---
    node_coords = {}
    mode = None
    for line in lines:
        s = line.strip()
        su = s.upper()
        if su == "*NODE":
            mode = "node"
            continue
        if su.startswith("*") and mode == "node":
            mode = None
        if mode == "node" and s:
            parts = [p.strip() for p in s.split(",")]
            tag = int(parts[0])
            if tag in elem_nodes:
                node_coords[tag] = parts[1:4]

    print(f"Found coordinates for {len(node_coords)} / 10 nodes")
    if len(node_coords) != 10:
        print("ERROR: not all 10 nodes found -- check parsing.")
        sys.exit(1)

    for n in elem_nodes:
        print(f"  node {n}: {node_coords[n]}")

    # --- Write minimal standalone .inp ---
    with open(OUT_INP, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Isolated single-element test -- element {TARGET_ELEM} from {SOURCE_INP}\n")
        f.write("*NODE\n")
        for n in elem_nodes:
            c = node_coords[n]
            f.write(f"{n}, {c[0]}, {c[1]}, {c[2]}\n")
        f.write(f"*ELEMENT, TYPE=C3D10, ELSET=ETEST\n")
        f.write(f"{TARGET_ELEM}, " + ", ".join(str(n) for n in elem_nodes) + "\n")
        f.write("*MATERIAL, NAME=AL7075T6\n")
        f.write("*ELASTIC\n71700.0, 0.33\n")
        f.write("*SOLID SECTION, ELSET=ETEST, MATERIAL=AL7075T6\n")
        f.write("*STEP\n*STATIC\n")
        # Fix corner node 1 fully, corner node 2 in y/z, corner node 3 in z
        # (minimal constraint to prevent rigid body motion for a single element)
        f.write(f"*BOUNDARY\n{elem_nodes[0]}, 1, 3, 0.0\n")
        f.write(f"{elem_nodes[1]}, 2, 3, 0.0\n")
        f.write(f"{elem_nodes[2]}, 3, 3, 0.0\n")
        f.write(f"*CLOAD\n{elem_nodes[3]}, 1, -1.0\n")
        f.write("*NODE FILE\nU\n")
        f.write("*END STEP\n")

    print(f"\nWrote standalone test: {OUT_INP}")

if __name__ == "__main__":
    main()
