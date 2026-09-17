"""
Diagnose the ~462 nonpositive-Jacobian elements (tags 170259-170720)
reported by ccx. Read their corner-node coordinates directly from the
solver .inp (the exact file that was fed to ccx) and check proximity
to each of the three refined features, to test the hypothesis that
these are curvilinear order-2 elements with midside nodes over-pulled
onto the tight R1mm toe-fillet curvature.
"""

import re
import numpy as np

INP_FILE = "simulation/stage5_level1_v3_solve.inp"

BAD_ELEM_START = 170259
BAD_ELEM_END = 170720

TOE_APPROX = np.array([38.0, 20.0, 2.0])
R3_APPROX = np.array([56.5, 20.0, 3.5])
HOLE_APPROX = np.array([20.0, 20.0, 0.0])

def main():
    nodes = {}
    elements = {}

    with open(INP_FILE) as f:
        lines = f.readlines()

    mode = None
    for line in lines:
        s = line.strip()
        if s.upper() == "*NODE" or s.upper().startswith("*NODE, NSET"):
            mode = "node" if "NSET" not in s.upper() else None
            continue
        if s.upper().startswith("*ELEMENT"):
            mode = "elem"
            continue
        if s.startswith("*"):
            mode = None
            continue
        if mode == "node" and s:
            parts = s.split(",")
            tag = int(parts[0])
            coords = [float(p) for p in parts[1:4]]
            nodes[tag] = np.array(coords)
        elif mode == "elem" and s:
            parts = [p.strip() for p in s.split(",")]
            tag = int(parts[0])
            if BAD_ELEM_START <= tag <= BAD_ELEM_END:
                corner_nodes = [int(p) for p in parts[1:5]]  # first 4 = corners
                elements[tag] = corner_nodes

    print(f"Parsed {len(nodes)} nodes, {len(elements)} target elements "
          f"(expected {BAD_ELEM_END - BAD_ELEM_START + 1})")

    if len(elements) == 0:
        print("ERROR: no matching elements found -- check tag range/parsing.")
        return

    dists_toe, dists_r3, dists_hole = [], [], []
    for tag, corner_ids in elements.items():
        coords = np.array([nodes[n] for n in corner_ids])
        centroid = coords.mean(axis=0)
        dists_toe.append(np.linalg.norm(centroid - TOE_APPROX))
        dists_r3.append(np.linalg.norm(centroid - R3_APPROX))
        dists_hole.append(np.linalg.norm(centroid - HOLE_APPROX))

    dists_toe = np.array(dists_toe)
    dists_r3 = np.array(dists_r3)
    dists_hole = np.array(dists_hole)

    print(f"\nDistance from bad-element centroids to TOE fillet: "
          f"min={dists_toe.min():.3f}, mean={dists_toe.mean():.3f}, max={dists_toe.max():.3f}")
    print(f"Distance from bad-element centroids to R3 fillet:  "
          f"min={dists_r3.min():.3f}, mean={dists_r3.mean():.3f}, max={dists_r3.max():.3f}")
    print(f"Distance from bad-element centroids to hole:       "
          f"min={dists_hole.min():.3f}, mean={dists_hole.mean():.3f}, max={dists_hole.max():.3f}")

    nearest_toe = (dists_toe < dists_r3) & (dists_toe < dists_hole)
    nearest_r3 = (dists_r3 < dists_toe) & (dists_r3 < dists_hole)
    nearest_hole = (dists_hole < dists_toe) & (dists_hole < dists_r3)
    print(f"\nElements nearest to TOE fillet:  {nearest_toe.sum()} / {len(elements)} "
          f"({100*nearest_toe.sum()/len(elements):.1f}%)")
    print(f"Elements nearest to R3 fillet:   {nearest_r3.sum()} / {len(elements)} "
          f"({100*nearest_r3.sum()/len(elements):.1f}%)")
    print(f"Elements nearest to hole:        {nearest_hole.sum()} / {len(elements)} "
          f"({100*nearest_hole.sum()/len(elements):.1f}%)")

if __name__ == "__main__":
    main()
