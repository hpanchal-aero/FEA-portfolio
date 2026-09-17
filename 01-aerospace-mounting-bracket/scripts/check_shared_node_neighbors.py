"""
For failing element 170300, find every OTHER element in the full mesh
that shares at least one node with it, and check each neighbor's own
linear-corner-volume sign. If a neighbor has a very different volume
scale or a topology suggesting inconsistent orientation, that's a
concrete candidate for the context-dependent failure -- since element
170300 was proven valid in total isolation.
"""

import numpy as np

INP_FILE = "simulation/stage5_level1_v3_solve.inp"
TARGET_ELEM = 170300

def signed_volume(p0, p1, p2, p3):
    return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

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
            nodes[int(parts[0])] = np.array([float(p) for p in parts[1:4]])
        elif mode == "elem" and s:
            parts = [p.strip() for p in s.split(",")]
            elements[int(parts[0])] = [int(p) for p in parts[1:11]]

    print(f"Parsed {len(nodes)} nodes, {len(elements)} elements")

    target_nodes = set(elements[TARGET_ELEM])
    print(f"\nTarget element {TARGET_ELEM} nodes: {sorted(target_nodes)}")

    neighbors = []
    for etag, enodes in elements.items():
        if etag == TARGET_ELEM:
            continue
        shared = target_nodes & set(enodes)
        if shared:
            corners = np.array([nodes[n] for n in enodes[:4]])
            vol = signed_volume(corners[0], corners[1], corners[2], corners[3])
            neighbors.append((etag, len(shared), vol))

    print(f"\nFound {len(neighbors)} elements sharing at least 1 node with {TARGET_ELEM}")

    neighbors.sort(key=lambda x: -x[1])  # most shared nodes first
    print(f"\nTop 15 neighbors by shared node count:")
    for etag, nshared, vol in neighbors[:15]:
        flag = "NEGATIVE VOLUME <-- SUSPECT" if vol < 0 else ""
        print(f"  elem {etag}: shares {nshared} nodes, linear volume = {vol:.6f}  {flag}")

    neg_vol_neighbors = [n for n in neighbors if n[2] < 0]
    print(f"\nTotal neighbors with NEGATIVE linear volume: {len(neg_vol_neighbors)}")
    if neg_vol_neighbors:
        print("These would be genuinely inverted/degenerate elements adjacent")
        print("to our test element -- worth checking if any of them are also")
        print("in the original 462-element failure list.")

if __name__ == "__main__":
    main()
