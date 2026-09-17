"""
Stage 5 -- locate and characterize the global peak von Mises stress node
in the L3 solve, to determine whether it's a real physical feature
(hole edge, root boundary layer, fillet) or a mesh-transition artifact.
"""

import numpy as np

INP_FILE = "simulation/stage5_netgen_L3_solve.inp"
FRD_FILE = "simulation/stage5_netgen_L3_solve.frd"

TOE_POINT = np.array([38.0, 20.0, 2.0])
HOLE_CENTER = np.array([20.0, 20.0, 4.0])  # flange A hole, approx z at surface
R3_REGION = np.array([55.0, 20.0, 2.0])    # approx R3 fillet junction area

def parse_inp_nodes_elements(inp_path):
    nodes = {}
    elements = []
    mode = None
    with open(inp_path) as f:
        for line in f:
            s = line.strip()
            su = s.upper()
            if su == "*NODE":
                mode = "node"; continue
            if su.startswith("*ELEMENT") and "C3D10" in su:
                mode = "elem"; continue
            if su.startswith("*"):
                mode = None; continue
            if not s:
                continue
            parts = [p.strip() for p in s.split(",")]
            if mode == "node":
                tag = int(parts[0])
                nodes[tag] = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            elif mode == "elem":
                eid = int(parts[0])
                conn = [int(x) for x in parts[1:]]
                elements.append((eid, conn))
    return nodes, elements

def parse_frd_stress(frd_path):
    stress = {}
    with open(frd_path) as f:
        lines = f.readlines()
    mode = False
    for line in lines:
        if line.startswith(" -4") and "STRESS" in line:
            mode = True; continue
        if mode and line.startswith(" -1"):
            node = int(line[3:13])
            stress[node] = [
                float(line[13:25]), float(line[25:37]), float(line[37:49]),
                float(line[49:61]), float(line[61:73]), float(line[73:85]),
            ]
            continue
        if mode and not line.startswith(" -1") and not line.startswith(" -5"):
            mode = False
    return stress

def von_mises(s):
    s11, s22, s33, s12, s13, s23 = s
    return np.sqrt(0.5 * ((s11-s22)**2 + (s22-s33)**2 + (s33-s11)**2
                          + 6*(s12**2 + s13**2 + s23**2)))

def main():
    nodes, elements = parse_inp_nodes_elements(INP_FILE)
    stress = parse_frd_stress(FRD_FILE)

    vm = {n: von_mises(s) for n, s in stress.items()}
    ranked = sorted(vm.items(), key=lambda kv: kv[1], reverse=True)

    print(f"{'Rank':<5}{'Node':<10}{'VM(MPa)':<12}{'x':<10}{'y':<10}{'z':<10}"
          f"{'dist_toe':<10}{'dist_hole':<10}{'dist_x0(root)':<14}")
    print("-" * 90)
    for i, (n, v) in enumerate(ranked[:20], start=1):
        if n not in nodes:
            print(f"{i:<5}{n:<10}{v:<12.2f}  (node not in .inp coordinate table)")
            continue
        c = nodes[n]
        d_toe = np.linalg.norm(c - TOE_POINT)
        d_hole = np.linalg.norm(c - HOLE_CENTER)
        d_root = abs(c[0] - 0.0)
        print(f"{i:<5}{n:<10}{v:<12.2f}{c[0]:<10.3f}{c[1]:<10.3f}{c[2]:<10.3f}"
              f"{d_toe:<10.3f}{d_hole:<10.3f}{d_root:<14.3f}")

    # Elements touching the single top-ranked node, for aspect-ratio inspection
    top_node = ranked[0][0]
    print(f"\n--- Elements containing top node {top_node} (VM={ranked[0][1]:.2f} MPa) ---")
    touching = [(eid, conn) for eid, conn in elements if top_node in conn]
    print(f"{len(touching)} elements touch this node")
    for eid, conn in touching[:5]:
        corner_nodes = conn[:4]  # C3D10 corners are first 4 entries
        coords = [nodes[c] for c in corner_nodes if c in nodes]
        if len(coords) == 4:
            edge_lens = []
            for i in range(4):
                for j in range(i+1, 4):
                    edge_lens.append(np.linalg.norm(coords[i] - coords[j]))
            print(f"  elem {eid}: corner edge lengths min={min(edge_lens):.4f} "
                  f"max={max(edge_lens):.4f} ratio={max(edge_lens)/min(edge_lens):.2f}")

if __name__ == "__main__":
    main()
