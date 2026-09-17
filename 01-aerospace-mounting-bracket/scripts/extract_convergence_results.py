"""
Stage 5 -- extract convergence metrics from all 3 Netgen mesh levels:
  1. Global equilibrium check (sanity/verification)
  2. Tip deflection (global stiffness metric, cross-check vs 4a/4b)
  3. Peak von Mises stress near the toe fillet (primary metric of interest)

Uses the same fixed-width .frd parsing validated earlier this session.
"""

import numpy as np

LEVELS = [
    ("Level 1", "simulation/stage5_netgen_L1_solve"),
    ("Level 2", "simulation/stage5_netgen_L2_solve"),
    ("Level 3", "simulation/stage5_netgen_L3_solve"),
]

TOE_POINT = np.array([38.0, 20.0, 2.0])
TOE_RADIUS = 3.0  # mm, region considered "near toe"

def parse_frd_disp_forc(frd_path):
    disp, rf = {}, {}
    with open(frd_path) as f:
        lines = f.readlines()
    mode = None
    for line in lines:
        if line.startswith(" -4") and "DISP" in line:
            mode = "disp"; continue
        if line.startswith(" -4") and "FORC" in line:
            mode = "forc"; continue
        if mode and line.startswith(" -3"):
            mode = None; continue
        if mode == "disp" and line.startswith(" -1"):
            node = int(line[3:13])
            disp[node] = [float(line[13:25]), float(line[25:37]), float(line[37:49])]
        if mode == "forc" and line.startswith(" -1"):
            node = int(line[3:13])
            rf[node] = [float(line[13:25]), float(line[25:37]), float(line[37:49])]
    return disp, rf

def parse_frd_stress(frd_path):
    """Parse nodal stress tensor components (SXX,SYY,SZZ,SXY,SYZ,SZX).

    CONFIRMED FORMAT (verified against raw file dump, 2026-09-16):
    All 6 components are written on a single ' -1' data row, NOT split
    across a '-1' + '-2' continuation pair. Layout per row:
        cols [0:3]   = " -1" row marker
        cols [3:13]  = node ID (10-char field)
        cols [13:25], [25:37], [37:49], [49:61], [61:73], [73:85]
            = SXX, SYY, SZZ, SXY, SYZ, SZX (12-char fields each)
    The block ends at the next non-'-1' line (typically ' -3').
    """
    stress = {}  # node -> [sxx,syy,szz,sxy,syz,szx]
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
            # first non-data line after the block's data rows ends the block
            # (handles '-3' terminator or straight into the next '-4' block)
            mode = False
    return stress

def von_mises(s):
    s11, s22, s33, s12, s13, s23 = s
    return np.sqrt(0.5 * ((s11-s22)**2 + (s22-s33)**2 + (s33-s11)**2
                          + 6*(s12**2 + s13**2 + s23**2)))

def get_node_coords(inp_path):
    coords = {}
    mode = None
    with open(inp_path) as f:
        for line in f:
            s = line.strip()
            su = s.upper()
            if su == "*NODE":
                mode = "node"; continue
            if su.startswith("*"):
                mode = None; continue
            if mode == "node" and s:
                p = s.split(",")
                coords[int(p[0])] = np.array([float(x) for x in p[1:4]])
    return coords

def get_root_nodes(inp_path):
    nodes = []
    capture = False
    with open(inp_path) as f:
        for line in f:
            if line.strip().startswith("*NSET, NSET=NROOT"):
                capture = True; continue
            if capture:
                if line.strip().startswith("*"):
                    break
                nodes.extend(int(x) for x in line.strip().split(",") if x.strip())
    return nodes

print(f"{'Level':<10} {'Equil.err%':<12} {'MaxUx(mm)':<12} {'PeakVM near toe(MPa)':<22} {'#nodes near toe'}")
print("-" * 80)

results = []
for label, base in LEVELS:
    frd = base + ".frd"
    inp = base + ".inp"

    disp, rf = parse_frd_disp_forc(frd)
    coords = get_node_coords(inp)
    root_nodes = get_root_nodes(inp)

    total_rf_x = sum(rf[n][0] for n in root_nodes if n in rf)
    equil_err = abs(total_rf_x - 235.44) / 235.44 * 100

    ux_vals = [v[0] for v in disp.values()]
    max_ux = min(ux_vals)  # most negative

    stress = parse_frd_stress(frd)
    near_toe_vm = []
    for node, s in stress.items():
        if node not in coords:
            continue
        if np.linalg.norm(coords[node] - TOE_POINT) < TOE_RADIUS:
            if len(s) >= 6:
                near_toe_vm.append(von_mises(s[:6]))
    peak_vm = max(near_toe_vm) if near_toe_vm else float("nan")

    print(f"{label:<10} {equil_err:<12.4f} {max_ux:<12.6f} {peak_vm:<22.4f} {len(near_toe_vm)}")
    results.append((label, equil_err, max_ux, peak_vm))

print("\n=== Convergence deltas ===")
for i in range(1, len(results)):
    prev, curr = results[i-1], results[i]
    ux_delta = 100 * abs(curr[2] - prev[2]) / abs(prev[2])
    vm_delta = 100 * abs(curr[3] - prev[3]) / abs(prev[3]) if prev[3] == prev[3] else float("nan")
    print(f"{prev[0]} -> {curr[0]}: Ux change = {ux_delta:.2f}%, Peak VM (near toe) change = {vm_delta:.2f}%")
