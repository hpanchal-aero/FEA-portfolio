"""
Project 01 - Aerospace Mounting Bracket
Extract tip displacement and root-region von Mises stress from the
structured hex mesh .frd file, using the same fixed-width parsing
approach validated for the tet mesh (check_frd_displacement.py,
check_frd_stress.py), extended for C3D20's 20-node connectivity.
"""

import os

FRD_FILE = "../mesh/hex_test/analysis_hex.frd"
NSET_FILE = "../mesh/hex_test/nsets_hex.inp"

ROOT_X_MIN = 3.0
ROOT_X_MAX = 9.0

def parse_nset(name):
    nodes = set()
    with open(NSET_FILE) as f:
        lines = f.readlines()
    capture = False
    for line in lines:
        if line.strip().upper() == f"*NSET,NSET={name.upper()}":
            capture = True
            continue
        if capture:
            if line.strip().startswith("*"):
                break
            for tok in line.strip().rstrip(",").split(","):
                tok = tok.strip()
                if tok.isdigit():
                    nodes.add(int(tok))
    return nodes

tip_nodes = parse_nset("tip")
print(f"Parsed {len(tip_nodes)} tip node IDs")

with open(FRD_FILE) as f:
    frd_lines = f.readlines()

def parse_frd_block(marker, ncols):
    data = {}
    in_block = False
    for line in frd_lines:
        if marker in line and line.rstrip("\n").lstrip().startswith("-4"):
            in_block = True
            continue
        if in_block and line.rstrip("\n").lstrip().startswith("-3"):
            break
        if in_block and line[:3].strip() == "-1":
            nid = int(line[3:13])
            vals = []
            for c in range(ncols):
                start = 13 + c * 12
                end = start + 12
                vals.append(float(line[start:end]))
            data[nid] = vals
    return data

def parse_node_coords():
    coords = {}
    in_block = False
    for line in frd_lines:
        if line.strip().startswith("2C"):
            in_block = True
            continue
        if in_block and line[:3].strip() == "-3":
            break
        if in_block and line[:3].strip() == "-1":
            nid = int(line[3:13])
            x = float(line[13:25])
            y = float(line[25:37])
            z = float(line[37:49])
            coords[nid] = (x, y, z)
    return coords

disp = parse_frd_block("DISP", 3)
stress = parse_frd_block("STRESS", 6)
coords = parse_node_coords()

print(f"Parsed {len(disp)} displacement records")
print(f"Parsed {len(stress)} stress records")
print(f"Parsed {len(coords)} node coordinate records")

tip_uz = [disp[n][2] for n in tip_nodes if n in disp]
if tip_uz:
    avg_uz = sum(tip_uz) / len(tip_uz)
    print(f"\nTip Uz average: {avg_uz:.6f} mm")
    print(f"Tip Uz range: min={min(tip_uz):.6f} max={max(tip_uz):.6f}")
    print(f"Analytical target (Euler-Bernoulli): -1.108000 mm")
    pct = abs(avg_uz - (-1.108)) / 1.108 * 100
    print(f"Percent difference: {pct:.2f}%")
else:
    print("NO TIP DISPLACEMENT DATA FOUND")

def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = s
    return (0.5 * ((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2
            + 6*(sxy**2 + syz**2 + szx**2))) ** 0.5

root_band_nodes = [nid for nid, (x,y,z) in coords.items()
                   if ROOT_X_MIN <= x <= ROOT_X_MAX]
root_vm = [von_mises(stress[n]) for n in root_band_nodes if n in stress]
if root_vm:
    max_vm = max(root_vm)
    print(f"\nMax von Mises in root band [{ROOT_X_MIN},{ROOT_X_MAX}] mm: "
          f"{max_vm:.3f} MPa")
    print(f"Analytical target: 132.4 MPa")
    pct = abs(max_vm - 132.4) / 132.4 * 100
    print(f"Percent difference: {pct:.2f}%")
else:
    print("NO ROOT-BAND STRESS DATA FOUND")
