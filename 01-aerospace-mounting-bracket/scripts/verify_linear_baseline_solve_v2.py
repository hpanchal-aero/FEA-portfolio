"""
Stage 5 -- verify linear baseline solve, v2: fixed-width column parsing
instead of whitespace-split, to correctly handle negative values that
run together with the preceding field in CalculiX's Fortran fixed-format
.frd output (e.g. "1-2.46119E-02" with no separating space).
"""

FRD_FILE = "simulation/stage5_level1_LINEAR_TEMP_solve.frd"

def parse_frd_block(frd_path, block_keyword):
    """Parse a -4 block (DISP or FORC) using fixed-width column slicing."""
    records = {}
    with open(frd_path) as f:
        lines = f.readlines()

    in_block = False
    for line in lines:
        if line.startswith(" -4") and block_keyword in line:
            in_block = True
            continue
        if in_block and line.startswith(" -3"):
            break
        if in_block and line.startswith(" -1"):
            node = int(line[3:13])
            v1 = float(line[13:25])
            v2 = float(line[25:37])
            v3 = float(line[37:49])
            records[node] = [v1, v2, v3]

    return records

disp = parse_frd_block(FRD_FILE, "DISP")
rf = parse_frd_block(FRD_FILE, "FORC")

print(f"Parsed {len(disp)} displacement records, {len(rf)} reaction force records")

# --- Equilibrium check ---
total_rf_x = sum(v[0] for v in rf.values())
print(f"\nSum of Fx reactions (all nodes with RF recorded): {total_rf_x:.4f} N")
print(f"Applied load: -235.44 N; expected reaction sum: +235.44 N")
if abs(total_rf_x) > 1e-6:
    err = abs(total_rf_x - 235.44) / 235.44 * 100
    print(f"Equilibrium error: {err:.4f}%")

# --- Reference/tip node displacement ---
REF_NODE = 28563
if REF_NODE in disp:
    print(f"\nReference node {REF_NODE} displacement (Ux,Uy,Uz): {disp[REF_NODE]}")
else:
    print(f"\nReference node {REF_NODE} not found in DISP block.")

# --- Global max |Ux| ---
if disp:
    ux_vals = [v[0] for v in disp.values()]
    print(f"\nMax |Ux| in model: {max(abs(min(ux_vals)), abs(max(ux_vals))):.6f} mm")
    print(f"Ux range: [{min(ux_vals):.6f}, {max(ux_vals):.6f}] mm")
