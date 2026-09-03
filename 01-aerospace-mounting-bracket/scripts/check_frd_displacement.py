"""
Project 01 - Aerospace Mounting Bracket
Diagnostic: read nodal displacements directly from the .frd result file
using CalculiX's fixed-width column format, independent of the
*NODE PRINT reference-node output.

CalculiX .frd nodal result line format (Fortran: (1X,I2,I10,6E12.5)):
  columns  1- 3 : "-1" flag (right-justified in 3 chars incl. leading space)
  columns  4-13 : node id (10-char field)
  columns 14-25 : value 1 (12-char E-format)
  columns 26-37 : value 2
  columns 38-49 : value 3
  ...
Fields are NOT reliably whitespace-separated because negative signs can
sit flush against the previous field with no space. A naive .split()
parser silently drops most negative-valued lines - this is the corrected
fixed-width version.
"""

FRD_FILE = "../simulation/01_verification_static.frd"
TIP_NSET_FILE = "../mesh/verification_plate_nsets.inp"

# --- Parse tip node IDs from the nset file ---
tip_nodes = set()
with open(TIP_NSET_FILE) as f:
    lines = f.readlines()

capture = False
for line in lines:
    if line.strip().upper().startswith("*NSET,NSET=TIP"):
        capture = True
        continue
    if capture:
        if line.strip().startswith("*"):
            break
        for tok in line.strip().rstrip(",").split(","):
            tok = tok.strip()
            if tok.isdigit():
                tip_nodes.add(int(tok))

print(f"Parsed {len(tip_nodes)} tip node IDs")

# --- Parse displacement block using fixed-width slicing ---
with open(FRD_FILE) as f:
    frd_lines = f.readlines()

disp = {}
in_disp = False
parse_failures = 0

for line in frd_lines:
    if " DISP" in line and line.rstrip("\n").lstrip().startswith("-4"):
        in_disp = True
        continue
    if in_disp and line.rstrip("\n").lstrip().startswith("-3"):
        break
    if in_disp and line[:3].strip() == "-1":
        try:
            nid = int(line[3:13])
            ux = float(line[13:25])
            uy = float(line[25:37])
            uz = float(line[37:49])
            disp[nid] = (ux, uy, uz)
        except ValueError:
            parse_failures += 1
            continue

print(f"Parsed {len(disp)} nodal displacement records")
print(f"Fixed-width parse failures: {parse_failures}")

if not disp:
    print("NO DISPLACEMENT DATA FOUND IN .frd")
else:
    max_mag = 0.0
    max_node = None
    for nid, (ux, uy, uz) in disp.items():
        mag = (ux**2 + uy**2 + uz**2) ** 0.5
        if mag > max_mag:
            max_mag = mag
            max_node = nid
    print(f"\nMax |U| over entire model: {max_mag:.6f} mm at node {max_node}")

    tip_disps = [disp[n] for n in tip_nodes if n in disp]
    print(f"Tip nodes found in displacement data: {len(tip_disps)} / {len(tip_nodes)}")
    if tip_disps:
        uz_vals = [d[2] for d in tip_disps]
        print(f"Tip Uz range: min={min(uz_vals):.6f}  max={max(uz_vals):.6f}")
        avg_uz = sum(uz_vals) / len(uz_vals)
        print(f"Tip Uz average: {avg_uz:.6f} mm")
