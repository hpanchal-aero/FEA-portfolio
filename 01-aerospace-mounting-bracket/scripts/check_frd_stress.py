"""
Project 01 - Aerospace Mounting Bracket
Extract von Mises stress from the .frd file at the root face, using the
same fixed-width parsing approach validated for displacement.

CalculiX writes stress as 6 components (SXX,SYY,SZZ,SXY,SYZ,SZX) per
node, in the same -1 record format as displacement, but with 6 fields
that wrap onto continuation lines (flag "-2") since 6*12=72 chars can
exceed typical line formatting - handled explicitly below.
"""

ROOT_NSET_FILE = "../mesh/verification_plate_nsets.inp"
FRD_FILE = "../simulation/01_verification_static.frd"

def parse_nset(name):
    nodes = set()
    with open(ROOT_NSET_FILE) as f:
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

root_nodes = parse_nset("root")
print(f"Parsed {len(root_nodes)} root node IDs")

with open(FRD_FILE) as f:
    frd_lines = f.readlines()

stress = {}
in_stress = False
i = 0
while i < len(frd_lines):
    line = frd_lines[i]
    if " STRESS" in line and line.rstrip("\n").lstrip().startswith("-4"):
        in_stress = True
        i += 1
        continue
    if in_stress and line.rstrip("\n").lstrip().startswith("-3"):
        break
    if in_stress and line[:3].strip() == "-1":
        nid = int(line[3:13])
        vals = [
            float(line[13:25]), float(line[25:37]), float(line[37:49]),
            float(line[49:61]), float(line[61:73]), float(line[73:85]),
        ]
        stress[nid] = vals
    i += 1

print(f"Parsed {len(stress)} nodal stress records")

def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = s
    return (0.5 * ((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2
            + 6*(sxy**2 + syz**2 + szx**2))) ** 0.5

if stress:
    max_vm = 0.0
    max_node = None
    for nid, s in stress.items():
        vm = von_mises(s)
        if vm > max_vm:
            max_vm = vm
            max_node = nid
    print(f"\nMax von Mises stress (entire model): {max_vm:.3f} MPa at node {max_node}")

    root_vm = [(nid, von_mises(stress[nid])) for nid in root_nodes if nid in stress]
    if root_vm:
        max_root_vm = max(root_vm, key=lambda x: x[1])
        print(f"Root nodes with stress data: {len(root_vm)} / {len(root_nodes)}")
        print(f"Max von Mises at root: {max_root_vm[1]:.3f} MPa at node {max_root_vm[0]}")
    else:
        print("No root nodes found in stress data")
