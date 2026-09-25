"""
Config A -- load-application test for the Configuration B study.

Converts a Config A CalculiX input that uses uniform top-face pressure
(*DSLOAD) into one that uses a common top-face axial displacement: every node
in NTOP has Uz tied to a reference node by *EQUATION, and the total force is
applied as a concentrated force on that reference node. The load then splits
between members by stiffness (equal-strain sharing), instead of assuming equal
stress.

Usage: /usr/bin/python3 scripts/make_equation_load_case.py <src.inp> <dst.inp>

The reference node is placed at (0, 0, 150), off the structure, and its Ux, Uy
are fixed so the system is not singular. Boundary conditions on the base are
left exactly as in the source file.
"""
import sys

FORCE = 156.5779                 # N = 0.1643 MPa x 953 mm^2 (same total as pressure case)
REF_XYZ = (0.0, 0.0, 150.0)


def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        print(f"ERROR: expected exactly 1 occurrence of {label}, found {n}")
        sys.exit(1)
    return text.replace(old, new)


def main():
    if len(sys.argv) != 3:
        print("Usage: make_equation_load_case.py <src.inp> <dst.inp>")
        sys.exit(1)
    src, dst = sys.argv[1], sys.argv[2]
    text = open(src).read()

    node_ids, ntop, mode = [], [], None
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("*"):
            u = s.upper().replace(" ", "")
            if u.startswith("*NODE") and not u.startswith("*NODEFILE"):
                mode = "node"
            elif u.startswith("*NSET") and "NSET=NTOP" in u:
                mode = "ntop"
            else:
                mode = None
            continue
        if mode == "node":
            node_ids.append(int(s.split(",")[0]))
        elif mode == "ntop":
            ntop += [int(v) for v in s.split(",") if v.strip()]

    n_nodes = len(node_ids)
    if node_ids != list(range(1, n_nodes + 1)):
        print("ERROR: node ids are not contiguous 1..N.")
        sys.exit(1)
    if not ntop or len(set(ntop)) != len(ntop):
        print("ERROR: NTOP is empty or has duplicates.")
        sys.exit(1)
    ref = n_nodes + 1
    print(f"Source: {src}")
    print(f"Nodes: {n_nodes}; NTOP nodes: {len(ntop)}; reference node id: {ref}")

    # 1. replace the pressure block with a concentrated force on the reference node
    i = text.find("*DSLOAD\n")
    j = text.find("*NODE FILE\n")
    if text.count("*DSLOAD\n") != 1 or text.count("*NODE FILE\n") != 1 or not (0 <= i < j):
        print("ERROR: could not isolate the *DSLOAD block.")
        sys.exit(1)
    text = text[:i] + f"*CLOAD\n{ref}, 3, {-FORCE}\n" + text[j:]

    # 2. fix the reference node laterally (first lines of the *BOUNDARY block)
    text = replace_once(text, "*BOUNDARY\n", f"*BOUNDARY\n{ref}, 1, 2, 0.0\n", "*BOUNDARY header")

    # 3. tie Uz of every top node to the reference node (model data, before *STEP)
    eq = "*EQUATION\n" + "".join(f"2\n{k}, 3, 1., {ref}, 3, -1.\n" for k in ntop)
    text = replace_once(text, "*STEP\n", eq + "*STEP\n", "*STEP header")

    # 4. add the reference node to the node list
    text = replace_once(
        text, "*ELEMENT, TYPE=C3D10, ELSET=EALL\n",
        f"{ref}, {REF_XYZ[0]:.10f}, {REF_XYZ[1]:.10f}, {REF_XYZ[2]:.10f}\n"
        "*ELEMENT, TYPE=C3D10, ELSET=EALL\n", "*ELEMENT header")

    open(dst, "w").write(text)
    print(f"Wrote {dst}: {len(ntop)} equations, force {-FORCE} N on node {ref} (Z)")


if __name__ == "__main__":
    main()
