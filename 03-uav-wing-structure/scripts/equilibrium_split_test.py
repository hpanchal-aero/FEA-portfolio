"""
Equilibrium split test (read-only w.r.t. existing results). Builds two
variants of a Stage 2 .inp -- gravity only (CLOAD removed) and lift only
(DLOAD GRAV removed) -- and, after they are solved, compares each root
reaction (sum of RF_Y over the root face) to its own expected value.

Usage:
  python3 scripts/equilibrium_split_test.py make  simulation/stage2_ribs_h2p5
  (solve <job>_gonly and <job>_lonly with run_ccx_watchdog.py)
  python3 scripts/equilibrium_split_test.py check simulation/stage2_ribs_h2p5
"""
import re
import sys
import numpy as np

RHO = 2780e-12
G = 9810.0
ASSEMBLY_VOLUME_MM3 = 117951.9990
TOTAL_LIFT_N = 24.5
TOL = 1e-3


def read_frd_block(path, name, ncomp, n_nodes):
    data = None
    in_block = False
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                in_block = (line[5:13].strip() == name)
                if in_block:
                    data = np.zeros((n_nodes + 1, ncomp))
                continue
            if in_block:
                if line.startswith(" -3"):
                    in_block = False
                elif line.startswith(" -1"):
                    node = int(line[3:13])
                    for k in range(ncomp):
                        data[node, k] = float(line[13 + 12*k:25 + 12*k])
    return data


def read_inp_nodes(path):
    xyz, mode = [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s: continue
            if s.startswith("*"):
                u = s.upper()
                mode = "node" if (u.startswith("*NODE") and not u.startswith("*NODE FILE")) else None
                continue
            if mode == "node":
                p = s.split(",")
                xyz.append([float(v) for v in p[1:4]])
    return np.array(xyz)


def make(job):
    text = open(job + ".inp").read()
    grav_line = "*DLOAD\nEALL, GRAV, 9810.0, 0.0, -1.0, 0.0\n"
    if text.count(grav_line) != 1:
        print("ERROR: gravity block not found exactly once."); sys.exit(1)
    if text.count("*CLOAD\n") != 1:
        print("ERROR: *CLOAD not found exactly once."); sys.exit(1)
    lift_only = text.replace(grav_line, "")
    grav_only = re.sub(r"\*CLOAD\n.*?(?=\*NODE FILE)", "", text, flags=re.S)
    if len(lift_only) >= len(text) or len(grav_only) >= len(text):
        print("ERROR: a variant was not modified."); sys.exit(1)
    open(job + "_lonly.inp", "w").write(lift_only)
    open(job + "_gonly.inp", "w").write(grav_only)
    print(f"Wrote {job}_gonly.inp and {job}_lonly.inp")


def check(job):
    weight = ASSEMBLY_VOLUME_MM3 * RHO * G
    cases = [("gravity only", job + "_gonly", +weight),
             ("lift only", job + "_lonly", -TOTAL_LIFT_N),
             ("combined", job, -(TOTAL_LIFT_N - weight))]
    print(f"Expected weight {weight:.6f} N, applied lift {TOTAL_LIFT_N} N")
    for label, j, expected in cases:
        xyz = read_inp_nodes(j + ".inp")
        n = len(xyz)
        coords = np.vstack([np.zeros((1, 3)), xyz])
        RF = read_frd_block(j + ".frd", "FORC", 3, n)
        z = coords[1:, 2]
        root = np.nonzero(np.abs(z - z.min()) < TOL)[0] + 1
        ry = RF[root, 1].sum()
        print(f"{label:13s}: root RF_Y={ry:+.6f} N, expected {expected:+.6f} N, "
              f"diff={ry - expected:+.6f} N ({100*(ry - expected)/abs(expected):+.4f}%)")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("make", "check"):
        print(__doc__); sys.exit(1)
    (make if sys.argv[1] == "make" else check)(sys.argv[2])
