"""
Read-only check of the systematic root-reaction residual.

Hypothesis: the reported root RF (sum of FORC over the fixed nodes)
excludes external load applied DIRECTLY on those fixed nodes. If so:
    reported RF_Y = R_true + F_ext_on_fixed_nodes,   R_true = -F_ext_total
For a straight-sided C3D10 element with uniform body force, the
consistent nodal share is -V/20 per corner node and +V/5 per midside
node (sum = V). CLOAD entries on fixed nodes are added directly.

Usage: python3 scripts/check_equilibrium_constrained_loads.py
"""
import numpy as np

RHO = 2780e-12
G = 9810.0
JOBS = ["simulation/stage1_spar_h2p0", "simulation/stage1_spar_h1p5",
        "simulation/stage2_ribs_h2p0", "simulation/stage2_ribs_h2p5"]
COEF = np.array([-1/20.0] * 4 + [1/5.0] * 6)


def parse_inp(path):
    nodes, elems, nroot, cload = [], [], [], []
    mode = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                if u.startswith("*NODE") and not u.startswith("*NODE FILE"):
                    mode = "node"
                elif u.startswith("*ELEMENT"):
                    mode = "elem"
                elif u.startswith("*NSET") and "NROOT" in u:
                    mode = "nroot"
                elif u.startswith("*CLOAD"):
                    mode = "cload"
                else:
                    mode = None
                continue
            p = [t.strip() for t in s.split(",") if t.strip()]
            if mode == "node":
                nodes.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                elems.append([int(v) for v in p[1:11]])
            elif mode == "nroot":
                nroot.extend(int(v) for v in p)
            elif mode == "cload":
                cload.append((int(p[0]), int(p[1]), float(p[2])))
    return np.array(nodes), np.array(elems), np.array(nroot), cload


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


for job in JOBS:
    xyz, elems, nroot, cload = parse_inp(job + ".inp")
    n = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    P = coords[elems[:, :4]]
    vol = np.abs(np.einsum("ij,ij->i",
                           np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                           P[:, 3] - P[:, 0])) / 6.0
    is_root = np.zeros(n + 1, dtype=bool)
    is_root[nroot] = True

    weight = RHO * G * vol.sum()
    grav_share_vol = (vol * (is_root[elems] * COEF).sum(axis=1)).sum()
    F_grav_root_Y = -RHO * G * grav_share_vol
    F_lift_total = sum(f for nd, dof, f in cload if dof == 2)
    F_lift_root = sum(f for nd, dof, f in cload if dof == 2 and is_root[nd])
    F_ext_total = F_lift_total - weight
    R_true = -F_ext_total
    expected_naive = R_true
    expected_corr = R_true + F_grav_root_Y + F_lift_root

    RF = read_frd_block(job + ".frd", "FORC", 3, n)
    measured = RF[nroot, 1].sum()

    print(f"\n{job}")
    print(f"  weight {weight:.6f} N, lift {F_lift_total:.6f} N, "
          f"root nodes {len(nroot)}")
    print(f"  external load on fixed nodes: gravity share {F_grav_root_Y:+.6f} N, "
          f"CLOAD {F_lift_root:+.6f} N")
    print(f"  measured root RF_Y        : {measured:+.6f} N")
    print(f"  expected (naive)          : {expected_naive:+.6f} N  "
          f"(diff {measured - expected_naive:+.6f} N, "
          f"{100*(measured - expected_naive)/abs(expected_naive):+.4f}%)")
    print(f"  expected (fixed-node load excluded): {expected_corr:+.6f} N  "
          f"(diff {measured - expected_corr:+.6f} N, "
          f"{100*(measured - expected_corr)/abs(expected_corr):+.4f}%)")
