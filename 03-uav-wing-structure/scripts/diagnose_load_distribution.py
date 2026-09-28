"""
Read-only diagnostic: what load did the .inp files actually apply?
Compares the *CLOAD resultant, centroid, and spanwise distribution
(100mm bins) against the intended elliptical lift
q(z) = q0*sqrt(1-(z/b)^2), total 24.5 N over b=1000mm.
Stage 1 files are included as a control (Stage 1 converged).

Usage: python3 scripts/diagnose_load_distribution.py
"""
import numpy as np

SPAN = 1000.0
TOTAL = 24.5
JOBS = ["simulation/stage1_spar_h2p0", "simulation/stage1_spar_h1p5",
        "simulation/stage2_ribs_h2p0", "simulation/stage2_ribs_h2p5"]
EDGES = np.linspace(0.0, SPAN, 11)


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


def read_cload(path):
    nodes, forces, mode = [], [], False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s: continue
            if s.startswith("*"):
                mode = s.upper().startswith("*CLOAD")
                continue
            if mode:
                p = s.split(",")
                nodes.append(int(p[0]))
                forces.append(float(p[2]))
    return np.array(nodes), np.array(forces)


def cum_intended(z):
    u = np.clip(z / SPAN, 0.0, 1.0)
    q0 = TOTAL / (np.pi * SPAN / 4.0)
    return q0 * (SPAN / 2.0) * (u * np.sqrt(1 - u**2) + np.arcsin(u))


intended_bins = np.diff(cum_intended(EDGES))
M_intended = TOTAL * 4.0 * SPAN / (3.0 * np.pi)
print(f"Intended: total {TOTAL} N, centroid z={4*SPAN/(3*np.pi):.2f} mm, "
      f"moment about z=0: {M_intended:.1f} N*mm")
print("Intended % of load per 100mm bin: " +
      " ".join(f"{100*b/TOTAL:5.1f}" for b in intended_bins))

for job in JOBS:
    xyz = read_inp_nodes(job + ".inp")
    nodes, forces = read_cload(job + ".inp")
    z = xyz[nodes - 1, 2]
    F = forces.sum()
    M = (forces * z).sum()
    binned, _ = np.histogram(z, bins=EDGES, weights=forces)
    print(f"\n{job}")
    print(f"  load nodes: {len(nodes)}, total F={F:.4f} N, centroid z={M/F:.2f} mm, "
          f"moment={M:.1f} N*mm ({100*(M - M_intended)/M_intended:+.2f}% vs intended)")
    print("  % of load per 100mm bin:         " +
          " ".join(f"{100*b/F:5.1f}" for b in binned))
    below0 = (z < 0.0).sum()
    print(f"  load nodes at z<0 (inside root rib overhang): {below0}")
