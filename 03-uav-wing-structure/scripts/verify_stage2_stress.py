"""
Stress verification (read-only). Compares FEA axial stress sigma_zz on
the OUTER FLANGE FACES with Euler-Bernoulli bending stress
    sigma_zz(z, y) = -M(z) * y / I(z)
(sign: upward lift -> top flange y>0 in compression, bottom flange in
tension) at stations away from the fixed root and from rib junctions.

M(z) comes from the same tapered-beam integration used for the
deflection reference (elliptical lift + spar self-weight; in 'ribs'
mode the ribs enter as point weights, zero bending stiffness).

Selection per station z0: nodes with |z-z0|<5 mm, |x|<3.5 mm (inside
the 8 mm flange, away from the flange edges) and |y| within 0.15 mm of
H(z)/2. Reported: mean FEA stress, mean beam stress, mean and spread
of the node-by-node ratio. Criterion chosen in advance: mean deviation
within 5%.

Usage: python3 scripts/verify_stage2_stress.py <job_path_without_ext> <ribs|spar>
"""
import sys
import numpy as np

SPAN = 1000.0
H_ROOT, H_TIP = 28.8, 14.4
WEB_T, CAP_W, CAP_T = 2.0, 8.0, 3.0
RHO = 2780e-12
G = 9810.0
LIFT = 24.5
RIB_Z = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
RIB_VOL = [11083.3800, 8977.5378, 7093.3632, 5430.8562, 3990.0168, 2770.8450]
RIB_OV = [70.1919, 131.7600, 123.1200, 114.4800, 105.8400, 48.6081]
STATIONS = [100.0, 300.0, 500.0]
HALF_WIN, XLIM, YTOL, CRIT = 5.0, 3.5, 0.15, 5.0


def H(z):
    return H_ROOT + (H_TIP - H_ROOT) * z / SPAN


def I(z):
    h = H(z)
    return (CAP_W * h**3 - (CAP_W - WEB_T) * (h - 2*CAP_T)**3) / 12.0


def A(z):
    h = H(z)
    return (h - 2*CAP_T) * WEB_T + 2 * CAP_W * CAP_T


def moment_profile(rib_weights_on):
    n = 20001
    z = np.linspace(0.0, SPAN, n)
    dz = z[1] - z[0]
    q0 = LIFT / (np.pi * SPAN / 4.0)
    q = q0 * np.sqrt(np.clip(1 - (z / SPAN)**2, 0, None)) - RHO * A(z) * G
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q[i] + q[i+1]) * dz
    if rib_weights_on:
        for zr, v, o in zip(RIB_Z, RIB_VOL, RIB_OV):
            V[z < zr - 1e-9] -= (v - o) * RHO * G
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dz
    return z, M


def read_inp_nodes(path):
    xyz, mode = [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                mode = "node" if (u.startswith("*NODE") and not u.startswith("*NODE FILE")) else None
                continue
            if mode == "node":
                p = s.split(",")
                xyz.append([float(v) for v in p[1:4]])
    return np.array(xyz)


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
    if data is None:
        print(f"ERROR: block {name} not found"); sys.exit(1)
    return data


def main():
    if len(sys.argv) != 3 or sys.argv[2] not in ("ribs", "spar"):
        print(__doc__); sys.exit(1)
    job, mode = sys.argv[1], sys.argv[2]
    xyz = read_inp_nodes(job + ".inp")
    n = len(xyz)
    S = read_frd_block(job + ".frd", "STRESS", 6, n)
    szz = S[1:, 2]
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    zg, M = moment_profile(mode == "ribs")

    print(f"Job: {job} ({mode} mode), {n} nodes; criterion: mean deviation within {CRIT:.0f}%")
    worst = 0.0
    for z0 in STATIONS:
        for side, sgn in (("top   (y>0)", +1.0), ("bottom(y<0)", -1.0)):
            m = ((np.abs(z - z0) < HALF_WIN) & (np.abs(x) < XLIM)
                 & (np.abs(y - sgn * H(z) / 2.0) < YTOL))
            if m.sum() == 0:
                print(f"  z={z0:5.0f} {side}: no nodes selected"); continue
            Mn = np.interp(z[m], zg, M)
            beam = -Mn * y[m] / I(z[m])
            fea = szz[m]
            ratio = fea / beam
            dev = 100 * (ratio.mean() - 1.0)
            worst = max(worst, abs(dev))
            print(f"  z={z0:5.0f} {side}: {m.sum():4d} nodes  FEA {fea.mean():+8.3f}  "
                  f"beam {beam.mean():+8.3f} MPa  ratio {ratio.mean():.4f} "
                  f"(std {ratio.std():.4f})  deviation {dev:+.2f}%")
    print(f"Largest mean deviation: {worst:.2f}%  -> "
          f"{'within' if worst <= CRIT else 'OUTSIDE'} the {CRIT:.0f}% criterion")


if __name__ == "__main__":
    main()
