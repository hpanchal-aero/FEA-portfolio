"""
Stage 2 diagnostic (read-only, no changes to any solver input/output).

1) Tip-band Uy vs x linear regression per mesh level: the slope shows
   twist-like variation across the chord; the intercept is the
   spar-axis (x=0) deflection, which is what beam theory refers to.
2) Beam-theory reference for spar + ribs where ribs are treated as
   point WEIGHTS with zero bending stiffness (thin plates normal to the
   span). Isolates the effect of added self-weight from any stiffness
   effect.

Usage: python3 scripts/diagnose_stage2_tip.py
"""
import numpy as np

SPAN = 1000.0
H_ROOT, H_TIP = 28.8, 14.4
WEB_T, CAP_W, CAP_T = 2.0, 8.0, 3.0
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
TOTAL_LIFT_N = 24.5

RIB_Z = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
RIB_VOL = [11083.3800, 8977.5378, 7093.3632, 5430.8562, 3990.0168, 2770.8450]
RIB_OVERLAP = [70.1919, 131.7600, 123.1200, 114.4800, 105.8400, 48.6081]


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


def I_of_z(z):
    H = H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    return (1/12.0) * (CAP_W * H**3 - (CAP_W - WEB_T) * (H - 2*CAP_T)**3)


def A_of_z(z):
    H = H_ROOT + (H_TIP - H_ROOT) * (z / SPAN)
    return (H - 2*CAP_T) * WEB_T + 2 * CAP_W * CAP_T


def beam_tip(rib_weights_on):
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = TOTAL_LIFT_N / (np.pi * SPAN / 4.0)
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN)**2, 0, None)) - RHO * A_of_z(z) * G

    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q[i] + q[i+1]) * dz
    if rib_weights_on:
        for zr, vol, ov in zip(RIB_Z, RIB_VOL, RIB_OVERLAP):
            P = (vol - ov) * RHO * G          # N, downward
            V[z < zr - 1e-9] -= P             # shear jump inboard of the rib
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dz
    curv = M / (E_MOD * I_of_z(z))
    th = np.zeros(n); v = np.zeros(n)
    for i in range(1, n):
        th[i] = th[i-1] + 0.5 * (curv[i-1] + curv[i]) * dz
    for i in range(1, n):
        v[i] = v[i-1] + 0.5 * (th[i-1] + th[i]) * dz
    return v[-1]


print("=" * 70)
print("BEAM-THEORY REFERENCE (ribs as point weights, zero bending stiffness)")
print("=" * 70)
v_spar = beam_tip(False)
v_ribw = beam_tip(True)
print(f"Spar only:               {v_spar:.4f} mm")
print(f"Spar + rib weights only: {v_ribw:.4f} mm  "
      f"({100*(v_ribw - v_spar)/v_spar:+.2f}% from added self-weight alone)")

for label, job in [("2.0mm", "simulation/stage2_ribs_h2p0"),
                   ("2.5mm", "simulation/stage2_ribs_h2p5")]:
    xyz = read_inp_nodes(job + ".inp")
    n = len(xyz)
    U = read_frd_block(job + ".frd", "DISP", 3, n)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    band = np.nonzero(np.abs(coords[1:, 2] - SPAN) < 2.0)[0] + 1
    x, uy = coords[band, 0], U[band, 1]
    slope, intercept = np.polyfit(x, uy, 1)
    axis = band[np.abs(coords[band, 0]) < 4.5]
    print("\n" + "=" * 70)
    print(f"{label}: tip band, {len(band)} nodes, x range [{x.min():.1f}, {x.max():.1f}]")
    print("=" * 70)
    print(f"Tip-band mean Uy (old metric):        {uy.mean():.4f} mm")
    print(f"Regression Uy = a + b*x:  b={slope:+.6f} (mm/mm, twist-like), "
          f"a={intercept:.4f} mm (value at x=0)")
    print(f"Spar-axis nodes (|x|<4.5): {len(axis)} nodes, mean Uy={U[axis,1].mean():.4f} mm")
    print(f"Residual std about the fit: {np.std(uy - (slope*x + intercept)):.4f} mm")
