#!/usr/bin/env python3
"""
Stage 4 figures. Run with the conda env python (numpy + matplotlib):
  python3 scripts/plot_stage4_figures.py
Writes figures/stage4_twist_span.png and figures/stage4_load_distribution.png.

Twist: nose-up rotation of the rib planes from results/stage4_twist.csv (least-squares fit of Uy vs x on solid
original nodes; ribs are not rigid, rms fit residual <= 0.0073 mm in all cases). The dashed Bredt-Batho curve is a HAND
ESTIMATE (single closed cell, t = 1.0 mm, ribs and web ignored, shear centre ASSUMED at the spar axis,
torque arm 0.05 c); plausibility only, not a verification (same formulas as estimate_stage4_twist_bredt.py).
Load: the applied lift is an ANALYTIC PLACEHOLDER, not an aerodynamic prediction.
"""
import csv
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPAN = 1000.0
RIB_Z = [0, 200, 400, 600, 800, 1000]
INP = "simulation/stage4_s4_k8_t1_n3p8.inp"


def G(xi):
    xi = np.clip(xi, 0.0, 1.0)
    th = np.arcsin(np.sqrt(xi))
    return (2.0 / np.pi) * (th + np.sin(th) * np.cos(th))


def Qs(s):
    s = np.clip(s, 0.0, 1.0)
    return 0.5 * (s * np.sqrt(1.0 - s * s) + np.arcsin(s))


def bredt_curve(L=372.78, xsc=0.30):
    E, NU, T = 73100.0, 0.33, 1.0
    Gm = E / (2.0 * (1.0 + NU))
    q0 = L / (SPAN * np.pi / 4.0)
    x = np.linspace(0.0, 1.0, 200001)
    yt = 0.6 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    a = 2.0 * np.trapz(yt, x)
    p = 2.0 * np.sum(np.hypot(np.diff(x), np.diff(yt)))
    z = np.linspace(0.0, SPAN, 200001)
    dz = z[1] - z[0]
    c = 300.0 - 0.15 * z
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN) ** 2, 0.0, 1.0))
    GJ = Gm * 4.0 * (a * c ** 2) ** 2 * T / (p * c)
    m = q * (xsc - 0.25) * c
    Tq = (np.cumsum(m[::-1]) * dz)[::-1]
    th = np.cumsum(Tq / GJ) * dz
    return z, th


def read_twist(path):
    out = {}
    for r in csv.DictReader(open(path)):
        out[(r["mesh"], r["case"])] = [float(r[f"theta_z{z}_rad"]) for z in RIB_Z]
    return out


def parse_inp(path):
    ids, xyz, cl = [], [], []
    sec = None
    L = None
    for line in open(path):
        if line.startswith("** LOADCASE"):
            L = float(re.search(r"LIFT_N=([0-9.eE+-]+)", line).group(1))
        if line.startswith("*") and not line.startswith("**"):
            sec = line.strip().upper().split(",")[0]
            continue
        if line.startswith("**") or not line.strip():
            continue
        p = line.split(",")
        if sec == "*NODE":
            ids.append(int(p[0]))
            xyz.append([float(v) for v in p[1:4]])
        elif sec == "*CLOAD":
            cl.append((int(p[0]), int(p[1]), float(p[2])))
    ids = np.array(ids)
    coords = np.zeros((ids.max() + 1, 3))
    coords[ids] = np.array(xyz)
    F = np.zeros((ids.max() + 1, 3))
    for nd, dof, val in cl:
        F[nd, dof - 1] += val
    return L, coords, F


def fig_twist():
    tw = read_twist("results/stage4_twist.csv")
    zb, thb = bredt_curve()
    tip = thb[-1]
    assert abs(tip / 2.999e-4 - 1.0) < 1e-3, f"Bredt-Batho tip {tip:.4e} differs from the printed 2.999e-4"
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    z = np.array(RIB_Z, float)
    sc = 1e4
    ax.plot(z, np.array(tw[("s4_k8", "stage4_full")]) * sc, "o-", color="C3", label="Stage 4 full (skin pressure + gravity), 4/8 mm")
    ax.plot(z, np.array(tw[("s4_k8", "stage4_lift")]) * sc, "s-", color="C0", label="Stage 4 lift only")
    ax.plot(z, np.array(tw[("s4_k8", "stage4_grav")]) * sc, "^-", color="C2", label="Stage 4 gravity only")
    ax.plot(z, np.array(tw[("s4_k8", "stage3_spar_cap")]) * sc, "d-", color="0.35", label="Stage 3 (lift on spar-cap face), 4/8 mm")
    ax.plot(z, np.array(tw[("s4_k6", "stage4_full")]) * sc, "x", color="k", ms=8, label="Stage 4 full, 4/6 mm (mesh check)")
    ax.plot(zb, thb * sc, "--", color="C0", lw=1.2, label="hand estimate, lift only (single cell, shear centre at spar axis)")
    ax.set_xlabel("span position z [mm]")
    ax.set_ylabel(r"nose-up rotation of rib plane  $\theta$  [$10^{-4}$ rad]")
    ax.set_title("Wing twist at the rib planes, limit load (n = 3.8)\n(rms fit residual <= 0.0073 mm in all cases; applied lift is an analytic placeholder)", fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim(-0.1, 5.3)
    ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig("figures/stage4_twist_span.png", dpi=200)
    plt.close(fig)
    print(f"  wrote figures/stage4_twist_span.png (Bredt-Batho tip {tip:.4e} rad)")


def fig_load():
    L, coords, F = parse_inp(INP)
    mk = np.abs(F[:, 1]) > 0.0
    z, x, fy = coords[mk, 2], coords[mk, 0], F[mk, 1]
    tot = fy.sum()
    print(f"  load file {INP}: header lift {L:.4f} N, sum Fy of nodal loads {tot:.6f} N, loaded nodes {int(mk.sum())}")
    q0 = L / (SPAN * np.pi / 4.0)

    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.3))
    ax = axs[0]
    edges = np.arange(0.0, SPAN + 0.01, 25.0)
    zm = (z >= 0.0) & (z <= SPAN)
    hb, _ = np.histogram(z[zm], bins=edges, weights=fy[zm])
    dz = np.diff(edges)
    q_fea = hb / dz
    q_an = q0 * SPAN * (Qs(edges[1:] / SPAN) - Qs(edges[:-1] / SPAN)) / dz
    cen = 0.5 * (edges[1:] + edges[:-1])
    ax.plot(cen, q_an, "-", color="C0", label="analytic elliptical, bin mean")
    ax.plot(cen, q_fea, "o", ms=3.5, color="C3", label="nodal loads in the .inp (25 mm bins)")
    ax.set_xlabel("span position z [mm]")
    ax.set_ylabel("net lift per unit span [N/mm]")
    ax.set_title("Spanwise (half-wing 372.78 N)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    order = np.argsort(z)
    cum = np.cumsum(fy[order]) / tot
    d_span = float(np.max(np.abs(cum - Qs(z[order] / SPAN) / (np.pi / 4.0))))
    print(f"  spanwise: max |cumulative nodal - analytic| = {d_span:.5f}")

    ax = axs[1]
    c = 300.0 - 0.15 * z
    xi = np.clip((x + 0.3 * c) / c, 0.0, 1.0)
    e2 = np.linspace(0.0, 1.0, 51)
    hb2, _ = np.histogram(xi, bins=e2, weights=fy)
    dxi = np.diff(e2)
    f_fea = hb2 / tot / dxi
    f_an = (G(e2[1:]) - G(e2[:-1])) / dxi
    cx = 0.5 * (e2[1:] + e2[:-1])
    ax.plot(cx, f_an, "-", color="C0", label=r"thin-airfoil $(2/\pi)\sqrt{(1-\xi)/\xi}$, bin mean")
    ax.plot(cx, f_fea, "o", ms=3.5, color="C3", label="nodal loads in the .inp (0.02 bins)")
    ax.axvline(0.25, color="0.5", ls=":", lw=1)
    ax.axvline(0.30, color="0.2", ls="--", lw=1)
    ax.text(0.245, 7.5, "25% c\n(centre of\npressure)", fontsize=7, color="0.4", ha="right")
    ax.text(0.305, 5.0, "30% c\n(spar axis)", fontsize=7, color="0.2")
    ax.set_ylim(0, 10)
    ax.set_xlabel(r"chord fraction from leading edge  $\xi$")
    ax.set_ylabel("normalised net lift per unit chord fraction")
    ax.set_title("Chordwise (upper suction + lower pressure, each half)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    o2 = np.argsort(xi)
    cum2 = np.cumsum(fy[o2]) / tot
    d_ch = float(np.max(np.abs(cum2 - G(xi[o2]))))
    print(f"  chordwise: max |cumulative nodal - analytic| = {d_ch:.5f}")
    fig.suptitle("Stage 4 applied load: ANALYTIC PLACEHOLDER (not an aerodynamic prediction)", fontsize=10)
    fig.tight_layout()
    fig.savefig("figures/stage4_load_distribution.png", dpi=200)
    plt.close(fig)
    print("  wrote figures/stage4_load_distribution.png")


def main():
    os.makedirs("figures", exist_ok=True)
    for p in ("results/stage4_twist.csv", INP):
        if not os.path.isfile(p):
            print(f"ERROR: missing {p}")
            sys.exit(1)
    fig_twist()
    fig_load()


if __name__ == "__main__":
    main()
