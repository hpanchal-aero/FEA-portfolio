#!/usr/bin/python3
"""Project 05: figures from results/cases/*_history.csv (no solver run). Writes figures/fig0*.png.

Palette: reference categorical slots 1-4 (#2a78d6, #eb6834, #1baf7a, #eda100) on the light surface #fcfcfb. The dataviz validator
passes the lightness, chroma, CVD and normal-vision gates (adjacent pairs; slots 1-3 also all pairs); aqua and yellow are below
3:1 contrast on the surface, so every figure keeps a legend and distinct line styles, and results/study_summary.md is the table view.
Colour identity: edge condition rest = blue, free = orange, clamp = aqua; thickness 1 / 2 / 3 mm = blue / orange / aqua.
Usage: /usr/bin/python3 scripts/make_figures.py
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from material import MATERIAL

ROOT = os.path.dirname(HERE)
CASES = os.path.join(ROOT, "results", "cases")
FIG = os.path.join(ROOT, "figures")
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
LS = ["-", (0, (5, 2)), (0, (1, 1.6)), (0, (5, 1.5, 1, 1.5))]
MK = ["o", "s", "^", "D"]
LW = 1.8
BCS = (("rest", "Simply supported, in-plane restrained"), ("free", "Simply supported, in-plane free"), ("clamp", "Clamped, in-plane restrained"))
T0 = 20.0

plt.rcParams.update({"figure.facecolor": SURF, "savefig.facecolor": SURF, "axes.facecolor": SURF, "font.size": 10,
                     "font.family": "sans-serif", "text.color": INK})


def cname(h, bc, a, n=20):
    hs = ("%g" % h).replace(".", "p")
    at = ("m" if a < 0 else "") + ("%g" % abs(a)).replace(".", "p")
    return "h%s_%s_a%s_n%d" % (hs, bc, at, n)


def load(h, bc, a, n=20):
    p = os.path.join(CASES, cname(h, bc, a, n) + "_history.csv")
    if not os.path.exists(p):
        raise SystemExit("missing " + p)
    return np.genfromtxt(p, delimiter=",", names=True)


def at(d, col, time):
    k = int(np.argmin(np.abs(d["time_s"] - time)))
    if abs(d["time_s"][k] - time) > 1e-6:
        raise SystemExit("no output at t = %g s" % time)
    return float(d[col][k])


def style(ax, xlabel, ylabel, title=None):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(1.0)
    ax.grid(True, color=GRID, lw=0.7, ls="-")
    ax.set_axisbelow(True)
    ax.tick_params(colors=AXIS, labelcolor=INK2, labelsize=9, length=3)
    ax.set_xlabel(xlabel, color=INK2, fontsize=10)
    ax.set_ylabel(ylabel, color=INK2, fontsize=10)
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="left")


def legend(obj, **kw):
    leg = obj.legend(frameon=False, fontsize=9, **kw)
    for t in leg.get_texts():
        t.set_color(INK2)
    return leg


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print("wrote figures/%s (%d kB)" % (name, os.path.getsize(p) // 1024))


def fig_thermal():
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    for i, h in enumerate((1, 2, 3)):
        d = load(h, "rest", 0.1)
        m = np.isfinite(d["T_top_centre_C"])
        t, T = d["time_s"][m], d["T_top_centre_C"][m]
        ax.plot(t, T, color=C[i], ls=LS[i], lw=LW, label="h = %g mm" % h)
        ax.annotate("%.0f \u00b0C" % T[-1], (t[-1], T[-1]), xytext=(5, 0), textcoords="offset points", color=INK2, va="center", fontsize=9)
    ax.set_xlim(0, 34)
    style(ax, "Time [s]", "Top-face temperature at the plate centre [\u00b0C]", "Gaussian heating pulse, 30 kW/m\u00b2 peak flux")
    legend(ax, loc="upper left")
    fig.tight_layout()
    save(fig, "fig01_thermal_history.png")


def fig_deflection():
    fig, axs = plt.subplots(1, 3, figsize=(12.8, 4.1))
    for ax, h in zip(axs, (1, 2, 3)):
        for i, (bc, lab) in enumerate(BCS):
            d = load(h, bc, 0.1)
            ax.plot(d["time_s"], d["W_centre_over_h"], color=C[i], ls=LS[i], lw=LW, label=lab)
        style(ax, "Time [s]", "Centre deflection w/h [-]" if h == 1 else "", "h = %g mm" % h)
    hd, lb = axs[0].get_legend_handles_labels()
    legend(fig, handles=hd, labels=lb, loc="lower center", ncol=3)
    fig.subplots_adjust(bottom=0.27, wspace=0.27, left=0.06, right=0.99, top=0.9)
    save(fig, "fig02_deflection_history.png")


def fig_stress():
    E, al, nu = MATERIAL["E"], MATERIAL["alpha"], MATERIAL["nu"]
    fig, axs = plt.subplots(1, 2, figsize=(10.6, 4.3))
    print("peak von Mises / fully constrained biaxial thermal stress E alpha dT/(1-nu), at 30 s (dT = top-centre rise):")
    for i, (bc, lab) in enumerate(BCS):
        hs, vm, rel = [], [], []
        for h in (1, 2, 3):
            d = load(h, bc, 0.1)
            dT = at(d, "T_top_centre_C", 30.0) - T0
            v = float(np.max(d["max_vM_MPa"]))
            hs.append(h)
            vm.append(v)
            rel.append(v / (E * al * dT / (1.0 - nu) / 1e6))
        print("  %-6s h = 1, 2, 3 mm: %s" % (bc, ", ".join("%.3f" % r for r in rel)))
        for ax, y in ((axs[0], vm), (axs[1], rel)):
            ax.plot(hs, y, color=C[i], ls=LS[i], lw=LW, marker=MK[i], ms=6.5, mec=SURF, mew=1.3, label=lab)
    for ax in axs:
        ax.set_xticks([1, 2, 3])
        ax.set_xlim(0.8, 3.2)
        ax.set_ylim(bottom=0)
    style(axs[0], "Panel thickness [mm]", "Peak von Mises stress at 30 s [MPa]", "Peak stress")
    style(axs[1], "Panel thickness [mm]", "Peak stress / fully constrained biaxial value [-]", "Fraction of the constrained thermal stress")
    hd, lb = axs[0].get_legend_handles_labels()
    legend(fig, handles=hd, labels=lb, loc="lower center", ncol=3)
    fig.subplots_adjust(bottom=0.25, wspace=0.25, left=0.08, right=0.99, top=0.9)
    save(fig, "fig03_stress_vs_thickness.png")


def fig_edge_force():
    fig, axs = plt.subplots(1, 3, figsize=(12.8, 4.1))
    for ax, h in zip(axs, (1, 2, 3)):
        for i, (bc, lab) in enumerate(BCS):
            if bc == "free":
                continue
            d = load(h, bc, 0.1)
            ax.plot(d["time_s"], np.abs(d["Rx_face_x0_N"]) / 1e3, color=C[i], ls=LS[i], lw=LW, label=lab)
        style(ax, "Time [s]", "Edge force on the face x = 0 [kN]" if h == 1 else "", "h = %g mm" % h)
        ax.set_ylim(bottom=0)
    hd, lb = axs[0].get_legend_handles_labels()
    legend(fig, handles=hd, labels=lb, loc="lower center", ncol=2)
    fig.subplots_adjust(bottom=0.27, wspace=0.27, left=0.06, right=0.99, top=0.9)
    save(fig, "fig04_edge_force_history.png")


def curve(ax, t, w, a, i, lab):
    if a < 0 and np.any(w > 0):
        j = int(np.argmax(w > 0))
        ax.plot(t[:j], w[:j], color=C[i], ls=LS[i], lw=LW, label=lab)
        ax.plot(t[j - 1:j + 1], w[j - 1:j + 1], color=C[i], ls=(0, (1, 2.5)), lw=LW)
        ax.plot(t[j:], w[j:], color=C[i], ls=LS[i], lw=LW)
        return j
    ax.plot(t, w, color=C[i], ls=LS[i], lw=LW, label=lab)
    return None


def fig_imperfection():
    fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.4))
    for ax, tmax in zip(axs, (31.0, 6.2)):
        for i, a in enumerate((-0.1, 0.01, 0.1, 0.5)):
            d = load(2, "rest", a)
            m = d["time_s"] <= tmax
            t, w = d["time_s"][m], d["W_centre_over_h"][m]
            j = curve(ax, t, w, a, i, "A/h = %+g" % a)
            if j is not None and tmax < 10:
                ax.annotate("increment 3.5 \u2192 4.0 s:\nsnap-through, path\nnot resolved", (0.5 * (t[j - 1] + t[j]), 0.5 * (w[j - 1] + w[j])),
                            xytext=(4.15, -0.52), ha="left", va="bottom", color=INK2, fontsize=9, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
        ax.axhline(0, color=AXIS, lw=0.8)
    style(axs[0], "Time [s]", "Centre deflection w/h [-]", "h = 2 mm, restrained: full pulse")
    style(axs[1], "Time [s]", "", "First 6 s")
    axs[1].set_xlim(0, 6.2)
    axs[1].set_ylim(-0.6, 1.5)
    hd, lb = axs[0].get_legend_handles_labels()
    legend(fig, handles=hd, labels=lb, loc="lower center", ncol=4)
    fig.subplots_adjust(bottom=0.24, wspace=0.18, left=0.07, right=0.99, top=0.9)
    save(fig, "fig05_imperfection_sensitivity.png")


def fig_mesh():
    fig, axs = plt.subplots(2, 2, figsize=(10.8, 7.4), sharex=True)
    for col, h in enumerate((1, 2)):
        da, db = load(h, "rest", 0.1, 20), load(h, "rest", 0.1, 40)
        for i, (d, lab) in enumerate(((da, "20 x 20 in-plane elements"), (db, "40 x 40 in-plane elements"))):
            axs[0, col].plot(d["time_s"], d["W_centre_over_h"], color=C[i], ls=LS[i], lw=LW, label=lab)
            m = np.isfinite(d["T_top_centre_C"])
            axs[1, col].plot(d["time_s"][m], d["T_top_centre_C"][m], color=C[i], ls=LS[i], lw=LW, label=lab)
        dw = max(abs(at(db, "W_centre_over_h", x) / at(da, "W_centre_over_h", x) - 1) for x in (5, 10, 15, 20, 25, 30))
        dT = at(db, "T_top_centre_C", 30.0) - at(da, "T_top_centre_C", 30.0)
        rise = at(da, "T_top_centre_C", 30.0) - T0
        axs[0, col].text(0.04, 0.94, "largest difference in w/h (5 to 30 s): %.2f %%" % (100 * dw), transform=axs[0, col].transAxes, color=INK2, fontsize=9, va="top")
        axs[1, col].text(0.04, 0.94, "top-face temperature at 30 s: %+.1f K (%.1f %% of the rise)" % (dT, 100 * dT / rise), transform=axs[1, col].transAxes, color=INK2, fontsize=9, va="top")
        style(axs[0, col], "", "Centre deflection w/h [-]" if col == 0 else "", "h = %g mm" % h)
        style(axs[1, col], "Time [s]", "Top-face temperature at the centre [\u00b0C]" if col == 0 else "")
    hd, lb = axs[0, 0].get_legend_handles_labels()
    legend(fig, handles=hd, labels=lb, loc="lower center", ncol=2)
    fig.subplots_adjust(bottom=0.13, hspace=0.12, wspace=0.2, left=0.09, right=0.99, top=0.95)
    save(fig, "fig06_mesh_sensitivity.png")


def main():
    os.makedirs(FIG, exist_ok=True)
    fig_thermal()
    fig_deflection()
    fig_stress()
    fig_edge_force()
    fig_imperfection()
    fig_mesh()


if __name__ == "__main__":
    main()
