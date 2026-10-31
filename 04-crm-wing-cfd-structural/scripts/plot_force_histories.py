#!/usr/bin/python3
"""Force-coefficient histories (Cl, Cd, Cm vs iteration) of the Project 04 first-order upwind CFD runs.

usage: plot_force_histories.py

Reads simulation/cfd_v1/runs/<tag>/postProcessing/forceCoeffsCompressible/0/forceCoeffs.dat (columns Time Cm Cd Cl Cl(f) Cl(r)).
Writes figures/force_history_mesh_study.png (L8 / L9 / L10 at alpha 0.36 deg; common iteration 1148 and window (948, 1148])
and figures/force_history_coupling.png (undeflected alpha 0.36, deflected alpha 0.36, deflected alpha 2.07; load-source iterations).
Axis limits are set from iterations >= TCLIP, so the start-up transient is clipped by the axes.
Prints the coefficients at the marked iterations and the Cl range over the preceding 100 iterations (drift check).
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = "simulation/cfd_v1/runs"
TCLIP = 400
TAGS = ["u_L8_p0p36", "u_p0p36", "u_L10_p0p36", "u_defl_p0p36", "u_defl_a2p07"]


def load(tag):
    p = f"{ROOT}/{tag}/postProcessing/forceCoeffsCompressible/0/forceCoeffs.dat"
    if not os.path.exists(p):
        sys.exit(f"ABORT: {p} not found")
    rows = []
    for line in open(p):
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        if len(t) != 6:
            continue
        try:
            rows.append([float(v) for v in t[:4]])
        except ValueError:
            continue
    return np.array(rows)  # Time, Cm, Cd, Cl


DATA = {t: load(t) for t in TAGS}
for t in TAGS:
    d = DATA[t]
    print(f"{t:14s} {len(d)} rows, iterations {int(d[0, 0])}..{int(d[-1, 0])}")


def report(tag, its):
    d = DATA[tag]
    for it in its:
        i = np.where(d[:, 0] == it)[0]
        if i.size == 0:
            print(f"  {tag:14s} iteration {it} not in the file (last {int(d[-1, 0])})")
            continue
        i = int(i[0])
        w = d[max(0, i - 100):i + 1]
        print(f"  {tag:14s} it {it:5d}  Cl {d[i, 3]:.5f}  Cd {d[i, 2]:.5f}  Cm {d[i, 1]:+.5f}   "
              f"Cl range over previous 100 it {w[:, 3].min():.5f}..{w[:, 3].max():.5f}")


def make(out, title, runs, marks, shade, hline, note):
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    handles = []
    for ax, (col, name) in zip(axs, [(3, "Cl"), (2, "Cd"), (1, "Cm")]):
        vals = []
        for tag, label, color in runs:
            d = DATA[tag]
            h, = ax.plot(d[:, 0], d[:, col], color=color, lw=1.2, label=label)
            ax.plot(d[-1, 0], d[-1, col], "o", color=color, ms=4)
            vals.append(d[d[:, 0] >= TCLIP, col])
            if name == "Cl":
                handles.append(h)
        v = np.concatenate(vals)
        pad = 0.12 * (v.max() - v.min())
        ax.set_ylim(v.min() - pad, v.max() + pad)
        if shade:
            ax.axvspan(shade[0], shade[1], color="0.88", zorder=0)
        for it, lab, ls in marks:
            ax.axvline(it, color="k", ls=ls, lw=0.9)
            if name == "Cl":
                ax.text(it, 1.01, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)
        if name == "Cl" and hline is not None:
            ax.axhline(hline, color="0.4", ls="--", lw=1.0)
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
    axs[-1].set_xlabel("iteration")
    axs[-1].set_xlim(0, 1550)
    fig.suptitle(title, fontsize=12)
    fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, 0.955), frameon=False)
    fig.text(0.5, 0.008, note, ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    fig.savefig(out, dpi=200)
    print("saved", out)


print("mesh study, common iteration 1148:")
for t in ("u_L8_p0p36", "u_p0p36", "u_L10_p0p36"):
    report(t, [1148])
make("figures/force_history_mesh_study.png",
     "Force-coefficient histories, mesh study (first-order upwind, alpha = 0.36 deg)",
     [("u_L8_p0p36", "L8, 0.64 M cells (stopped at 1161)", "#1f4e79"),
      ("u_p0p36", "L9, 1.10 M cells (1500 iterations)", "#d98216"),
      ("u_L10_p0p36", "L10, 2.02 M cells (stopped at 1148)", "#2a8a5b")],
     [(1148, "common iteration 1148", ":")], (948, 1148), None,
     "Grey band: common window (948, 1148]. Dots: last iteration of each run. "
     "Axis limits from iterations >= 400: the start-up transient is clipped.")

print("coupling runs, load-source iterations:")
report("u_p0p36", [1500])
report("u_defl_p0p36", [1500])
report("u_defl_a2p07", [1300, 1500])
make("figures/force_history_coupling.png",
     "Force-coefficient histories, undeflected and deflected wing (first-order upwind)",
     [("u_p0p36", "undeflected wing, alpha 0.36 deg", "#1f4e79"),
      ("u_defl_p0p36", "deflected wing, alpha 0.36 deg", "#d98216"),
      ("u_defl_a2p07", "deflected wing, alpha 2.07 deg", "#2a8a5b")],
     [(1300, "it 1300", ":"), (1500, "it 1500", "--")],
     None, 0.5,
     "it 1300: load source of the alpha 2.07 deg run; it 1500: load source of the two alpha 0.36 deg runs. Dashed line: CL = 0.5.\nAxis limits from iterations >= 400: the start-up transient is clipped.")
