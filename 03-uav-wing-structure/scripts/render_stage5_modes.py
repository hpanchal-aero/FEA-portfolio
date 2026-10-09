#!/usr/bin/env python3
"""
Stage 5 figure: critical buckling mode (buckling mode 1 = block 2 of the buckle .frd; block 1 is the static preload),
upper skin, plan view (span to the right, leading edge at the bottom), for rib spacing 200 mm (top) and 100 mm
(bottom). Points colored by |U| relative to each mode's own maximum (same 0-1 scale); black lines = rib planes.
Needs the local buckle .frd files (not in git) and results/stage5_comparison.csv.
Run (conda env, pyvista):  python3 scripts/render_stage5_modes.py
"""
import csv
import os
import sys
import numpy as np
import pyvista as pv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_buckle_modes import read_modes

JOBS = [("200 mm", 200.0, "simulation/stage4_s4_k8_t1_n3p8_buckle"),
        ("100 mm", 100.0, "simulation/stage4_s4_k8_rs100_t1_n3p8_buckle")]


def load_mode(job):
    ids, xyz, blocks = read_modes(job + ".frd")
    blk = blocks[1]
    bid = np.array([b[0] for b in blk])
    U = np.array([b[1] for b in blk])
    o = np.argsort(ids)
    pos = np.searchsorted(ids[o], bid)
    P = xyz[o][pos]
    rel = np.linalg.norm(U, axis=1)
    rel = rel / rel.max()
    sel = (P[:, 1] > 0) & (P[:, 2] < 450.0)
    return P[sel], rel[sel]


def main():
    rows = {r["label"]: r for r in csv.DictReader(open("results/stage5_comparison.csv"))}
    os.makedirs("figures", exist_ok=True)
    p = pv.Plotter(shape=(2, 1), off_screen=True, window_size=(1600, 1300))
    p.set_background("white")
    for i, (label, spacing, job) in enumerate(JOBS):
        P, rel = load_mode(job)
        r = rows[label]
        print(f"{label}: {len(P)} upper-skin nodes with z < 450, peak at z = {P[np.argmax(rel), 2]:.1f} mm")
        p.subplot(i, 0)
        cloud = pv.PolyData(P)
        cloud["relative |U|"] = rel
        p.add_mesh(cloud, scalars="relative |U|", cmap="viridis", clim=(0.0, 1.0), point_size=4,
                   render_points_as_spheres=True, show_scalar_bar=(i == 1),
                   scalar_bar_args=dict(title="relative |U|", color="black", vertical=True,
                                    position_x=0.92, position_y=0.08, height=0.75, width=0.04))
        for zs in np.arange(0.0, 450.1, spacing):
            p.add_mesh(pv.Line((-95.0, 22.0, float(zs)), (215.0, 22.0, float(zs))), color="black", line_width=1.5)
        p.add_text(f"{label} rib spacing ({r['n_ribs']} ribs over the span): critical buckling factor "
                   f"{float(r['buckling_critical']):.2f}; upper skin, plan view, span to the right, LE at the bottom; "
                   f"black lines = rib planes", color="black", font_size=8, position="upper_left")
        p.camera_position = [(60.0, 800.0, 200.0), (60.0, 0.0, 200.0), (1.0, 0.0, 0.0)]
    p.show(screenshot="figures/stage5_buckling_mode1.png")
    print("wrote figures/stage5_buckling_mode1.png")


if __name__ == "__main__":
    main()
