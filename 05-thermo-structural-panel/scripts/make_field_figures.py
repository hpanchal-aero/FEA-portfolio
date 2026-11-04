#!/usr/bin/python3
"""Project 05: 3D field figure at t = 30 s (PyVista, off-screen). Writes figures/fig07_fields_30s.png.

Rows: h = 1, 2, 3 mm simply supported restrained and h = 3 mm clamped (A/h = +0.1, 20 x 20 mesh). Columns: top-face temperature
(plan view), out-of-plane deflection w (deformed shape) and von Mises stress (deformed shape). The stress shown at each plan position
is the largest von Mises value over the integration points of both element layers through the thickness (the visible top faces
alone would show only the top layer). One colour scale per column, shared by all rows (temperature 20-370 C, w 0-12.6 mm, von Mises
0-115 MPa); the scale bars are drawn once, under the bottom row. Lighting is off so that colour means value only. Display mesh: the
8 corner nodes of each element only. Deformed shape = original coordinates + DEF x (solver displacement + initial imperfection), DEF = 10.
Needs the kept .dat of each case (scripts/run_case.py --keep-dat) and the thermal frame file simulation/p5_thermal_hX_n20/temps/frame_060.inp.
Self-check: centre deflection, peak von Mises and peak temperature from the fields must agree with results/cases/*_history.csv and *.txt.
Also prints the largest element von Mises of the top and of the bottom element layer.
Usage: /usr/bin/python3 scripts/make_field_figures.py
"""
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pyvista as pv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio

pv.OFF_SCREEN = True
ROOT = os.path.dirname(HERE)
SIM = os.path.join(ROOT, "simulation")
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
DEF = 10.0
ROWS = ((1.0, "rest"), (2.0, "rest"), (3.0, "rest"), (3.0, "clamp"))
EDGE = {"rest": "simply supported, restrained", "clamp": "clamped"}
T_CLIM, W_CLIM, VM_CLIM = (20.0, 370.0), (0.0, 12.6), (0.0, 115.0)


def ramp(name):
    return LinearSegmentedColormap.from_list(name + "_t", plt.get_cmap(name)(np.linspace(0.15, 1.0, 256)))


CM_T, CM_B = ramp("Oranges"), ramp("Blues")


def von_mises(s):
    sxx, syy, szz, sxy, sxz, syz = [s[:, k] for k in range(6)]
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))


def read_temps(path):
    d = {}
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        a = s.split(",")
        if len(a) == 2:
            d[int(a[0])] = float(a[1])
    return d


def build(h, bc):
    hs = ("%g" % h).replace(".", "p")
    name = "h%s_%s_a0p1_n20" % (hs, bc)
    base = "p5_h%s_n20" % hs
    mb = os.path.join(ROOT, "mesh", base + ".inp")
    nb = ccxio.read_nodes(mb)
    ni = ccxio.read_nodes(os.path.join(ROOT, "mesh", base + "_a0p1.inp"))
    eids, conn = ccxio.read_elements(mb)
    ids = np.array(sorted(nb))
    pos = dict((int(i), k) for k, i in enumerate(ids))
    xyz = np.array([nb[int(i)] for i in ids])
    dz0 = np.array([ni[int(i)][2] - nb[int(i)][2] for i in ids])
    blocks = ccxio.parse_dat(os.path.join(SIM, name, name + ".dat"))
    ub = sorted([b for b in blocks if b["title"].lower().startswith("displacements") and b["set"] == "NALL"], key=lambda b: b["time"])[-1]
    sb = sorted([b for b in blocks if b["title"].lower().startswith("stresses") and b["set"] == "PANEL"], key=lambda b: b["time"])[-1]
    if abs(ub["time"] - 30.0) > 1e-6 or abs(sb["time"] - 30.0) > 1e-6:
        raise SystemExit("%s: last output is not at 30 s" % name)
    rows = np.array(ub["rows"])
    rows = rows[np.argsort(rows[:, 0])]
    if not np.array_equal(rows[:, 0].astype(int), ids):
        raise SystemExit("%s: node ids of the displacement block differ from the mesh" % name)
    disp = np.column_stack([rows[:, 1], rows[:, 2], dz0 + rows[:, 3]])
    r = np.array(sb["rows"])
    vm = von_mises(r[:, 2:8]) / 1e6
    emax = {}
    for e, v in zip(r[:, 0].astype(int), vm):
        if v > emax.get(int(e), -1.0):
            emax[int(e)] = float(v)
    cvm = np.array([emax[int(e)] for e in eids])
    # element layer (from the centroid height) and largest value through the thickness at each plan position
    ecen = np.array([np.mean([nb[int(c)] for c in list(cn)[:8]], axis=0) for cn in conn])
    top = ecen[:, 2] > h * 5e-4
    key = [(round(float(x), 6), round(float(y), 6)) for x, y in ecen[:, :2]]
    colmax = {}
    for k, v in zip(key, cvm):
        if v > colmax.get(k, -1.0):
            colmax[k] = float(v)
    cvm_col = np.array([colmax[k] for k in key])
    print("%-18s largest element von Mises: top layer %.3f MPa, bottom layer %.3f MPa (%d top, %d bottom elements)" % (
        name, cvm[top].max(), cvm[~top].max(), int(top.sum()), int((~top).sum())))
    tmap = read_temps(os.path.join(SIM, "p5_thermal_h%s_n20" % hs, "temps", "frame_060.inp"))
    if len(tmap) != len(ids):
        raise SystemExit("%s: %d temperatures for %d nodes" % (name, len(tmap), len(ids)))
    T = np.array([tmap[int(i)] for i in ids])
    cells = np.hstack([[8] + [pos[int(c)] for c in list(cn)[:8]] for cn in conn]).astype(np.int64)
    grid = pv.UnstructuredGrid(cells, np.full(len(eids), 12, dtype=np.uint8), xyz * 1e3)
    grid.point_data["T_C"] = T
    grid.point_data["w_mm"] = disp[:, 2] * 1e3
    grid.cell_data["vM_MPa"] = cvm_col
    gdef = grid.copy()
    gdef.points = (xyz + DEF * disp) * 1e3
    # self-check against the stored history
    cen = [i for i in ids if abs(nb[int(i)][0] - 0.25) < 1e-9 and abs(nb[int(i)][1] - 0.25) < 1e-9 and abs(nb[int(i)][2] - h * 5e-4) < 1e-9]
    wc = disp[pos[int(cen[0])], 2] * 1e3
    hist = np.genfromtxt(os.path.join(ROOT, "results", "cases", name + "_history.csv"), delimiter=",", names=True)
    txt = open(os.path.join(ROOT, "results", "cases", name + ".txt")).read()
    tpk = float(re.search(r"peak temperature at 30 s ([0-9.]+) C", txt).group(1))
    ref = (hist["W_centre_over_h"][-1] * h, float(np.max(hist["max_vM_MPa"])), tpk)
    got = (wc, float(cvm_col.max()), float(T.max()))
    print("%-18s centre w %.4f mm (history %.4f)   max vM %.3f MPa (history %.3f)   max T %.2f C (report %.2f)" % (name, got[0], ref[0], got[1], ref[1], got[2], ref[2]))
    if abs(got[0] / ref[0] - 1) > 5e-4 or abs(got[1] / ref[1] - 1) > 1e-4 or abs(got[2] - ref[2]) > 0.01:
        raise SystemExit("%s: field data do not agree with the stored history" % name)
    for label, val, clim in (("temperature", got[2], T_CLIM), ("deflection", float(disp[:, 2].max() * 1e3), W_CLIM), ("von Mises", got[1], VM_CLIM)):
        if val > clim[1]:
            raise SystemExit("%s: %s %.2f exceeds the colour scale %s" % (name, label, val, clim))
    return grid, gdef


def main():
    nrow = len(ROWS)
    pl = pv.Plotter(shape=(nrow, 3), off_screen=True, window_size=[1500, 500 * nrow])
    pl.set_background(SURF)
    cols = (("T_C", CM_T, T_CLIM, "Top-face temperature [\u00b0C], 30 s", "%.0f"),
            ("w_mm", CM_B, W_CLIM, "Deflection w [mm], shape x%g" % DEF, "%.1f"),
            ("vM_MPa", CM_B, VM_CLIM, "von Mises [MPa], max through thickness, shape x%g" % DEF, "%.0f"))
    for row, (h, bc) in enumerate(ROWS):
        grid, gdef = build(h, bc)
        s0 = grid.extract_surface(algorithm="dataset_surface")
        s1 = gdef.extract_surface(algorithm="dataset_surface")
        head = "h = %g mm, %s" % (h, EDGE[bc])
        last = row == nrow - 1
        for col, (arr, cmap, clim, label, fmt) in enumerate(cols):
            pl.subplot(row, col)
            surf = s0 if col == 0 else s1
            kw = {}
            if last:
                kw["scalar_bar_args"] = dict(title=" " * (col + 1), color=INK, n_labels=4, fmt=fmt, vertical=False, position_x=0.18,
                                             position_y=0.02, width=0.64, height=0.05, label_font_size=10)
            pl.add_mesh(surf, scalars=arr, cmap=cmap, clim=clim, lighting=False, show_scalar_bar=last, **kw)
            pl.add_mesh(surf.extract_feature_edges(feature_angle=60), color=INK2, line_width=1)
            if col == 0:
                pl.view_xy()
            else:
                pl.view_isometric()
                pl.camera.zoom(0.9)
            pl.add_text(head + "\n" + label, position="upper_left", font_size=9, color=INK)
    out = os.path.join(ROOT, "figures", "fig07_fields_30s.png")
    pl.screenshot(out)
    pl.close()
    print("wrote figures/fig07_fields_30s.png (%d kB)" % (os.path.getsize(out) // 1024))


if __name__ == "__main__":
    main()
