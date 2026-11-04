#!/usr/bin/python3
"""Project 05: thickness-comparison video of the 30 s heating pulse (PyVista, off-screen).

Writes figures/video_thickness_comparison.mp4 and three stills (figures/video_still_t05s.png, _t15s.png, _t30s.png).

Columns: h = 1, 2, 3 mm, simply supported restrained, A/h = +0.1, 20 x 20 x 2 mesh (cases results/cases/h*_rest_a0p1_n20*).
Top row: top-face temperature, plan view, 20-370 C. Bottom row: deformed panel coloured by the largest von Mises stress over both
element layers at each plan position, 0-115 MPa; shape = original coordinates + 10 x (solver displacement + initial imperfection).
The colour scales are those of figures/fig07_fields_30s.png (scripts/make_field_figures.py).
One video frame per solver output (steps 1 to 60, every 0.5 s); nothing is interpolated; the extra outputs written at automatic
cutbacks are ignored. The last frame is held for --hold seconds (repeated identical frames). 1 s of video = 0.5 s x fps of simulated time.
Heat flux readout: peak flux 30 kW/m2 reached by a linear 5 s ramp (P05-D8, Amendment 4), centre of the Gaussian spot.
Needs the kept .dat of each case (scripts/run_case.py --keep-dat) and simulation/p5_thermal_hX_n20/temps/frame_001..060.inp.
Self-check: for every frame the centre deflection, the peak von Mises stress and the top-centre temperature taken from the fields must
agree with results/cases/*_history.csv; the script stops at the first disagreement.
Usage: /usr/bin/python3 scripts/make_thickness_video.py [--fps 8] [--hold 2]
"""
import argparse
import os
import sys

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
import make_field_figures as mff

ROOT, SIM = mff.ROOT, mff.SIM
INK, INK2, SURF, DEF = mff.INK, mff.INK2, mff.SURF, mff.DEF
T_CLIM, VM_CLIM = mff.T_CLIM, mff.VM_CLIM
DT_FRAME, NFRAME = 0.5, 60
Q_PEAK_KW, T_RAMP = 30.0, 5.0
THICK = (1.0, 2.0, 3.0)
STILLS = {10: "t05s", 30: "t15s", 60: "t30s"}
DEG = "\u00b0C"


def element_max(r, eids):
    """Largest von Mises value [MPa] over the integration points of every element, in the order of eids."""
    vm = mff.von_mises(r[:, 2:8]) / 1e6
    e = r[:, 0].astype(int)
    order = np.argsort(e, kind="stable")
    uniq, start = np.unique(e[order], return_index=True)
    mx = np.maximum.reduceat(vm[order], start)
    idx = np.searchsorted(uniq, eids)
    if idx.max() >= len(uniq) or not np.array_equal(uniq[idx], eids):
        raise SystemExit("element ids of the stress block differ from the mesh")
    return mx[idx]


def load_case(h):
    hs = ("%g" % h).replace(".", "p")
    name = "h%s_rest_a0p1_n20" % hs
    base = "p5_h%s_n20" % hs
    hm = h * 1e-3
    mb = os.path.join(ROOT, "mesh", base + ".inp")
    nb = ccxio.read_nodes(mb)
    ni = ccxio.read_nodes(os.path.join(ROOT, "mesh", base + "_a0p1.inp"))
    eids, conn = ccxio.read_elements(mb)
    eids = np.array(eids, dtype=int)
    ids = np.array(sorted(nb))
    pos = dict((int(i), k) for k, i in enumerate(ids))
    xyz = np.array([nb[int(i)] for i in ids])
    dz0 = np.array([ni[int(i)][2] - nb[int(i)][2] for i in ids])
    ecen = np.array([np.mean([nb[int(c)] for c in list(cn)[:8]], axis=0) for cn in conn])
    _, inv = np.unique(np.round(ecen[:, :2], 6), axis=0, return_inverse=True)
    inv = np.asarray(inv).reshape(-1)

    def col_max(cvm):
        g = np.full(inv.max() + 1, -1.0)
        np.maximum.at(g, inv, cvm)
        return g[inv]

    blocks = ccxio.parse_dat(os.path.join(SIM, name, name + ".dat"))
    ub, sb, extra = {}, {}, 0
    for b in blocks:
        k = int(round(b["time"] / DT_FRAME))
        on_grid = 1 <= k <= NFRAME and abs(b["time"] - k * DT_FRAME) < 1e-6
        title = b["title"].lower()
        if title.startswith("displacements") and b["set"] == "NALL":
            if on_grid:
                ub[k] = b
            else:
                extra += 1
        elif title.startswith("stresses") and b["set"] == "PANEL" and on_grid:
            sb[k] = b
    del blocks
    missing = [k for k in range(1, NFRAME + 1) if k not in ub or k not in sb]
    if missing:
        raise SystemExit("%s: outputs missing at steps %s" % (name, missing))

    T = np.empty((NFRAME, len(ids)))
    D = np.empty((NFRAME, len(ids), 3))
    C = np.empty((NFRAME, len(eids)))
    for k in range(1, NFRAME + 1):
        rows = np.array(ub[k]["rows"])
        rows = rows[np.argsort(rows[:, 0])]
        if not np.array_equal(rows[:, 0].astype(int), ids):
            raise SystemExit("%s step %d: node ids of the displacement block differ from the mesh" % (name, k))
        D[k - 1] = np.column_stack([rows[:, 1], rows[:, 2], dz0 + rows[:, 3]])
        C[k - 1] = col_max(element_max(np.array(sb[k]["rows"]), eids))
        tm = mff.read_temps(os.path.join(SIM, "p5_thermal_h%s_n20" % hs, "temps", "frame_%03d.inp" % k))
        if len(tm) != len(ids):
            raise SystemExit("%s frame %d: %d temperatures for %d nodes" % (name, k, len(tm), len(ids)))
        T[k - 1] = [tm[int(i)] for i in ids]

    # self-check of every frame against the stored history
    hist = np.genfromtxt(os.path.join(ROOT, "results", "cases", name + "_history.csv"), delimiter=",", names=True)
    cmid = [i for i in ids if abs(nb[int(i)][0] - 0.25) < 1e-9 and abs(nb[int(i)][1] - 0.25) < 1e-9 and abs(nb[int(i)][2] - hm / 2) < 1e-9]
    ctop = [i for i in ids if abs(nb[int(i)][0] - 0.25) < 1e-9 and abs(nb[int(i)][1] - 0.25) < 1e-9 and abs(nb[int(i)][2] - hm) < 1e-9]
    if len(cmid) != 1 or len(ctop) != 1:
        raise SystemExit("%s: centre nodes not found (%d, %d)" % (name, len(cmid), len(ctop)))
    pm, pt = pos[int(cmid[0])], pos[int(ctop[0])]
    worst = [0.0, 0.0, 0.0]
    bad = []
    for k in range(1, NFRAME + 1):
        t = k * DT_FRAME
        j = int(np.argmin(np.abs(hist["time_s"] - t)))
        if abs(hist["time_s"][j] - t) > 1e-6:
            bad.append("t = %.1f s: no row in the history" % t)
            continue
        dT = abs(T[k - 1, pt] - hist["T_top_centre_C"][j])
        wref = hist["W_centre_over_h"][j] * h
        dw = abs(D[k - 1, pm, 2] * 1e3 - wref)
        dv = abs(C[k - 1].max() / hist["max_vM_MPa"][j] - 1.0)
        worst = [max(worst[0], dT), max(worst[1], dw), max(worst[2], dv)]
        if dT > 0.01 or dw > 1e-4 + 5e-4 * abs(wref) or dv > 1e-4:
            bad.append("t = %.1f s: dT %.4f K, dw %.5f mm, dvM %.2e" % (t, dT, dw, dv))
    if bad:
        raise SystemExit("%s: fields disagree with the stored history at %d frame(s), first: %s" % (name, len(bad), "; ".join(bad[:3])))
    print("%-18s 60 frames, %d extra outputs ignored; history agreement: max |dT| %.4f K, max |dw| %.5f mm, max |dvM|/vM %.1e" % (
        name, extra, worst[0], worst[1], worst[2]))
    if D[:, :, 2].max() * 1e3 * 1.0 > 12.6 + 1e-9 or T.max() > T_CLIM[1] or C.max() > VM_CLIM[1]:
        raise SystemExit("%s: a frame exceeds the colour scale" % name)

    cells = np.hstack([[8] + [pos[int(c)] for c in list(cn)[:8]] for cn in conn]).astype(np.int64)
    grid = pv.UnstructuredGrid(cells, np.full(len(eids), 12, dtype=np.uint8), xyz * 1e3)
    grid.point_data["T_C"] = T[0]
    grid.cell_data["vM_MPa"] = C[0]
    kmax = int(np.argmax(np.abs(D[:, :, 2]).max(axis=1)))
    g1 = grid.copy()
    g1.points = (xyz + DEF * D[kmax]) * 1e3
    s0 = grid.extract_surface(algorithm="dataset_surface")
    s1 = g1.extract_surface(algorithm="dataset_surface")
    pid0 = np.asarray(s0.point_data["vtkOriginalPointIds"])
    pid1 = np.asarray(s1.point_data["vtkOriginalPointIds"])
    cid1 = np.asarray(s1.cell_data["vtkOriginalCellIds"])
    if not np.array_equal(pid0, pid1):
        raise SystemExit("%s: surface point order differs between the flat and the deformed grid" % name)
    e0 = s0.extract_feature_edges(feature_angle=60)
    dist, eidx = cKDTree(s0.points).query(e0.points)
    if dist.max() > 1e-6:
        raise SystemExit("%s: edge points do not coincide with surface points" % name)
    e1 = e0.copy()
    e1.points = np.asarray(s1.points)[eidx]
    return dict(h=h, name=name, xyz=xyz, T=T, D=D, C=C, pm=pm, s0=s0, s1=s1, e0=e0, e1=e1, pid0=pid0, cid1=cid1, eidx=eidx, kmax=kmax)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--hold", type=float, default=2.0)
    a = ap.parse_args()
    cases = [load_case(h) for h in THICK]
    out = os.path.join(ROOT, "figures", "video_thickness_comparison.mp4")

    pl = pv.Plotter(shape=(2, 3), off_screen=True, window_size=[1600, 1008])
    pl.set_background(SURF)
    pl.enable_anti_aliasing("ssaa")
    # bar titles are blank (and distinct, PyVista keys bars by title); the quantity name is a separate label below the bar
    bar_t = dict(title=" ", color=INK, n_labels=5, fmt="%.0f", vertical=False, position_x=0.15, position_y=0.07,
                 width=0.7, height=0.05, label_font_size=16)
    bar_v = dict(title="  ", color=INK, n_labels=5, fmt="%.0f", vertical=False, position_x=0.15, position_y=0.07,
                 width=0.7, height=0.05, label_font_size=16)
    for col, c in enumerate(cases):
        pl.subplot(0, col)
        pl.add_mesh(c["s0"], scalars="T_C", cmap=mff.CM_T, clim=T_CLIM, lighting=False, show_scalar_bar=(col == 1),
                    scalar_bar_args=bar_t if col == 1 else None)
        pl.add_mesh(c["e0"], color=INK2, line_width=1)
        pl.view_xy()
        if col == 1:
            pl.add_text("Top-face temperature [%s]" % DEG, position="lower_edge", font_size=8, color=INK, name="lab_T")
        pl.subplot(1, col)
        pl.add_mesh(c["s1"], scalars="vM_MPa", cmap=mff.CM_B, clim=VM_CLIM, lighting=False, show_scalar_bar=(col == 1),
                    scalar_bar_args=bar_v if col == 1 else None)
        pl.add_mesh(c["e1"], color=INK2, line_width=1)
        pl.view_isometric()
        if col == 1:
            pl.add_text("von Mises, max through thickness [MPa]", position="lower_edge", font_size=8, color=INK, name="lab_V")

    def set_frame(k):
        t = k * DT_FRAME
        i = k - 1
        for col, c in enumerate(cases):
            c["s0"].point_data["T_C"] = c["T"][i][c["pid0"]]
            c["s0"].Modified()
            pts = ((c["xyz"] + DEF * c["D"][i]) * 1e3)[c["pid0"]]
            c["s1"].points = pts
            c["s1"].cell_data["vM_MPa"] = c["C"][i][c["cid1"]]
            c["s1"].Modified()
            c["e1"].points = pts[c["eidx"]]
            c["e1"].Modified()
            pl.subplot(0, col)
            pl.add_text("h = %g mm\nmax top-face T %.0f %s" % (c["h"], c["T"][i].max(), DEG), position="upper_left", font_size=12,
                        color=INK, name="hdr0_%d" % col)
            pl.subplot(1, col)
            pl.add_text("h = %g mm\ncentre deflection %.1f mm\npeak von Mises %.1f MPa" % (
                c["h"], c["D"][i, c["pm"], 2] * 1e3, c["C"][i].max()), position="upper_left", font_size=11, color=INK,
                name="hdr1_%d" % col)
        pl.subplot(0, 0)
        pl.add_text("t = %4.1f s" % t, position="lower_left", font_size=13, color=INK, name="clock")
        pl.subplot(0, 2)
        pl.add_text("heat flux at the spot centre\n%.1f kW/m\u00b2" % (Q_PEAK_KW * min(t / T_RAMP, 1.0)), position="upper_right",
                    font_size=11, color=INK, name="flux")
        pl.subplot(1, 0)
        pl.add_text("Flux ramps to its peak at 5 s, then soaks to 30 s.\n"
                    "Shape x%g (solver displacement + imperfection).\n"
                    "One frame per 0.5 s solver output, no interpolation.\n"
                    "Video: 1 s = %.1f s simulated." % (DEF, DT_FRAME * a.fps), position="lower_left", font_size=11, color=INK2,
                    name="note")

    set_frame(1)
    # camera: set from the frame of largest deflection of each case so that nothing is clipped later
    for col, c in enumerate(cases):
        pts = ((c["xyz"] + DEF * c["D"][c["kmax"]]) * 1e3)[c["pid0"]]
        c["s1"].points = pts
        c["s1"].Modified()
        pl.subplot(1, col)
        pl.view_isometric()
    pl.subplot(1, 0)
    cam = pl.camera_position
    for col in range(3):
        pl.subplot(1, col)
        pl.camera_position = cam
        pl.camera.zoom(0.8)
    pl.open_movie(out, framerate=a.fps, quality=8)
    pl.show(auto_close=False)
    for k in range(1, NFRAME + 1):
        set_frame(k)
        pl.write_frame()
        if k in STILLS:
            pl.screenshot(os.path.join(ROOT, "figures", "video_still_%s.png" % STILLS[k]))
    nhold = int(round(a.hold * a.fps))
    for _ in range(nhold):
        pl.write_frame()
    pl.close()
    total = NFRAME + nhold
    print("wrote figures/video_thickness_comparison.mp4: %d frames (%d solver frames + %d hold), %d fps, %.1f s, %d kB" % (
        total, NFRAME, nhold, a.fps, total / a.fps, os.path.getsize(out) // 1024))
    print("stills: " + ", ".join("figures/video_still_%s.png" % v for v in STILLS.values()))


if __name__ == "__main__":
    main()
