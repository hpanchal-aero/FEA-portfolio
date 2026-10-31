#!/usr/bin/python3
"""Tracer-particle animation of a steady CFD solution, two cases stacked, fixed face-on view at slice y = 6 m.

usage:
  make_particle_video.py preview <out_prefix> "label1=internal1.vtk" "label2=internal2.vtk"
  make_particle_video.py video   <out.mp4>    "label1=internal1.vtk" "label2=internal2.vtk"

The solution is steady, so a particle path is a streamline of the converged velocity field U.  Streamlines are
integrated forward from seed points 14 m upstream of the leading edge (61 heights x 3 spanwise offsets), each carrying
its integration time.  Particles are re-emitted periodically (NP phases per streamline), drawn black with short
fading trails over the Mach slice (coolwarm 0.6..1.4, thick line Ma 1.0, thin lines Ma 1.2 and 1.4).
Only particles within SLAB of the slice plane are drawn (projected onto it) and particles inside the wing section
are hidden.  Time scale: 1 s of video = SLOW s of flow.  The animation is a visualisation of a steady field, not a
time-accurate result.
"""
import os
import sys
import time

import imageio.v2 as imageio
import numpy as np
import pyvista as pv
from matplotlib.path import Path

pv.OFF_SCREEN = True
W, H = 1920, 1080
FPS = 30
DURATION = 12            # s of video
YS = 6.0                 # slice position, m
HALF_W = 20.0            # half panel width, m
HALF_H = HALF_W * (H / 2) / W
CLIM = (0.6, 1.4)
OFF = np.array([0.0, -0.05, 0.0])   # toward the camera (camera sits at y - 100)
SLOW = 0.02              # s of flow per s of video
TCYC = 0.16              # s, particle lifetime / re-emission period
MAXLEN = 70.0             # m, maximum streamline length
NP = 8                   # particles per streamline
NZ = 72                  # seed heights (spacing about 0.16 m)
ZSPAN = 5.0              # m, seeds span +-ZSPAN about the section centre
YOFFS = (-0.5, 0.0, 0.5)
UPSTREAM = 14.0          # m ahead of the leading edge
SLAB = 1.5               # m, half thickness of the displayed slab
XCLIP = 16.0             # m, particles drawn only left of xc + XCLIP (colour bar)
ZLOW = 1.3               # m, no particles within ZLOW of the panel bottom (legend)
ZHIGH = 1.7              # m, no particles within ZHIGH of the panel top (titles)
TRAIL_DT = 0.004         # s of flow between trail dots
NTRAIL = 5
SIZES = [7, 6, 5, 4, 3.5, 3]
OPAC = [1.0, 0.55, 0.40, 0.30, 0.20, 0.12]


def load(path):
    m = pv.read(path)
    for k in list(m.cell_data.keys()):
        del m.cell_data[k]
    for k in list(m.point_data.keys()):
        if k not in ("Ma", "U"):
            del m.point_data[k]
    base = os.path.dirname(path)
    t = os.path.basename(path)[:-4].rsplit("_", 1)[1]
    wing = pv.read(os.path.join(base, "wing", f"wing_{t}.vtk"))
    return m, wing


def prep(label, path):
    m, wing = load(path)
    wl = wing.slice(normal=(0.0, 1.0, 0.0), origin=(0.0, YS, 0.0))
    b = wl.bounds
    xc = 0.5 * (b[0] + b[1])
    zc0 = 0.5 * (b[4] + b[5])
    zc = zc0 + 0.2 * HALF_H
    p = np.asarray(wl.points)
    c = p.mean(axis=0)
    order = np.argsort(np.arctan2(p[:, 2] - c[2], p[:, 0] - c[0]))
    poly = Path(p[order][:, [0, 2]])
    sl = m.slice(normal=(0.0, 1.0, 0.0), origin=(0.0, YS, 0.0))
    sl = sl.clip_box(bounds=(xc - HALF_W - 3.0, xc + HALF_W + 3.0, YS - 1.0, YS + 1.0,
                             zc - HALF_H - 3.0, zc + HALF_H + 3.0), invert=False)
    son = sl.contour(isosurfaces=[1.0], scalars="Ma")
    hi = sl.contour(isosurfaces=[1.2, 1.4], scalars="Ma")
    for g in (son, hi):
        if g.n_points:
            g.points = g.points + OFF
    wl.points = wl.points + OFF
    zs = np.linspace(zc - HALF_H - 0.5, zc + HALF_H + 2.5, NZ)
    seeds = np.array([[b[0] - UPSTREAM, YS + dy, z] for dy in YOFFS for z in zs])
    st = m.streamlines_from_source(pv.PolyData(seeds), vectors="U", integration_direction="forward",
                                   max_length=MAXLEN, initial_step_length=0.5, max_steps=100000,
                                   terminal_speed=1.0, compute_vorticity=False)
    tt = np.asarray(st.point_data["IntegrationTime"])
    P = np.asarray(st.points)
    lines = np.asarray(st.lines)
    paths = []
    i = 0
    while i < len(lines):
        n = int(lines[i])
        ids = lines[i + 1:i + 1 + n]
        i += n + 1
        if n >= 3:
            t = tt[ids] - tt[ids][0]
            paths.append((t, P[ids].copy()))
    print(f"{label}: {len(paths)} streamlines of {len(seeds)} seeds, "
          f"median end time {np.median([q[0][-1] for q in paths]) * 1000:.0f} ms, "
          f"max Mach on slice {np.asarray(sl.point_data['Ma']).max():.2f}", flush=True)
    return dict(label=label, sl=sl, son=son, hi=hi, wl=wl, xc=xc, zc=zc, poly=poly, paths=paths,
                mamax=float(np.asarray(sl.point_data["Ma"]).max()))


def particles(c, tau, k):
    """Positions of the particles at trail level k (k = 0: current position)."""
    out = []
    for t, P in c["paths"]:
        ages = (tau + np.arange(NP) * TCYC / NP) % TCYC - k * TRAIL_DT
        ages = ages[(ages >= 0.0) & (ages <= t[-1])]
        if ages.size == 0:
            continue
        out.append(np.column_stack([np.interp(ages, t, P[:, 0]), np.interp(ages, t, P[:, 1]),
                                    np.interp(ages, t, P[:, 2])]))
    if not out:
        return np.zeros((0, 3))
    X = np.vstack(out)
    ok = ((np.abs(X[:, 0] - c["xc"]) < HALF_W) & (np.abs(X[:, 2] - c["zc"]) < HALF_H)
          & (np.abs(X[:, 1] - YS) < SLAB))
    X = X[ok]
    X = X[(X[:, 0] < c["xc"] + XCLIP) & (X[:, 2] > c["zc"] - HALF_H + ZLOW) & (X[:, 2] < c["zc"] + HALF_H - ZHIGH)]
    if len(X):
        X = X[~c["poly"].contains_points(X[:, [0, 2]])]
    X[:, 1] = YS - 0.2
    return X


def frame(cases, tau):
    n = len(cases)
    pl = pv.Plotter(shape=(n, 1), off_screen=True, window_size=(W, H), border=False)
    pl.set_background("white")
    counts = []
    for j, c in enumerate(cases):
        pl.subplot(j, 0)
        last = (j == n - 1)
        sb = dict(title="Mach", vertical=True, height=0.8, position_x=0.93, position_y=0.1, color="black",
                  title_font_size=18, label_font_size=14, n_labels=5) if last else None
        pl.add_mesh(c["sl"], scalars="Ma", cmap="coolwarm", lighting=False, clim=CLIM,
                    show_scalar_bar=last, scalar_bar_args=sb)
        if c["son"].n_points:
            pl.add_mesh(c["son"], color="black", line_width=2.5)
        if c["hi"].n_points:
            pl.add_mesh(c["hi"], color="black", line_width=1.0)
        pl.add_mesh(c["wl"], color="black", line_width=3.0)
        for k in range(NTRAIL, -1, -1):
            X = particles(c, tau, k)
            if k == 0:
                counts.append(len(X))
            if len(X):
                pl.add_points(X, color="black", point_size=SIZES[k], opacity=OPAC[k],
                              render_points_as_spheres=True)
        pl.add_text(f"{c['label']}\nflow tracers, slice y = {YS:.2f} m (particles within {SLAB:.1f} m)"
                    f"    max Mach on slice {c['mamax']:.2f}",
                    font_size=14, color="black", position="upper_left")
        if last:
            pl.add_text(f"flow time {tau * 1000:6.1f} ms    video is {1 / SLOW:.0f}x slower than the flow\n"
                        "steady solution: tracers follow streamlines    "
                        "thick line: Mach 1.0    thin lines: Mach 1.2 and 1.4",
                        font_size=11, color="black", position="lower_left")
        pl.camera.position = (c["xc"], YS - 100.0, c["zc"])
        pl.camera.focal_point = (c["xc"], YS, c["zc"])
        pl.camera.up = (0.0, 0.0, 1.0)
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = HALF_H
    img = pl.screenshot(return_img=True)
    pl.close()
    return img, counts


mode, out = sys.argv[1], sys.argv[2]
specs = [a.split("=", 1) for a in sys.argv[3:]]
if mode not in ("preview", "video") or len(specs) != 2 or any(len(s) != 2 for s in specs):
    sys.exit('usage: make_particle_video.py preview|video <out> "label=internal.vtk" "label=internal.vtk"')
t0 = time.time()
cases = [prep(label, path) for label, path in specs]
print(f"prepared 2 cases in {time.time() - t0:.1f} s", flush=True)

if mode == "preview":
    for tv in (0.0, 0.2, 0.4):
        t1 = time.time()
        img, counts = frame(cases, tv * SLOW)
        name = f"{out}_t{int(round(tv * 10)):02d}.png"
        imageio.imwrite(name, img)
        print(f"video time {tv:3.1f} s (flow {tv * SLOW * 1000:5.1f} ms): {time.time() - t1:4.1f} s, "
              f"particles drawn {counts}  -> {name}", flush=True)
else:
    total = DURATION * FPS
    tmp = out[:-4] + ".tmp.mp4"
    wr = imageio.get_writer(tmp, fps=FPS, codec="libx264", quality=8, macro_block_size=1, pixelformat="yuv420p")
    t1 = time.time()
    for f in range(total):
        img, _ = frame(cases, f / FPS * SLOW)
        wr.append_data(img)
        if (f + 1) % 20 == 0:
            print(f"frame {f + 1}/{total}  elapsed {time.time() - t1:7.1f} s", flush=True)
    wr.close()
    os.replace(tmp, out)
    print(f"done: {total} frames, {out}, {os.path.getsize(out) / 1e6:.1f} MB, {time.time() - t1:.0f} s", flush=True)
