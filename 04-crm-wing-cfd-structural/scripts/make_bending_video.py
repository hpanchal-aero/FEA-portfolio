#!/usr/bin/python3
"""Wing bending animation: two structural cases stacked, rear view (looking along -x), true scale.

usage:
  make_bending_video.py preview <out_prefix> "label1=<job prefix>" "label2=<job prefix>"
  make_bending_video.py video   <out.mp4>    "label1=<job prefix>" "label2=<job prefix>"

<job prefix> is the path of a CalculiX job without extension (needs <prefix>.inp and <prefix>.frd).
Deck nodes and S4 connectivity are read from the .inp.  The .frd of a shell model is the EXPANDED 3D model
with renumbered nodes, so each frd node is assigned to its nearest deck shell node by position and the deck-node
displacement is the mean of the assigned frd nodes (same method as contour_displacement.py).
The structure is linear static: the displacement field is scaled by a load factor f in [0, 1]
(shape = jig shape + f * displacement, 1x: no exaggeration).  The animation is NOT a transient.
Colour = f * uz (m) on one common scale 0..umax; grey = unloaded jig shape.
Schedule (30 fps, 360 frames): hold at f = 0 for 0.5 s, ramp up 3.5 s, hold at f = 1 for 2.5 s, ramp down 3.5 s,
hold at f = 0 for 2 s (smooth-step ramps).
"""
import os
import sys
import time

import imageio.v2 as imageio
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

pv.OFF_SCREEN = True
W, H = 1920, 1080
FPS = 30
RMAX = 0.05  # m, maximum distance between a frd node and its deck node
N0, NUP, N1, NDN, N2 = 15, 105, 75, 105, 60


def read_deck(path):
    nodes, elems = {}, []
    section, is_shell = None, False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                parts = [p.strip() for p in line.split(",")]
                section = parts[0].upper()
                is_shell = False
                if section == "*ELEMENT":
                    typ = ""
                    for p in parts[1:]:
                        if p.upper().startswith("TYPE="):
                            typ = p.split("=", 1)[1].upper()
                    is_shell = typ.startswith("S4")
                continue
            tok = [t.strip() for t in line.split(",") if t.strip() != ""]
            if section == "*NODE" and len(tok) >= 4:
                nodes[int(tok[0])] = (float(tok[1]), float(tok[2]), float(tok[3]))
            elif is_shell and len(tok) >= 5:
                elems.append([int(t) for t in tok[1:5]])
    return nodes, np.array(elems, dtype=np.int64)


def read_frd(path):
    coords, disp, block = {}, {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                block = "coord"
            elif line.startswith(" -4  DISP"):
                block = "disp"
            elif line.startswith(" -4  "):
                block = None
            elif line.startswith(" -3"):
                block = None
            elif line.startswith(" -1") and block:
                node = int(line[3:13])
                vals = [float(line[13 + 12 * i:25 + 12 * i]) for i in range(3)]
                (coords if block == "coord" else disp)[node] = vals
    return coords, disp


def build(prefix):
    nodes, E = read_deck(prefix + ".inp")
    coords, disp = read_frd(prefix + ".frd")
    used = np.unique(E)
    miss = [int(n) for n in used if int(n) not in nodes]
    if miss:
        sys.exit(f"ABORT {prefix}: {len(miss)} shell nodes missing in the deck, e.g. {miss[:5]}")
    nod = [i for i in coords if i not in disp]
    if nod:
        sys.exit(f"ABORT {prefix}: {len(nod)} frd nodes without displacement, e.g. {nod[:5]}")
    P0 = np.array([nodes[int(n)] for n in used], dtype=float)
    fid = np.array(sorted(coords), dtype=np.int64)
    FX = np.array([coords[int(i)] for i in fid], dtype=float)
    FU = np.array([disp[int(i)] for i in fid], dtype=float)
    tree = cKDTree(P0)
    d, k = tree.query(FX)
    near = d <= RMAX
    cnt = np.bincount(k[near], minlength=len(P0))
    if (cnt == 0).any():
        bad = np.where(cnt == 0)[0]
        sys.exit(f"ABORT {prefix}: {bad.size} deck shell nodes have no frd node within {RMAX} m, e.g. deck IDs {used[bad][:5].tolist()}")
    U = np.column_stack([np.bincount(k[near], weights=FU[near, c], minlength=len(P0)) / cnt for c in range(3)])
    stats = dict(assigned=int(near.sum()), other=int((~near).sum()), dmax=float(d[near].max()),
                 cmin=int(cnt.min()), cmax=int(cnt.max()))
    lut = np.full(int(used.max()) + 1, -1, dtype=np.int64)
    lut[used] = np.arange(used.size)
    Ei = lut[E]
    faces = np.hstack([np.full((Ei.shape[0], 1), 4, dtype=np.int64), Ei]).ravel()
    return P0, U, faces, E.shape[0], stats


def smooth(s):
    return s * s * (3.0 - 2.0 * s)


def factors():
    f = [0.0] * N0
    f += [smooth((i + 1) / NUP) for i in range(NUP)]
    f += [1.0] * N1
    f += [smooth(1.0 - (i + 1) / NDN) for i in range(NDN)]
    f += [0.0] * N2
    return f


def frame(cases, f, g):
    n = len(cases)
    pl = pv.Plotter(shape=(n, 1), off_screen=True, window_size=(W, H), border=False)
    pl.set_background("white")
    for j, c in enumerate(cases):
        pl.subplot(j, 0)
        last = (j == n - 1)
        mesh = pv.PolyData(c["P0"] + f * c["U"], c["faces"])
        mesh.point_data["uz"] = f * c["U"][:, 2]
        pl.add_mesh(c["ghost"], color="#c9c9c9", lighting=False)
        sb = dict(title="uz (m)", vertical=True, height=0.8, position_x=0.93, position_y=0.1, color="black",
                  title_font_size=18, label_font_size=14, n_labels=6) if last else None
        pl.add_mesh(mesh, scalars="uz", lighting=False, cmap="viridis", clim=(0.0, g["umax"]),
                    show_scalar_bar=last, scalar_bar_args=sb)
        pl.add_text(f"{c['label']}\ntip uz {f * c['tip']:5.3f} m", font_size=14, color="black",
                    position="upper_left")
        pl.add_text(f"load factor {f:4.2f}", font_size=12, color="black", position="upper_right", font="courier")
        if last:
            pl.add_text("static linear solution scaled by the load factor (not a transient)    "
                        "true scale, no exaggeration    grey: unloaded jig shape",
                        font_size=11, color="black", position="lower_left")
        pl.camera.position = (g["xc"] + 150.0, g["yc"], g["zv"])
        pl.camera.focal_point = (g["xc"], g["yc"], g["zv"])
        pl.camera.up = (0.0, 0.0, 1.0)
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = g["half_h"]
    img = pl.screenshot(return_img=True)
    pl.close()
    return img


mode, out = sys.argv[1], sys.argv[2]
specs = [a.split("=", 1) for a in sys.argv[3:]]
if mode not in ("preview", "video") or len(specs) != 2 or any(len(s) != 2 for s in specs):
    sys.exit('usage: make_bending_video.py preview|video <out> "label=job_prefix" "label=job_prefix"')

t0 = time.time()
cases = []
for label, prefix in specs:
    P0, U, faces, ne, s = build(prefix)
    y = P0[:, 1]
    tip = y >= y.max() - 0.001
    tipuz = float(U[tip, 2].mean())
    print(f"{label}: shell nodes {len(P0)}  S4 elements {ne}  y max {y.max():.3f} m  tip uz {tipuz:.4f} m ({int(tip.sum())} nodes)")
    print(f"   frd nodes assigned {s['assigned']}, not assigned {s['other']}, max assignment distance {s['dmax']:.4f} m, "
          f"frd nodes per deck node {s['cmin']}..{s['cmax']}", flush=True)
    cases.append(dict(label=label, P0=P0, U=U, faces=faces, tip=tipuz,
                      ghost=pv.PolyData(P0 + np.array([-0.05, 0.0, 0.0]), faces)))

P0a = cases[0]["P0"]
ymin, ymax = float(P0a[:, 1].min()), float(P0a[:, 1].max())
zlo = min(float((c["P0"][:, 2] + np.minimum(c["U"][:, 2], 0.0)).min()) for c in cases)
zhi = max(float((c["P0"][:, 2] + c["U"][:, 2]).max()) for c in cases)
n = len(cases)
g = dict(xc=0.5 * float(P0a[:, 0].min() + P0a[:, 0].max()), yc=0.5 * (ymin + ymax), zv=0.5 * (zlo + zhi),
         umax=float(np.ceil(max(c["U"][:, 2].max() for c in cases) / 0.1) * 0.1))
g["half_h"] = max(0.5 * (zhi - zlo) + 1.7, 0.5 * (ymax - ymin) / 0.80 * (H / n) / W)
print(f"z extent jig+loaded {zlo:.2f}..{zhi:.2f} m, y {ymin:.2f}..{ymax:.2f} m, colour scale 0..{g['umax']:.1f} m, "
      f"panel half height {g['half_h']:.2f} m (half width {g['half_h'] * W / (H / n):.2f} m)")
print(f"loaded {n} cases in {time.time() - t0:.1f} s", flush=True)

if mode == "preview":
    for f in (0.0, 0.5, 1.0):
        t1 = time.time()
        img = frame(cases, f, g)
        name = f"{out}_f{int(round(f * 10)):02d}.png"
        imageio.imwrite(name, img)
        print(f"load factor {f:3.1f}: {time.time() - t1:4.1f} s  image {img.shape}  -> {name}", flush=True)
else:
    fl = factors()
    total = len(fl)
    tmp = out[:-4] + ".tmp.mp4"
    wr = imageio.get_writer(tmp, fps=FPS, codec="libx264", quality=7, macro_block_size=1, pixelformat="yuv420p")
    t1 = time.time()
    for i, f in enumerate(fl):
        wr.append_data(frame(cases, f, g))
        if (i + 1) % 20 == 0:
            print(f"frame {i + 1}/{total}  load factor {f:4.2f}  elapsed {time.time() - t1:7.1f} s", flush=True)
    wr.close()
    os.replace(tmp, out)
    print(f"done: {total} frames, {out}, {os.path.getsize(out) / 1e6:.1f} MB, {time.time() - t1:.0f} s", flush=True)
