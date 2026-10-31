#!/usr/bin/python3
"""Vertical-displacement contours on the deformed S4 wing-box shell mesh (three structural cases).

usage: contour_displacement.py <out.png> "label1=<job prefix>" ["label2=<job prefix>" ...]
<job prefix> is the path of a CalculiX job without extension (needs <prefix>.inp and <prefix>.frd).

Deck nodes and S4 connectivity come from the .inp.  The .frd of a shell model is the EXPANDED 3D model
with renumbered nodes (IDs unrelated to the deck), so each frd node is assigned to its nearest deck
shell node by position (undeflected coordinates) and the deck-node displacement is the mean of the
assigned frd nodes (top/bottom offsets are antisymmetric, so the mean is the mid-surface value).
Row 1: planform view (from above); row 2: rear view (looking along -x, y to the right), true scale.
Colour = uz in metres on one common scale; black lines every 0.25 m.  The shape is the jig shape plus
the computed displacement (1x, no exaggeration).
"""
import sys

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

pv.OFF_SCREEN = True
STEP = 0.25
RMAX = 0.05  # m, maximum distance between a frd node and its deck node


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
    dev = np.linalg.norm(FU[near] - U[k[near]], axis=1)
    nn = tree.query(P0, k=2)[0][:, 1]
    stats = dict(assigned=int(near.sum()), other=int((~near).sum()), dmax=float(d[near].max()),
                 cmin=int(cnt.min()), cmed=float(np.median(cnt)), cmax=int(cnt.max()),
                 devp99=float(np.percentile(dev, 99)), devmax=float(dev.max()), nnmin=float(nn.min()))
    lut = np.full(int(used.max()) + 1, -1, dtype=np.int64)
    lut[used] = np.arange(used.size)
    Ei = lut[E]
    faces = np.hstack([np.full((Ei.shape[0], 1), 4, dtype=np.int64), Ei]).ravel()
    return P0, U, faces, E.shape[0], stats


out = sys.argv[1]
cases = [a.split("=", 1) for a in sys.argv[2:]]
if not cases or any(len(c) != 2 for c in cases):
    sys.exit('usage: contour_displacement.py <out.png> "label=job_prefix" ...')

built = []
for label, prefix in cases:
    P0, U, faces, ne, s = build(prefix)
    y = P0[:, 1]
    tip = y >= y.max() - 0.001
    k = int(np.argmax(U[:, 2]))
    print(f"{label}: shell nodes {len(P0)}  S4 elements {ne}  y max {y.max():.3f} m")
    print(f"   frd nodes assigned {s['assigned']}, not assigned (>{RMAX} m from any shell node) {s['other']}, max assignment distance {s['dmax']:.4f} m")
    print(f"   frd nodes per deck node: min {s['cmin']}  median {s['cmed']:.0f}  max {s['cmax']}   displacement spread within a node: p99 {s['devp99']:.2e} m  max {s['devmax']:.2e} m")
    print(f"   smallest spacing between two deck shell nodes {s['nnmin']:.4f} m")
    print(f"   tip uz mean {U[tip, 2].mean():.4f} m ({int(tip.sum())} nodes)  max uz {U[k, 2]:.4f} m at y {y[k]:.2f}  "
          f"max|ux| {np.abs(U[:, 0]).max():.4f}  max|uy| {np.abs(U[:, 1]).max():.4f}")
    built.append((label, P0, U, faces, float(U[tip, 2].mean())))

umax = float(np.ceil(max(b[2][:, 2].max() for b in built) / 0.1) * 0.1)
levels = list(np.arange(STEP, umax, STEP))
P0a = built[0][1]
xc = 0.5 * (P0a[:, 0].min() + P0a[:, 0].max())
yc = 0.5 * (P0a[:, 1].min() + P0a[:, 1].max())
zc = 0.5 * (P0a[:, 2].min() + P0a[:, 2].max()) + 0.5 * umax
n = len(built)
print(f"common colour scale 0..{umax:.1f} m, isolines at {[round(v, 2) for v in levels]}")

pl = pv.Plotter(shape=(2, n), off_screen=True, window_size=(1000 * n, 1500), border=False)
pl.set_background("white")
for j, (label, P0, U, faces, tipuz) in enumerate(built):
    mesh = pv.PolyData(P0 + U, faces)
    mesh.point_data["uz"] = U[:, 2]
    lines = mesh.contour(isosurfaces=levels, scalars="uz")
    for i in range(2):
        pl.subplot(i, j)
        sb = dict(title="uz (m)", vertical=True, height=0.6, position_x=0.90, position_y=0.2, color="black",
                  title_font_size=18, label_font_size=16, n_labels=6) if (i == 0 and j == n - 1) else None
        pl.add_mesh(mesh, scalars="uz", lighting=False, cmap="viridis", clim=(0.0, umax), show_scalar_bar=sb is not None,
                    scalar_bar_args=sb, smooth_shading=False)
        if lines.n_points:
            ln = lines.copy()
            ln.points = ln.points + (np.array([0.0, 0.0, 0.05]) if i == 0 else np.array([0.05, 0.0, 0.0]))
            pl.add_mesh(ln, color="black", line_width=1.0)
        if i == 0:
            pl.add_text(f"{label}\ntip uz {tipuz:.3f} m", font_size=11, color="black", position="upper_left")
            pl.camera.position = (xc, yc, 150.0)
            pl.camera.focal_point = (xc, yc, 0.0)
            pl.camera.up = (-1.0, 0.0, 0.0)
            pl.camera.parallel_scale = 12.5
        else:
            pl.add_text(f"{label} - rear view, true scale", font_size=11, color="black", position="upper_left")
            pl.camera.position = (xc + 150.0, yc, zc)
            pl.camera.focal_point = (xc, yc, zc)
            pl.camera.up = (0.0, 0.0, 1.0)
            pl.camera.parallel_scale = 12.5
        pl.camera.parallel_projection = True
pl.screenshot(out)
print("saved", out)
