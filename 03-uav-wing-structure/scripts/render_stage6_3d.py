#!/usr/bin/env python3
"""
Project 03, Stage 6: 3D renders for L000 (baseline), L006 (best buckling), L050 (lowest tip deflection at 1.5 mm).
Solid structure (ribs, spar, rib faces): von Mises stress, one colour scale shared by all cases. Free skin: Uy only
(skin stress is not shown). Deformation scale 50x, linear ramp 0..100% of limit load (linear static: a proportional
ramp, not dynamics). Load = ANALYTIC PLACEHOLDER pressure, not CFD. 20 kg and 290 MPa yield are assumptions.
Outputs (figures/3d): <case>_stress.png, <case>_uy.png, <case>_ramp.mp4, <case>.glb (binary glTF, per-vertex colours:
solid = stress on the shared scale, skin = Uy, translucent) ; L006_buckling_mode1.png and a zoom on the point of
maximum modal displacement (L006_buckling_mode1_zoom.png), only with RENDER_BUCKLE=1. These two stills are
NOT published: their peak location (TE rib corner, z = 201 mm) disagrees with modal_span_distribution.py
(x = 85, z = 260 mm), and the cause (probably the skin vertex-to-node mapping in render_buckle) is unresolved.
Usage (system python): /usr/bin/python3 scripts/render_stage6_3d.py [L000,L006,L050]
"""
import json
import os
import struct
import sys
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
from matplotlib import colormaps
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import parse_inp, read_frd, read_frd_nodes, von_mises
from analyze_buckle_modes import read_modes

SCALE = 50.0
OUT = "figures/3d"
os.makedirs(OUT, exist_ok=True)
NAMES = {"L000": "L000 baseline (uniform 200 mm ribs)", "L006": "L006 best buckling (ribs 60/130/200/400)",
         "L050": "L050 lowest tip deflection (ribs 168/799/848/927)"}
CAM = [(-900.0, 700.0, 1900.0), (60.0, 0.0, 450.0), (0.0, 1.0, 0.0)]
SB = dict(width=0.5, height=0.08, position_x=0.25, position_y=0.03, title_font_size=11, label_font_size=10)
ZOOM_OFFSET = np.array([-120.0, 380.0, 260.0])   # camera offset from the focal point for the buckling zoom (first guess)


def bar_label(p, text):
    p.add_text(text, position=(int(0.25 * p.window_size[0]), 6), font_size=10, color="black")


def load_case(lid):
    job = "simulation/stage4_s4_k8_%s_closeedgefacnone_t1_n3p8" % lid
    viewer = "mesh/stage3_skin/stage3_skin_s4_k8_%s_closeedgefacnone_t1_n3p8_viewer.npz" % lid
    xyz, solid, shell, root, cload, t = parse_inp(job + ".inp")
    d = np.load(viewer)
    P, tris, cls = d["P"][1:], d["tris"] - 1, d["tri_cls"]
    dist, idx = cKDTree(xyz).query(P)
    assert dist.max() < 1e-6 and len(np.unique(idx)) == len(P), "vertex mapping gate failed"
    nid = idx + 1
    U = read_frd(job + ".frd", "DISP", 3)
    S = read_frd(job + ".frd", "STRESS", 6)
    vm = von_mises(S)
    fid, fxyz = read_frd_nodes(job + ".frd")
    pos = np.full(max(fid.max(), nid.max()) + 1, -1)
    pos[fid] = np.arange(len(fid))
    is_solid = np.isin(nid, np.unique(solid))
    assert (pos[nid[is_solid]] >= 0).all(), "solid vertex IDs missing from the .frd"
    gap = np.abs(fxyz[pos[nid[is_solid]]] - P[is_solid]).max()
    assert gap < 1e-2, ".frd coordinates do not match the viewer vertices (max %.3e)" % gap
    # skin-only vertices: CalculiX wrote them under expanded-node IDs; use the nearest .frd node
    use = nid.copy()
    sd, si = cKDTree(fxyz).query(P[~is_solid])
    assert sd.max() < 1e-2, "skin vertex has no .frd node within 1e-2 mm (max %.3e)" % sd.max()
    use[~is_solid] = fid[si]
    print("[%s] mapping OK: %d vertices, %d solid (frd coord gap %.2e), %d skin via nearest node (max %.2e mm)"
          % (lid, len(P), is_solid.sum(), gap, (~is_solid).sum(), sd.max()))
    return job, P, tris, cls, nid, U[use], np.where(is_solid, vm[nid], 0.0), is_solid


def surface(P, tris, disp, scale):
    return pv.PolyData(P + scale * disp, np.hstack([np.full((len(tris), 1), 3), tris]).ravel())


def stress_cap(lids):
    """One colour cap for all cases: the largest solid-vertex p99, rounded up to a multiple of 5 MPa."""
    p99 = {}
    for lid in lids:
        _, _, _, _, _, _, vm, is_solid = load_case(lid)
        p99[lid] = float(np.percentile(vm[is_solid], 99))
    cap = float(np.ceil(max(p99.values()) / 5.0) * 5.0)
    print("solid-vertex p99 per case: %s -> shared colour cap %.0f MPa" % (
        ", ".join("%s %.1f" % (k, v) for k, v in p99.items()), cap))
    return cap


def _rgba(values, vmin, vmax, cmap, alpha):
    t = np.clip((values - vmin) / (vmax - vmin), 0.0, 1.0)
    c = colormaps[cmap](t)
    c[:, 3] = alpha
    return c.astype(np.float32)


def write_glb(path, parts):
    """Binary glTF 2.0. parts: list of (positions (n,3), triangles (m,3) local indices, rgba (n,4), alphaMode)."""
    blob = bytearray()
    views, accs, prims, mats = [], [], [], []

    def add(arr, target):
        while len(blob) % 4:
            blob.append(0)
        off = len(blob)
        blob.extend(arr.tobytes())
        views.append(dict(buffer=0, byteOffset=off, byteLength=int(arr.nbytes), target=target))
        return len(views) - 1

    for P, tris, rgba, mode in parts:
        P = np.ascontiguousarray(P, dtype=np.float32)
        rgba = np.ascontiguousarray(rgba, dtype=np.float32)
        idx = np.ascontiguousarray(tris, dtype=np.uint32).ravel()
        a_pos = len(accs)
        accs.append(dict(bufferView=add(P, 34962), componentType=5126, count=len(P), type="VEC3",
                         min=[float(v) for v in P.min(axis=0)], max=[float(v) for v in P.max(axis=0)]))
        a_col = len(accs)
        accs.append(dict(bufferView=add(rgba, 34962), componentType=5126, count=len(rgba), type="VEC4"))
        a_idx = len(accs)
        accs.append(dict(bufferView=add(idx, 34963), componentType=5125, count=len(idx), type="SCALAR"))
        mats.append(dict(pbrMetallicRoughness=dict(baseColorFactor=[1.0, 1.0, 1.0, 1.0], metallicFactor=0.0,
                                                   roughnessFactor=1.0), doubleSided=True, alphaMode=mode))
        prims.append(dict(attributes=dict(POSITION=a_pos, COLOR_0=a_col), indices=a_idx, mode=4,
                          material=len(mats) - 1))
    gltf = dict(asset=dict(version="2.0", generator="FEA-portfolio scripts/render_stage6_3d.py"), scene=0,
                scenes=[dict(nodes=[0])], nodes=[dict(mesh=0)], meshes=[dict(primitives=prims)], materials=mats,
                accessors=accs, bufferViews=views, buffers=[dict(byteLength=len(blob))])
    js = json.dumps(gltf, separators=(",", ":")).encode()
    js += b" " * (-len(js) % 4)
    bn = bytes(blob) + b"\0" * (-len(blob) % 4)
    total = 12 + 8 + len(js) + 8 + len(bn)
    with open(path, "wb") as f:
        f.write(struct.pack("<4sII", b"glTF", 2, total))
        f.write(struct.pack("<I4s", len(js), b"JSON"))
        f.write(js)
        f.write(struct.pack("<I4s", len(bn), b"BIN\0"))
        f.write(bn)
    # self-check: re-read what was written
    b = open(path, "rb").read()
    magic, ver, length = struct.unpack("<4sII", b[:12])
    clen, ctype = struct.unpack("<I4s", b[12:20])
    j = json.loads(b[20:20 + clen])
    assert magic == b"glTF" and ver == 2 and length == len(b) and ctype == b"JSON", "GLB header check failed"
    assert len(j["meshes"][0]["primitives"]) == len(parts), "GLB primitive count mismatch"
    print("   GLB check OK: %s, %.2f MB, %d primitives, %d vertices" % (
        path, len(b) / 1e6, len(parts), sum(a["count"] for a in j["accessors"] if a["type"] == "VEC3")))


def glb_parts(P, tris, U, vm, is_solid, tri_solid, cap):
    Pd = P + SCALE * U
    umax = float(U[:, 1].max())
    parts = []
    for mask, kind in ((tri_solid, "solid"), (~tri_solid, "skin")):
        sub = tris[mask]
        if len(sub) == 0:
            continue
        u, inv = np.unique(sub.ravel(), return_inverse=True)
        loc = inv.reshape(-1, 3)
        if kind == "solid":
            parts.append((Pd[u], loc, _rgba(vm[u], 0.0, cap, "turbo", 1.0), "OPAQUE"))
        else:
            parts.append((Pd[u], loc, _rgba(U[u, 1], 0.0, umax, "coolwarm", 0.35), "BLEND"))
    return parts


def render_case(lid, cap):
    job, P, tris, cls, nid, U, vm, is_solid = load_case(lid)
    tri_solid = is_solid[tris].all(axis=1)
    print("[%s] solid-vertex stress: max %.1f MPa, p99 %.1f; shared colour cap %.0f MPa; tip Uy max %.3f mm"
          % (lid, vm[is_solid].max(), np.percentile(vm[is_solid], 99), cap, U[:, 1].max()))
    dte = 0.7 * (300.0 - 0.15 * P[:, 2]) - P[:, 0]
    ok = is_solid & (dte >= 3.0) & (P[:, 2] > 10.0)
    print("[%s] solid-vertex stress excluding the TE zone (dTE < 3 mm) and root (z <= 10 mm): %d vertices, max %.1f MPa, p99 %.1f MPa"
          % (lid, ok.sum(), vm[ok].max(), np.percentile(vm[ok], 99)))
    sk = np.where(is_solid, vm, np.nan)

    def mesh_for(frac):
        m = surface(P, tris, U * frac, SCALE)
        m.point_data["stress"] = sk
        m.point_data["uy"] = U[:, 1] * frac
        return m

    def scene(p, frac, field):
        m = mesh_for(frac)
        if field == "stress":
            sol = m.extract_cells(np.where(tri_solid)[0])
            skin = m.extract_cells(np.where(~tri_solid)[0])
            p.add_mesh(skin, scalars="uy", cmap="coolwarm", clim=(0.0, float(U[:, 1].max())), opacity=0.35, show_scalar_bar=False)
            p.add_mesh(sol, scalars="stress", cmap="turbo", clim=(0, cap), nan_color="lightgrey",
                       scalar_bar_args=dict(SB, title=" "))
            bar_label(p, "von Mises, solid parts [MPa] (capped at %.0f, same scale for all cases)" % cap)
        else:
            p.add_mesh(m, scalars="uy", cmap="coolwarm", clim=(0.0, float(U[:, 1].max())),
                       scalar_bar_args=dict(SB, title=" "))
            bar_label(p, "Uy [mm] (limit load)")
        p.add_text("%s\n%gx deformation, %d%% of limit load; placeholder load; skin stress not shown\ntrailing-edge rib-tip peaks exceed the scale (mesh-dependent singularity, excluded from design stress)"
                   % (NAMES[lid], SCALE, round(frac * 100)), font_size=9, position="upper_left", color="black")
        p.camera_position = CAM
        p.camera.view_angle = 20.0

    for field in ("stress", "uy"):
        p = pv.Plotter(off_screen=True, window_size=(1600, 900))
        p.set_background("white")
        scene(p, 1.0, field)
        p.show(screenshot="%s/%s_%s.png" % (OUT, lid, field))
        print("[%s] wrote %s/%s_%s.png" % (lid, OUT, lid, field))

    p = pv.Plotter(off_screen=True, window_size=(1280, 720))
    p.open_movie("%s/%s_ramp.mp4" % (OUT, lid), framerate=24, quality=8)
    for k in range(48):
        frac = min(k / 36.0, 1.0)
        p.clear()
        p.set_background("white")
        scene(p, frac, "stress")
        p.write_frame()
    p.close()
    print("[%s] wrote %s/%s_ramp.mp4" % (lid, OUT, lid))

    write_glb("%s/%s.glb" % (OUT, lid), glb_parts(P, tris, U, vm, is_solid, tri_solid, cap))


def render_buckle():
    job = "simulation/stage4_s4_k8_L006_closeedgefacnone_t1_n3p8"
    factor = json.load(open("results/stage6/records/L006.json"))["metrics"]["buckling_crit"]
    ids, xyz, blocks = read_modes(job + "_buckle.frd")
    blk = blocks[1]
    bid = np.array([b[0] for b in blk])
    bu = np.array([b[1] for b in blk])
    d = np.load("mesh/stage3_skin/stage3_skin_s4_k8_L006_closeedgefacnone_t1_n3p8_viewer.npz")
    P, tris = d["P"][1:], d["tris"] - 1
    _, idx = cKDTree(parse_inp(job + ".inp")[0]).query(P)
    nid = idx + 1
    full = np.zeros((nid.max() + 1 if nid.max() > bid.max() else bid.max() + 1, 3))
    full[bid] = bu
    u = full[nid]
    amp = np.linalg.norm(u, axis=1).max()
    m = surface(P, tris, u * (50.0 / amp), 1.0)
    m.point_data["mag"] = np.linalg.norm(u, axis=1) / amp
    title = "L006 buckling mode 1, factor %.4f; mode shape, arbitrary amplitude (max 50 mm shown)" % factor
    imax = int(np.argmax(np.linalg.norm(u, axis=1)))
    focal = P[imax]
    print("L006 buckling: factor %.4f (from record); max modal displacement at undeformed x=%.1f y=%.1f z=%.1f mm"
          % (factor, focal[0], focal[1], focal[2]))
    for name, cam, txt in (
            ("L006_buckling_mode1.png", CAM, title),
            ("L006_buckling_mode1_zoom.png", [tuple(focal + ZOOM_OFFSET), tuple(focal), (0.0, 1.0, 0.0)],
             title + "\nzoom on the point of maximum modal displacement (x=%.0f, z=%.0f mm)" % (focal[0], focal[2]))):
        p = pv.Plotter(off_screen=True, window_size=(1600, 900))
        p.set_background("white")
        p.add_mesh(m, scalars="mag", cmap="viridis", scalar_bar_args=dict(SB, title=" "))
        bar_label(p, "displacement magnitude / its maximum (buckling mode 1)")
        p.add_text(txt, font_size=9, position="upper_left", color="black")
        p.camera_position = cam
        if cam is CAM:
            p.camera.view_angle = 20.0
        p.show(screenshot="%s/%s" % (OUT, name))
        print("wrote %s/%s" % (OUT, name))


def main():
    cases = sys.argv[1].split(",") if len(sys.argv) > 1 else ["L000", "L006", "L050"]
    cap = stress_cap(cases)
    for lid in cases:
        render_case(lid, cap)
    if "L006" in cases and os.environ.get("RENDER_BUCKLE") == "1":
        render_buckle()


if __name__ == "__main__":
    main()
