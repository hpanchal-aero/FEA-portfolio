#!/usr/bin/env python3
"""
Project 03, Stage 6: 3D renders for L000 (baseline), L006 (best buckling), L050 (lowest tip deflection at 1.5 mm).
Solid structure (ribs, spar, rib faces): von Mises stress. Free skin: Uy only (skin stress is not shown).
Deformation scale 50x, linear ramp 0..100% of limit load (linear static: a proportional ramp, not dynamics).
Load = ANALYTIC PLACEHOLDER pressure, not CFD. 20 kg and 290 MPa yield are assumptions.
Outputs: figures/3d/<case>_stress.png, _uy.png, <case>_ramp.mp4, <case>.glb ; L006 buckling mode still.
Usage (system python): /usr/bin/python3 scripts/render_stage6_3d.py [L000,L006,L050]
"""
import os
import sys
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import parse_inp, read_frd, read_frd_nodes, von_mises
from analyze_buckle_modes import read_modes

SCALE = 50.0
OUT = "figures/3d"
os.makedirs(OUT, exist_ok=True)
NAMES = {"L000": "L000 baseline (uniform 200 mm ribs)", "L006": "L006 best buckling (ribs 60/130/200/400)",
         "L050": "L050 lowest tip deflection (ribs 168/799/848/927)"}


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


def render_case(lid):
    job, P, tris, cls, nid, U, vm, is_solid = load_case(lid)
    tri_solid = is_solid[tris].all(axis=1)
    cap = float(np.percentile(vm[is_solid], 99))
    cap = round(cap, 0) if cap > 10 else cap
    print("[%s] solid-vertex stress: max %.1f MPa, p99 %.1f -> colour cap %.0f MPa; tip Uy max %.3f mm"
          % (lid, vm[is_solid].max(), np.percentile(vm[is_solid], 99), cap, U[:, 1].max()))
    sk = np.where(is_solid, vm, np.nan)
    cam = [(-900.0, 700.0, 1900.0), (60.0, 0.0, 450.0), (0.0, 1.0, 0.0)]

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
                       scalar_bar_args=dict(title="von Mises, solid parts [MPa] (capped at %.0f)" % cap))
        else:
            p.add_mesh(m, scalars="uy", cmap="coolwarm", clim=(0.0, float(U[:, 1].max())),
                       scalar_bar_args=dict(title="Uy [mm] (limit load)"))
        p.add_text("%s\n%gx deformation, %d%% of limit load; placeholder load; skin stress not shown"
                   % (NAMES[lid], SCALE, round(frac * 100)), font_size=9, position="upper_left", color="black")
        p.camera_position = cam

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

    m = mesh_for(1.0)
    pl = pv.Plotter(off_screen=True)
    pl.add_mesh(m, scalars="uy", cmap="coolwarm")
    try:
        pl.export_gltf("%s/%s.glb" % (OUT, lid))
        print("[%s] wrote %s/%s.glb" % (lid, OUT, lid))
    except Exception as e:
        print("[%s] GLB export FAILED: %s" % (lid, e))
    pl.close()


def render_buckle():
    job = "simulation/stage4_s4_k8_L006_closeedgefacnone_t1_n3p8"
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
    amp = np.abs(u).max()
    m = surface(P, tris, u * (50.0 / amp), 1.0)
    m.point_data["mag"] = np.linalg.norm(u, axis=1) / amp
    p = pv.Plotter(off_screen=True, window_size=(1600, 900))
    p.set_background("white")
    p.add_mesh(m, scalars="mag", cmap="viridis", scalar_bar_args=dict(title="|U| / max (buckling mode 1)"))
    p.add_text("L006 buckling mode 1, factor 4.7174; mode shape, arbitrary amplitude (max 50 mm shown)",
               font_size=9, position="upper_left", color="black")
    p.camera_position = [(-900.0, 700.0, 1900.0), (60.0, 0.0, 450.0), (0.0, 1.0, 0.0)]
    p.show(screenshot="%s/L006_buckling_mode1.png" % OUT)
    print("wrote %s/L006_buckling_mode1.png" % OUT)


def main():
    cases = sys.argv[1].split(",") if len(sys.argv) > 1 else ["L000", "L006", "L050"]
    for lid in cases:
        render_case(lid)
    if "L006" in cases:
        render_buckle()


if __name__ == "__main__":
    main()
