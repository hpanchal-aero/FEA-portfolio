#!/usr/bin/env python3
"""
Project 03, Stage 3 -- build the CalculiX .inp for spar v2 + ribs (C3D10) + skin (S6).

Input : raw Netgen mesh .npz (mesh_stage3_skin_raw.py) and the skin thickness [mm].
Steps : classify surface tris (match face-descriptor idx to stage3_skin_faces.json by area +
        centroid) -> weld coincident nodes (1e-6 mm; tet node ID survives) -> verify skin shell
        mesh (no degenerates, edge-manifold, boundary only at z-ends, area == lateral area,
        quality) -> outward winding (n . (p - (0,0,z)) > 0) + directed-edge consistency ->
        contiguous renumbering -> one global midside-node table for C3D10 + S6 ->
        counts vs ceilings -> root BC (all nodes at z_min), slab-by-slab placeholder lift ->
        write .inp (+ viewer .npz).

Skin convention: nodes lie ON the airfoil envelope, material extends OUTWARD
(*SHELL SECTION OFFSET=-0.5, outward normals) => outer mold line = nominal + t.
Placeholder load (24.5 N elliptical lift on the spar top-cap face, + self-weight) as in
Stages 1/2; NOT a sourced design requirement.

Usage (conda env ccx223):
  python3 scripts/build_stage3_skin_inp.py mesh/stage3_skin/raw_s2p5_k8.npz 1.0
"""
import json
import math
import os
import sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

SPAN = 1000.0
ROOT_CHORD, TIP_CHORD = 300.0, 150.0
SPAR_CHORD_FRAC = 0.30
N_AIRFOIL_PTS = 40
CAP_W = 8.0
import re as _re_rs
import sys as _sys_rs
_base_rs = _sys_rs.argv[1].split("/")[-1] if len(_sys_rs.argv) > 1 else ""
_m_tag = _re_rs.search(r"((?:_rs\d+|_L[0-9A-Za-z]+)(?:_rt\d+(?:p\d+)?)?|_rt\d+(?:p\d+)?)", _base_rs)
RS_SUFFIX = _m_tag.group(1) if _m_tag else ""
_m_rt = _re_rs.search(r"_rt(\d+)(?:p(\d+))?", RS_SUFFIX)
RIB_T = float(_m_rt.group(1) + ("." + _m_rt.group(2) if _m_rt.group(2) else "")) if _m_rt else 1.5
E_MOD, NU, RHO, G = 73100.0, 0.33, 2780e-12, 9810.0
TOTAL_LIFT_N = 24.5            # PLACEHOLDER
TOL = 1e-3
WELD = 1e-6
EQ_CEILING, ELEM_CEILING = 1_180_000, 230_000
SOLID_VOL = 128582.1216
Z_LO, Z_HI = -RIB_T / 2.0, SPAN + RIB_T / 2.0
PIECES = "mesh/stage3_skin/stage3_skin_faces" + RS_SUFFIX + ".json"
if RS_SUFFIX:
    try:
        SOLID_VOL = float(open("mesh/stage3_skin/stage3_solid_volume" + RS_SUFFIX + ".txt").read())
    except OSError:
        raise SystemExit("missing mesh/stage3_skin/stage3_solid_volume" + RS_SUFFIX +
                         ".txt (run build_stage3_skin_geom.py with RIB_SPACING set)")


def die(msg):
    print(f"ERROR: {msg} -- STOP."); sys.exit(1)


# ---- airfoil helpers (copied unchanged from build_stage2/3 scripts) ----
def chord_at(z):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * (z / SPAN)


def naca0012_half_thickness(x_over_c, t=0.12):
    return 5 * t * (0.2969 * math.sqrt(x_over_c) - 0.1260 * x_over_c
                     - 0.3516 * x_over_c**2 + 0.2843 * x_over_c**3
                     - 0.1015 * x_over_c**4)


def airfoil_outline(chord, x_shift):
    xs = [chord * (1 - math.cos(math.pi * i / (2 * (N_AIRFOIL_PTS - 1))))
          for i in range(N_AIRFOIL_PTS)]
    upper = [(x - x_shift, chord * naca0012_half_thickness(x / chord)) for x in xs]
    lower = [(x - x_shift, -chord * naca0012_half_thickness(x / chord)) for x in reversed(xs[:-1])]
    return upper + lower[:-1]


def h_at(z):
    c = chord_at(z)
    up = airfoil_outline(c, SPAR_CHORD_FRAC * c)[:N_AIRFOIL_PTS]
    return float(np.interp(0.0, [p[0] for p in up], [p[1] for p in up]))


ZT = np.linspace(0.0, SPAN, 201)
HT = np.array([h_at(z) for z in ZT])


def edge_checks(cs, P, label):
    e = np.concatenate([np.sort(cs[:, [0, 1]], axis=1), np.sort(cs[:, [1, 2]], axis=1),
                        np.sort(cs[:, [0, 2]], axis=1)])
    ue, ecnt = np.unique(e, axis=0, return_counts=True)
    print(f"{label} edge-sharing histogram: {dict(zip(*np.unique(ecnt, return_counts=True)))}")
    bnd = ue[ecnt == 1]
    zb = P[bnd][:, :, 2]
    on_lo = np.all(np.abs(zb - Z_LO) < TOL, axis=1)
    on_hi = np.all(np.abs(zb - Z_HI) < TOL, axis=1)
    stray = int((~(on_lo | on_hi)).sum())
    print(f"{label} boundary edges: {len(bnd)} (root end {int(on_lo.sum())}, tip end "
          f"{int(on_hi.sum())}, elsewhere {stray})")
    return stray, int((ecnt > 2).sum())


def main():
    global TOTAL_LIFT_N
    if len(sys.argv) not in (3, 5):
        print("Usage: build_stage3_skin_inp.py <raw.npz> <skin_thickness_mm> "
              "[<lift_N_per_half_wing> <load_factor_n>]"); sys.exit(1)
    raw, t_skin = sys.argv[1], float(sys.argv[2])
    ng = 1.0
    if len(sys.argv) == 5:
        TOTAL_LIFT_N, ng = float(sys.argv[3]), float(sys.argv[4])
    base = os.path.splitext(os.path.basename(raw))[0]
    tag = (f"stage3_skin_{base[4:] if base.startswith('raw_') else base}_t{t_skin:g}"
           + (f"_n{ng:g}" if ng != 1.0 else "")).replace(".", "p")
    inp_out = f"simulation/{tag}.inp"
    os.makedirs("simulation", exist_ok=True)

    d = np.load(raw)
    P, tris, tri_idx, tets, inb = d["P"], d["tris"], d["tri_idx"], d["tets"], d["in_tet_bnd"]
    pieces = json.load(open(PIECES))
    print("=" * 70)
    print(f"STAGE 3 SKIN .inp: raw={raw}, skin thickness={t_skin} mm, "
          f"lift={TOTAL_LIFT_N} N, load factor n={ng}")
    print("=" * 70)

    # ---------- classify ----------
    p0, p1, p2 = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    cen = (p0 + p1 + p2) / 3.0
    nid = int(tri_idx.max()) + 1
    asum = np.bincount(tri_idx, weights=area, minlength=nid)
    csum = np.stack([np.bincount(tri_idx, weights=area * cen[:, k], minlength=nid)
                     for k in range(3)], axis=1)
    jc = np.array([[p["cx"], p["cy"], p["cz"]] for p in pieces])
    ja = np.array([p["area"] for p in pieces])
    jshared = np.array([p["shared"] for p in pieces])
    tree = cKDTree(jc)
    idx_cls = np.zeros(nid, dtype=int)
    nmatch = 0
    for i in np.nonzero(asum > 0)[0]:
        dist, j = tree.query(csum[i] / asum[i])
        if dist < 1e-3 and abs(asum[i] - ja[j]) / ja[j] < 1e-5:
            idx_cls[i] = 1 if jshared[j] else 2
            nmatch += 1
    tri_cls = idx_cls[tri_idx]
    print(f"Matched {nmatch}/{len(pieces)} skin pieces; tris: solid-nonskin {int((tri_cls == 0).sum())}, "
          f"shared skin {int((tri_cls == 1).sum())}, free skin {int((tri_cls == 2).sum())}")
    if nmatch != len(pieces):
        die("not all skin pieces matched")
    if int(((tri_cls == 1) & ~inb).sum()) or int(((tri_cls == 2) & inb).sum()):
        die("classification inconsistent with tet boundary")

    # ---------- weld ----------
    used = np.unique(np.concatenate([tets.ravel(), tris.ravel()]))
    pairs = cKDTree(P[used]).query_pairs(WELD, output_type="ndarray")
    n = len(used)
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    ncomp, comp = connected_components(g, directed=False)
    is_tet = np.isin(used, np.unique(tets))
    if np.bincount(comp[is_tet], minlength=ncomp).max() > 1:
        die("a weld group contains two solid nodes (duplicate solid nodes)")
    BIG = 10 ** 9
    key = np.where(is_tet, used, used + BIG).astype(np.int64)
    mn = np.full(ncomp, 2 ** 62, dtype=np.int64)
    np.minimum.at(mn, comp, key)
    cmap = np.zeros(int(used.max()) + 1, dtype=np.int64)
    cmap[used] = mn[comp] % BIG
    if not np.array_equal(cmap[tets], tets):
        die("weld changed a tet node ID")
    print(f"Weld: {n} node IDs -> {ncomp} nodes")

    # ---------- verify skin shell mesh ----------
    cs = cmap[tris[tri_cls > 0]]
    if ((cs[:, 0] == cs[:, 1]) | (cs[:, 1] == cs[:, 2]) | (cs[:, 0] == cs[:, 2])).any():
        die("degenerate skin triangle after weld")
    if len(np.unique(np.sort(cs, axis=1), axis=0)) != len(cs):
        die("duplicate skin triangles after weld")
    stray, nonmanifold = edge_checks(cs, P, "Skin")
    if stray or nonmanifold:
        die("skin mesh not watertight (stray boundary edges or edges shared by >2 triangles)")
    q0, q1, q2 = P[cs[:, 0]], P[cs[:, 1]], P[cs[:, 2]]
    ar = 0.5 * np.linalg.norm(np.cross(q1 - q0, q2 - q0), axis=1)
    skin_area = float(ar.sum())
    area_exp = float(ja.sum())
    print(f"Skin triangles: {len(cs)}; mesh area {skin_area:.4f} vs geometry {area_exp:.4f} mm^2 "
          f"({100*(skin_area-area_exp)/area_exp:+.5f}%)")
    if abs(skin_area - area_exp) / area_exp > 1e-4:
        die("skin mesh area disagrees with geometry")
    ea, eb, ec = (np.linalg.norm(q1 - q0, axis=1), np.linalg.norm(q2 - q1, axis=1),
                  np.linalg.norm(q0 - q2, axis=1))

    def angle(a, b, c):
        return np.degrees(np.arccos(np.clip((b * b + c * c - a * a) / (2 * b * c), -1, 1)))
    amin = np.minimum(np.minimum(angle(ea, eb, ec), angle(eb, ea, ec)), angle(ec, ea, eb))
    aspect = np.maximum(np.maximum(ea, eb), ec) ** 2 / ar / (4.0 / math.sqrt(3.0))
    print(f"Skin quality: min angle {amin.min():.2f} deg (1st pct {np.percentile(amin, 1):.2f}, "
          f"<10 deg: {int((amin < 10).sum())}); aspect metric (1=equilateral): "
          f"95th pct {np.percentile(aspect, 95):.2f}, max {aspect.max():.2f}")

    # ---------- outward winding ----------
    nrm = np.cross(q1 - q0, q2 - q0)
    vout = (q0 + q1 + q2) / 3.0
    vout[:, 2] = 0.0                      # from the spar axis point (0,0,z)
    dot = np.einsum("ij,ij->i", nrm, vout)
    flip = dot < 0
    cs = cs.copy()
    cs[flip] = cs[flip][:, [0, 2, 1]]
    print(f"Winding: {int(flip.sum())} of {len(cs)} skin triangles flipped to point outward")
    qq0, qq1, qq2 = P[cs[:, 0]], P[cs[:, 1]], P[cs[:, 2]]
    nn = np.cross(qq1 - qq0, qq2 - qq0)
    nn /= np.linalg.norm(nn, axis=1)[:, None]
    vv = (qq0 + qq1 + qq2) / 3.0
    vv[:, 2] = 0.0
    vv /= np.linalg.norm(vv, axis=1)[:, None]
    cosmin = float(np.einsum("ij,ij->i", nn, vv).min())
    print(f"  min cos(normal, outward vector) after flip: {cosmin:.4f} (must be > 0)")
    if cosmin <= 0:
        die("some skin normal still not outward")
    dir_e = np.concatenate([cs[:, [0, 1]], cs[:, [1, 2]], cs[:, [2, 0]]])
    M = int(cs.max()) + 1
    keys = dir_e[:, 0].astype(np.int64) * M + dir_e[:, 1]
    rev = dir_e[:, 1].astype(np.int64) * M + dir_e[:, 0]
    has_rev = np.isin(rev, keys)
    und = np.sort(dir_e, axis=1)
    _, inv, cnt = np.unique(und, axis=0, return_inverse=True, return_counts=True)
    ecount = cnt[inv.ravel()]
    n_inconsistent = int(((~has_rev) & (ecount == 2)).sum())
    print(f"  directed-edge consistency: edges whose two triangles disagree on direction: "
          f"{n_inconsistent} (must be 0)")
    if n_inconsistent:
        die("inconsistent skin winding")

    # ---------- renumber contiguously ----------
    V = np.unique(np.concatenate([tets.ravel(), cs.ravel()]))
    renum = np.zeros(int(P.shape[0]) + 1, dtype=np.int64)
    renum[V] = np.arange(1, len(V) + 1)
    Nv = len(V)
    C = np.vstack([np.zeros((1, 3)), P[V]])
    T = renum[tets][:, [0, 1, 3, 2]]
    S = renum[cs]
    nt, ns = len(T), len(S)
    pt = C[T]
    vols = np.einsum("ij,ij->i", np.cross(pt[:, 1] - pt[:, 0], pt[:, 2] - pt[:, 0]),
                     pt[:, 3] - pt[:, 0]) / 6.0
    print(f"\nTets {nt}: min volume {vols.min():.6f} (bad {int((vols <= 0).sum())}), "
          f"total {vols.sum():.4f} vs {SOLID_VOL} ({100*(vols.sum()-SOLID_VOL)/SOLID_VOL:+.4f}%)")
    if (vols <= 0).any() or abs(vols.sum() - SOLID_VOL) / SOLID_VOL > 0.005:
        die("tet volume check failed")

    # ---------- global midside table ----------
    tp = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    sp = [(0, 1), (1, 2), (2, 0)]
    te = np.stack([np.sort(np.stack([T[:, a], T[:, b]], axis=1), axis=1) for a, b in tp], axis=1)
    se = np.stack([np.sort(np.stack([S[:, a], S[:, b]], axis=1), axis=1) for a, b in sp], axis=1)
    alle = np.concatenate([te.reshape(-1, 2), se.reshape(-1, 2)])
    ue, inv = np.unique(alle, axis=0, return_inverse=True)
    mid_ids = Nv + 1 + inv.ravel()
    tmid = mid_ids[:nt * 6].reshape(nt, 6)
    smid = mid_ids[nt * 6:].reshape(ns, 3)
    allC = np.vstack([C, 0.5 * (C[ue[:, 0]] + C[ue[:, 1]])])
    N = len(allC) - 1
    conn_t = np.hstack([T, tmid])
    conn_s = np.hstack([S, smid])
    solid_nodes = np.unique(np.concatenate([T.ravel(), tmid.ravel()]))
    skin_nodes = np.unique(np.concatenate([S.ravel(), smid.ravel()]))
    n_shared_nodes = len(np.intersect1d(solid_nodes, skin_nodes))
    print(f"Nodes: {N} total ({Nv} vertex + {N - Nv} midside); solid {len(solid_nodes)}, "
          f"skin {len(skin_nodes)}, shared (solid & skin) {n_shared_nodes}")
    n_elem = nt + ns
    eq_est = 3 * N + 3 * len(skin_nodes)
    print(f"Elements: {nt} C3D10 + {ns} S6 = {n_elem} ({100*n_elem/ELEM_CEILING:.1f}% of ceiling)")
    print(f"Equation estimate: {eq_est} ({100*eq_est/EQ_CEILING:.1f}% of ceiling; "
          f"3/node + 3 extra rotational DOF per skin node -- an ESTIMATE)")
    if n_elem > ELEM_CEILING or eq_est > EQ_CEILING:
        die("size exceeds ceiling")

    # ---------- BC and load ----------
    z_all = allC[1:, 2]
    z_min = z_all.min()
    root_nodes = (np.nonzero(np.abs(z_all - z_min) < TOL)[0] + 1)
    n_root_skin = int(np.isin(root_nodes, skin_nodes).sum())
    n_root_solid = int(np.isin(root_nodes, solid_nodes).sum())
    print(f"\nRoot face z_min={z_min:.6f} (expected {Z_LO}); root nodes {len(root_nodes)} "
          f"(solid {n_root_solid}, skin {n_root_skin})")
    if abs(z_min - Z_LO) > TOL:
        die("root face not at expected z")

    x_all, y_all = allC[1:, 0], allC[1:, 1]
    hz = np.interp(z_all, ZT, HT)
    mask = (np.abs(x_all) <= CAP_W / 2.0 + TOL) & (np.abs(y_all - hz) < 0.15) & \
           (z_all >= 0.0) & (z_all <= SPAN)
    load_nodes = np.nonzero(mask)[0] + 1
    n_nonsolid = int((~np.isin(load_nodes, solid_nodes)).sum())
    print(f"Load nodes: {len(load_nodes)}; of which NOT solid nodes: {n_nonsolid} (must be 0)")
    if n_nonsolid or len(load_nodes) == 0:
        die("load-node selection contaminated or empty")
    SLAB = 5.0
    q0c = TOTAL_LIFT_N / (math.pi * SPAN / 4.0)

    def cum_q(zz):
        u = min(max(zz / SPAN, 0.0), 1.0)
        return q0c * (SPAN / 2.0) * (u * math.sqrt(1.0 - u * u) + math.asin(u))

    zl = allC[load_nodes, 2]
    n_slabs = int(round(SPAN / SLAB))
    sidx = np.minimum((zl / SLAB).astype(int), n_slabs - 1)
    counts = np.bincount(sidx, minlength=n_slabs)
    if (counts == 0).any():
        die("empty load slab")
    forces = np.zeros(len(load_nodes))
    for k in range(n_slabs):
        forces[sidx == k] = (cum_q((k + 1) * SLAB) - cum_q(k * SLAB)) / counts[k]
    M_app = float(np.sum(forces * zl))
    M_int = q0c * SPAN ** 2 / 3.0
    print(f"Total force {forces.sum():.6f} N (intended {TOTAL_LIFT_N}); load moment "
          f"{M_app:.4f} vs {M_int:.4f} N*mm ({100*(M_app-M_int)/M_int:+.4f}%)")
    if abs(M_app - M_int) / M_int > 0.005:
        die("load moment off by >0.5%")

    # ---------- write .inp ----------
    with open(inp_out, "w") as f:
        f.write("*HEADING\n")
        f.write(f"Project 03 Stage 3: spar v2 + ribs (C3D10) + skin (S6, t={t_skin} mm), {base}\n")
        f.write("** Skin nodes on the airfoil envelope, material outward (OFFSET=-0.5).\n")
        f.write(f"** LOADCASE LIFT_N={TOTAL_LIFT_N:.6f} NG={ng:.6f} (elliptical lift on spar cap face + self-weight x NG; root z_min fixed)\n")
        f.write("*NODE\n")
        f.write("\n".join(f"{i}, {allC[i,0]:.10f}, {allC[i,1]:.10f}, {allC[i,2]:.10f}"
                          for i in range(1, N + 1)) + "\n")
        f.write("*ELEMENT, TYPE=C3D10, ELSET=ESOLID\n")
        f.write("\n".join(f"{i+1}, " + ", ".join(map(str, conn_t[i])) for i in range(nt)) + "\n")
        f.write("*ELEMENT, TYPE=S6, ELSET=ESKIN\n")
        f.write("\n".join(f"{nt+i+1}, " + ", ".join(map(str, conn_s[i])) for i in range(ns)) + "\n")
        f.write("*NSET, NSET=NROOT\n")
        for i in range(0, len(root_nodes), 10):
            f.write(", ".join(map(str, root_nodes[i:i + 10])) + "\n")
        f.write("*MATERIAL, NAME=AL2024T3\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=ESOLID, MATERIAL=AL2024T3\n")
        f.write("*SHELL SECTION, ELSET=ESKIN, MATERIAL=AL2024T3, OFFSET=-0.5\n")
        f.write(f"{t_skin}\n")
        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNROOT, 1, 6, 0.0\n")
        f.write(f"*DLOAD\nESOLID, GRAV, {9810.0 * ng:.4f}, 0.0, -1.0, 0.0\nESKIN, GRAV, {9810.0 * ng:.4f}, 0.0, -1.0, 0.0\n")
        f.write("*CLOAD\n")
        f.write("\n".join(f"{nd}, 2, {fo:.8f}" for nd, fo in zip(load_nodes, forces)) + "\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")
    print(f"\nWrote: {inp_out}")

    # viewer file (final numbering; vertex nodes only)
    tri_all = renum[cmap[tris]]
    vout_npz = f"mesh/stage3_skin/{tag}_viewer.npz"
    np.savez(vout_npz, P=C, tris=tri_all, tri_cls=tri_cls, tets=T)
    print(f"Wrote: {vout_npz}")


if __name__ == "__main__":
    main()
