#!/usr/bin/env python3
"""Map the CFD wing-surface pressure load onto the structural nodes of the uCRM-9 coarse deck (D20, D21).

usage: map_cfd_loads.py <time> <case> <alpha_deg> <deck_in.inp> <out_dir> <tag> [--drop-fx] [--strip 1.0] [--pullback kin.csv]

Targets : the dataset-loaded nodes (y >= 2.94 m, same convention as the dataset) plus the free stub skin nodes (y < 2.94 m,
          shell elements with |n_z| >= 0.7, not constrained).
Baseline: each CFD face force goes to its nearest target node.
Correction (per spanwise strip): minimum-norm force-only correction that restores the exact strip moment about the strip centroid
          (total force unchanged).
Pullback: for a CFD wing that was deformed with the section kinematics of a structural result (scripts/deform_wing_stl.py),
          the face centres are first mapped back to the undeflected (jig) frame of the structural nodes (solving c0 + D(c0) = c);
          the force vectors are kept as computed (first-order, linear structure, loads applied in the jig frame).
Output  : <out_dir>/uCRM9_cfd_<tag>.inp (deck with the *CLOAD block replaced), <out_dir>/nodal_loads_<tag>.csv
"""
import os
import sys
import glob
import argparse
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import inspect_deck_loads as idl
import inspect_deck_skin2 as isk

PINF = 21662.7


def cfd_faces(time, case, alpha):
    files = glob.glob("%s/VTK/wing/wing_%s.vtk" % (case, time))
    if len(files) != 1:
        sys.exit("ABORT: wing VTK not found for %s %s" % (case, time))
    m = pv.read(files[0])
    m = m.compute_cell_sizes(length=False, area=True, volume=False)
    m = m.compute_normals(cell_normals=True, point_normals=False, consistent_normals=False,
                          auto_orient_normals=False, split_vertices=False)
    area = np.asarray(m.cell_data["Area"])
    nrm = np.asarray(m.cell_data["Normals"])
    p = np.asarray(m.cell_data["p"])
    C = np.asarray(m.cell_centers().points)
    Fv = ((p - PINF) * area)[:, None] * nrm
    a = np.radians(alpha)
    tot = Fv.sum(axis=0)
    lift = -tot[0] * np.sin(a) + tot[2] * np.cos(a)
    flipped = lift < 0
    if flipped:
        Fv = -Fv
    return C, Fv, flipped, files[0]


def pull_back_to_jig(C, csv_path, iters=10):
    """Points of the deflected wing -> undeflected frame: solve c0 + D(c0) = c by fixed-point iteration.
    D is the section-kinematics displacement of scripts/deform_wing_stl.py (zero-motion root point, constant extrapolation)."""
    K = np.genfromtxt(csv_path, delimiter=",", names=True)
    Y = np.concatenate([[0.0], K["y_c"]])
    UX = np.concatenate([[0.0], K["ux"]])
    UY = np.concatenate([[0.0], K["uy"]])
    UZ = np.concatenate([[0.0], K["uz"]])
    B = np.concatenate([[0.0], K["b_slope"]])
    XB = np.concatenate([[K["xbar"][0]], K["xbar"]])

    def disp(P):
        y = np.abs(P[:, 1])
        s = np.sign(P[:, 1])
        return np.column_stack([np.interp(y, Y, UX), s * np.interp(y, Y, UY),
                                np.interp(y, Y, UZ) + np.interp(y, Y, B) * (P[:, 0] - np.interp(y, Y, XB))])

    P = C.copy()
    for _ in range(iters):
        P = C - disp(P)
    resid = np.abs(P + disp(P) - C).max()
    moved = np.linalg.norm(C - P, axis=1)
    print("pullback to jig frame (%s): displacement removed mean %.3f m, max %.3f m, fixed-point residual %.1e m" % (
        csv_path, moved.mean(), moved.max(), resid))
    return P


def skew(r):
    return np.array([[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]])


def correct_strip(R, Mres):
    n = len(R)
    A = np.zeros((6, 3 * n))
    for i in range(n):
        A[0:3, 3 * i:3 * i + 3] = np.eye(3)
        A[3:6, 3 * i:3 * i + 3] = skew(R[i])
    b = np.concatenate([np.zeros(3), Mres])
    lam = np.linalg.lstsq(A @ A.T, b, rcond=1e-13)[0]
    return (A.T @ lam).reshape(n, 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("time")
    ap.add_argument("case")
    ap.add_argument("alpha", type=float)
    ap.add_argument("deck")
    ap.add_argument("outdir")
    ap.add_argument("tag")
    ap.add_argument("--drop-fx", action="store_true")
    ap.add_argument("--strip", type=float, default=1.0)
    ap.add_argument("--pullback", default=None)
    a = ap.parse_args()

    out_inp = os.path.join(a.outdir, "uCRM9_cfd_%s.inp" % a.tag)
    out_csv = os.path.join(a.outdir, "nodal_loads_%s.csv" % a.tag)
    if os.path.exists(out_inp):
        sys.exit("ABORT: %s already exists" % out_inp)

    nodes, block_of, cload, _, boundary = idl.parse(a.deck)
    loaded = set(i for i in cload if i in nodes)
    constrained = set(k for k, _, _ in boundary if isinstance(k, int))
    E, _ = isk.parse_shells(a.deck)
    X = np.array([[nodes[n] for n in e] for e in E])
    nv = np.cross(X[:, 2, :] - X[:, 0, :], X[:, 3, :] - X[:, 1, :])
    nl = np.linalg.norm(nv, axis=1)
    nz = np.abs(nv[:, 2]) / np.where(nl > 0, nl, 1.0)
    skin_nodes = set(int(n) for n in E[nz >= 0.7].ravel())
    stub = set(n for n in skin_nodes if nodes[n][1] < 2.94 and n not in constrained)
    targets = sorted(loaded | stub)
    T = np.array([nodes[n] for n in targets])
    print("deck in:", a.deck)
    print("targets: dataset-loaded %d + free stub skin %d = %d nodes (overlap %d)" % (len(loaded), len(stub), len(targets), len(loaded & stub)))

    C, Fv, flipped, wfile = cfd_faces(a.time, a.case, a.alpha)
    print("CFD wing patch:", wfile, " faces", len(C), " normals flipped:", flipped)
    if a.pullback:
        C = pull_back_to_jig(C, a.pullback)
    if a.drop_fx:
        Fv[:, 0] = 0.0
        print("--drop-fx: CFD x-force component set to zero before mapping")
    F_cfd = Fv.sum(axis=0)
    M_cfd_o = np.cross(C, Fv).sum(axis=0)
    print("CFD force  Fx %.1f  Fy %.1f  Fz %.1f N" % tuple(F_cfd))

    dy = a.strip
    s_node = np.floor(T[:, 1] / dy).astype(int)
    tree = cKDTree(T)
    dist, j = tree.query(C)
    s_face = s_node[j]
    print("nearest-node distance of CFD faces: median %.3f  90%% %.3f  max %.3f m" % (np.median(dist), np.percentile(dist, 90), dist.max()))
    f0 = np.zeros((len(targets), 3))
    np.add.at(f0, j, Fv)
    f = f0.copy()

    print("\nper-strip moment correction (strip width %.2f m)" % dy)
    print("  y0   nodes  faces   Fz[kN]  |Mres|before[kNm]  Mres after[Nm]  force err[N]  max|df|[N]  max|f0|[N]  Fz reversed[%]")
    total_df = 0.0
    total_f0 = 0.0
    for s in sorted(set(s_node.tolist())):
        ni = np.where(s_node == s)[0]
        fi = np.where(s_face == s)[0]
        if len(fi) == 0:
            continue
        c = T[ni].mean(axis=0)
        R = T[ni] - c
        Fs = Fv[fi].sum(axis=0)
        Mcfd = np.cross(C[fi] - c, Fv[fi]).sum(axis=0)
        Mbase = np.cross(R, f0[ni]).sum(axis=0)
        Mres = Mcfd - Mbase
        df = correct_strip(R, Mres)
        f[ni] = f0[ni] + df
        Mafter = np.cross(R, f[ni]).sum(axis=0)
        ferr = np.linalg.norm(f[ni].sum(axis=0) - Fs)
        rev = 100.0 * np.mean(np.sign(f[ni, 2]) == -np.sign(Fs[2])) if abs(Fs[2]) > 0 else 0.0
        total_df += np.abs(df).sum()
        total_f0 += np.abs(f0[ni]).sum()
        print("%5.1f %6d %6d %9.2f %14.2f %16.3e %12.3e %10.1f %10.1f %10.1f" % (
            s * dy, len(ni), len(fi), Fs[2] / 1e3, np.linalg.norm(Mres) / 1e3,
            np.linalg.norm(Mcfd - Mafter), ferr, np.abs(df).max(), np.abs(f0[ni]).max(), rev))

    F_map = f.sum(axis=0)
    M_map_o = np.cross(T, f).sum(axis=0)
    print("\n=== totals ===")
    print("CFD    force  Fx %.3f  Fy %.3f  Fz %.3f N" % tuple(F_cfd))
    print("mapped force  Fx %.3f  Fy %.3f  Fz %.3f N" % tuple(F_map))
    print("relative force difference  %.3e" % (np.linalg.norm(F_map - F_cfd) / np.linalg.norm(F_cfd)))
    print("CFD    moment about origin  Mx %.6e  My %.6e  Mz %.6e N m" % tuple(M_cfd_o))
    print("mapped moment about origin  Mx %.6e  My %.6e  Mz %.6e N m" % tuple(M_map_o))
    print("relative moment difference %.3e" % (np.linalg.norm(M_map_o - M_cfd_o) / np.linalg.norm(M_cfd_o)))
    print("sum|correction| / sum|baseline| = %.3f" % (total_df / total_f0))
    print("largest nodal force magnitude: baseline %.1f N, final %.1f N" % (np.linalg.norm(f0, axis=1).max(), np.linalg.norm(f, axis=1).max()))

    print("\n=== mapped vs dataset load, y bins (kN) ===")
    DF = np.array([cload[i][:3] for i in sorted(loaded)])
    DX = np.array([nodes[i] for i in sorted(loaded)])
    print("   y0    y1 | mapped Fx   Fy   Fz | dataset Fx   Fy   Fz | Fz ratio")
    for lo, hi in zip(np.arange(0.0, 30.0, 3.0), np.arange(3.0, 33.0, 3.0)):
        mm = (T[:, 1] >= lo) & (T[:, 1] < hi)
        md = (DX[:, 1] >= lo) & (DX[:, 1] < hi)
        a1 = f[mm].sum(axis=0) / 1e3
        a2 = DF[md].sum(axis=0) / 1e3
        ratio = a1[2] / a2[2] if abs(a2[2]) > 1e-9 else float("nan")
        print("%5.1f %5.1f | %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f | %8.3f" % (lo, hi, a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], ratio))

    os.makedirs(a.outdir, exist_ok=True)
    with open(out_csv, "w") as g:
        g.write("node,x,y,z,fx,fy,fz\n")
        for n, p, v in zip(targets, T, f):
            g.write("%d,%.6f,%.6f,%.6f,%.9e,%.9e,%.9e\n" % (n, p[0], p[1], p[2], v[0], v[1], v[2]))
    with open(a.deck) as g:
        lines = g.readlines()
    i0 = None
    i1 = None
    for k, ln in enumerate(lines):
        s = ln.strip().upper()
        if i0 is None and s.startswith("*CLOAD"):
            i0 = k
        elif i0 is not None and ln.startswith("*") and not ln.startswith("**"):
            i1 = k
            break
    if i0 is None or i1 is None:
        sys.exit("ABORT: could not locate the *CLOAD block in the deck")
    new = ["** aero load replaced by the CFD-mapped load (scripts/map_cfd_loads.py, tag %s, time %s, alpha %.3f deg%s%s)\n" % (
        a.tag, a.time, a.alpha, ", Fx dropped" if a.drop_fx else "", ", pulled back to jig frame" if a.pullback else ""), "*CLOAD\n"]
    nl_ = 0
    for n, v in zip(targets, f):
        for d in range(3):
            if v[d] != 0.0:
                new.append("%d, %d, %.9e\n" % (n, d + 1, v[d]))
                nl_ += 1
    with open(out_inp, "w") as g:
        g.writelines(lines[:i0] + new + lines[i1:])
    print("\nreplaced old *CLOAD block (%d lines) by %d lines" % (i1 - i0 - 1, nl_))
    print("wrote", out_inp)
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
