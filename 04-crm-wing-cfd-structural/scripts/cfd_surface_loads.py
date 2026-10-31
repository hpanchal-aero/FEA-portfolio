#!/usr/bin/env python3
"""CFD wing-surface pressure loads vs the dataset aero loads of the structural deck (read-only).
usage: cfd_surface_loads.py <time> <case> <alpha_deg> <deck.inp>
Force on the wing = sum (p - pInf) * Sf, Sf = face area vector pointing out of the fluid domain.
"""
import os
import sys
import glob
import numpy as np
import pyvista as pv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import inspect_deck_loads as idl

PINF = 21662.7
Q = 10955.9
AREF = 191.870


def main(time, case, alpha, deck):
    files = glob.glob("%s/VTK/wing/wing_%s.vtk" % (case, time))
    if len(files) != 1:
        print("ABORT: wing VTK not found for", case, time)
        sys.exit(1)
    m = pv.read(files[0])
    print("wing patch:", files[0], " cells", m.n_cells, " cell arrays", list(m.cell_data.keys()))
    if "p" not in m.cell_data:
        print("ABORT: no cell array p on the wing patch")
        sys.exit(1)
    m = m.compute_cell_sizes(length=False, area=True, volume=False)
    m = m.compute_normals(cell_normals=True, point_normals=False, consistent_normals=False,
                          auto_orient_normals=False, split_vertices=False)
    area = np.asarray(m.cell_data["Area"])
    nrm = np.asarray(m.cell_data["Normals"])
    p = np.asarray(m.cell_data["p"])
    C = np.asarray(m.cell_centers().points)
    Fv = ((p - PINF) * area)[:, None] * nrm
    tot = Fv.sum(axis=0)
    a = np.radians(alpha)
    lift = -tot[0] * np.sin(a) + tot[2] * np.cos(a)
    drag = tot[0] * np.cos(a) + tot[2] * np.sin(a)
    flipped = False
    if lift < 0:
        Fv = -Fv
        tot = -tot
        lift = -lift
        drag = -drag
        flipped = True
    print("normals flipped to give positive lift:", flipped)
    print("\n=== verification against the solver (alpha %.2f deg) ===" % alpha)
    print("total force  Fx %.1f  Fy %.1f  Fz %.1f N" % tuple(tot))
    print("Cl from surface pressure %.5f   Cd %.5f   (solver forceCoeffs at 1500 is in the log)" % (lift / (Q * AREF), drag / (Q * AREF)))
    print("lift %.1f N  drag %.1f N" % (lift, drag))

    nodes, block_of, cload, _, _ = idl.parse(deck)
    ids = [i for i in sorted(cload) if i in nodes]
    DF = np.array([cload[i][:3] for i in ids])
    DX = np.array([nodes[i] for i in ids])
    dtot = DF.sum(axis=0)

    edges = np.arange(0.0, 33.0, 3.0)
    print("\n=== spanwise bins: CFD vs dataset (kN) ===")
    print("   y0    y1 |  CFD Fx   CFD Fy   CFD Fz |  DS Fx    DS Fy    DS Fz | CFD/DS Fz")
    for lo, hi in zip(edges[:-1], edges[1:]):
        mc = (C[:, 1] >= lo) & (C[:, 1] < hi)
        md = (DX[:, 1] >= lo) & (DX[:, 1] < hi)
        fc = Fv[mc].sum(axis=0) / 1000.0
        fd = DF[md].sum(axis=0) / 1000.0
        ratio = fc[2] / fd[2] if abs(fd[2]) > 1e-6 else float("nan")
        print("%5.1f %5.1f | %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f | %8.3f" % (lo, hi, fc[0], fc[1], fc[2], fd[0], fd[1], fd[2], ratio))

    # footprint classification: dataset-loaded chord range at each span station
    yb = np.arange(3.0, 30.0, 1.0)
    cen, xlo, xhi = [], [], []
    for b in yb:
        md = (DX[:, 1] >= b) & (DX[:, 1] < b + 1.0)
        if md.sum() > 0:
            cen.append(b + 0.5)
            xlo.append(DX[md, 0].min())
            xhi.append(DX[md, 0].max())
    cen, xlo, xhi = np.array(cen), np.array(xlo), np.array(xhi)
    YD = 2.94
    inb = C[:, 1] < YD
    lo_x = np.interp(C[:, 1], cen, xlo)
    hi_x = np.interp(C[:, 1], cen, xhi)
    ahead = (~inb) & (C[:, 0] < lo_x)
    behind = (~inb) & (C[:, 0] > hi_x)
    within = (~inb) & (~ahead) & (~behind)
    print("\n=== where the CFD load acts relative to the dataset load footprint ===")
    print("(inboard = y < %.2f m; ahead/behind = outside the loaded-node chord range at that y)" % YD)
    for name, mk in (("inboard y<2.94", inb), ("ahead of footprint", ahead), ("within footprint", within), ("behind footprint", behind)):
        f = Fv[mk].sum(axis=0)
        print("%-20s faces %7d  area %8.2f m2  Fx %9.2f  Fy %9.2f  Fz %10.2f kN   Fz share %5.1f %%" % (
            name, mk.sum(), area[mk].sum(), f[0] / 1e3, f[1] / 1e3, f[2] / 1e3, 100.0 * f[2] / tot[2]))

    outb = ~inb
    fo = Fv[outb].sum(axis=0)
    mo = np.cross(C[outb], Fv[outb]).sum(axis=0)
    md_ = np.cross(DX, DF).sum(axis=0)
    print("\n=== y >= 2.94 m only, CFD vs dataset ===")
    print("CFD     force  Fx %.1f  Fy %.1f  Fz %.1f N" % tuple(fo))
    print("dataset force  Fx %.1f  Fy %.1f  Fz %.1f N" % tuple(dtot))
    print("Fz ratio CFD/dataset: %.4f" % (fo[2] / dtot[2]))
    print("CFD     moment about origin  Mx %.4e  My %.4e  Mz %.4e N m" % tuple(mo))
    print("dataset moment about origin  Mx %.4e  My %.4e  Mz %.4e N m" % tuple(md_))
    print("CFD     centre of pressure x %.3f  y %.3f m" % ((C[outb, 0] * Fv[outb, 2]).sum() / fo[2], (C[outb, 1] * Fv[outb, 2]).sum() / fo[2]))
    print("dataset centre of pressure x %.3f  y %.3f m" % ((DX[:, 0] * DF[:, 2]).sum() / dtot[2], (DX[:, 1] * DF[:, 2]).sum() / dtot[2]))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4])
