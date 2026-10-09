#!/usr/bin/env python3
"""
Fully effective skin-ring beam reference (as extract_stage3_skin.beam_tip) with an optional rib-section term:
over |z - z_rib| < RIB_T/2 the section is the SOLID airfoil (inside the envelope) plus the skin ring (outside the
envelope): I = I_solid(z) + I_ring(z) instead of I_spar(z) + I_ring(z). The no-rib variant is validated against
extract_stage3_skin.beam_tip. Run with the SAME RIB_SPACING env as the job (system python).
Usage: RIB_SPACING=100 /usr/bin/python3 scripts/beam_with_rib_inertia.py [t_skin_mm]
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_stage3_skin as X

b = X.b
LIFT, NG = 372.78, 3.8


def solid_I(z):
    pts = b.outline_at(z)
    n = len(pts)
    s = 0.0
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        s += (x0 * y1 - x1 * y0) * (y0 * y0 + y0 * y1 + y1 * y1)
    return abs(s) / 12.0


def beam(t, lift, ng, rib_inertia, ptab, isolid_tab):
    perim_tab, Iring_tab = ptab[:, 0], ptab[:, 1]
    n = 20000
    z = np.linspace(0, X.SPAN, n)
    dz = z[1] - z[0]
    q0 = lift / (np.pi * X.SPAN / 4.0)
    A_skin = np.interp(z, X.ZTAB, perim_tab) * t
    q = q0 * np.sqrt(np.clip(1.0 - (z / X.SPAN) ** 2, 0, None)) - X.RHO * X.G * ng * (X.A_of_z(z) + A_skin)
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i + 1] + 0.5 * (q[i] + q[i + 1]) * dz
    for zr, _, _, P in X.RIBW:
        V[z < zr - 1e-9] -= P * ng
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i + 1] + 0.5 * (V[i] + V[i + 1]) * dz
    ring = np.interp(z, X.ZTAB, Iring_tab)
    Itot = X.I_of_z(z) + ring
    if rib_inertia:
        isol = np.interp(z, X.ZTAB, isolid_tab)
        for zr in b.RIB_STATIONS:
            m = np.abs(z - zr) < b.RIB_T / 2.0
            Itot[m] = isol[m] + ring[m]
    curv = M / (X.E_MOD * Itot)
    th = np.zeros(n)
    v = np.zeros(n)
    for i in range(1, n):
        th[i] = th[i - 1] + 0.5 * (curv[i - 1] + curv[i]) * dz
    for i in range(1, n):
        v[i] = v[i - 1] + 0.5 * (th[i - 1] + th[i]) * dz
    return v[-1]


def main():
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    ptab = np.array([X.ring_props(z, t) for z in X.ZTAB])
    isolid_tab = np.array([solid_I(z) for z in X.ZTAB])
    print(f"RIB_SPACING = {b.RIB_SPACING:g} mm, {len(b.RIB_STATIONS)} rib stations; t_skin = {t} mm")
    print(f"I_solid at root {isolid_tab[0]:.1f}, mid {isolid_tab[100]:.1f}, tip {isolid_tab[-1]:.1f} mm^4; "
          f"I_ring at root {ptab[0, 1]:.1f}, mid {ptab[100, 1]:.1f}, tip {ptab[-1, 1]:.1f}; "
          f"I_spar at root {float(X.I_of_z(np.array([0.0]))[0]):.1f}")
    ref = X.beam_tip(t, True, LIFT, NG)
    mine = beam(t, LIFT, NG, False, ptab, isolid_tab)
    print(f"validation: extract_stage3_skin.beam_tip {ref:.6f} mm, this script (no rib inertia) {mine:.6f} mm, "
          f"diff {mine - ref:+.2e}")
    for lab, lift, ng in (("lift only", LIFT, 0.0), ("gravity only", 0.0, NG), ("full", LIFT, NG)):
        a = beam(t, lift, ng, False, ptab, isolid_tab)
        c = beam(t, lift, ng, True, ptab, isolid_tab)
        print(f"{lab:13s}: no rib inertia {a:10.6f} mm | with rib-section inertia {c:10.6f} mm | change {100*(c/a-1):+.3f}%")


if __name__ == "__main__":
    main()
