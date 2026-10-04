#!/usr/bin/env python3
"""
Indirect equilibrium / load-application check for the skin model, using three solves of the
same mesh (lift only, lift + solid gravity, lift + all gravity).

For each variant: spar-axis tip Uy (SOLID nodes only, so no missing-shell-node zero rows)
vs tapered beam (spar v2 + fully effective skin ring "B_hi" and spar only "B_lo") loaded with
the SAME load set. Then the EFFECT of each added load (solid gravity, skin gravity) in FEA vs
beam. Reaction forces are NOT used (not interpretable for expanded shells).

Usage (system Python):
  /usr/bin/python3 scripts/check_skin_load_variants.py simulation/stage3_skin_s2p5_k8_t1
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import (parse_inp, read_frd, ring_props, ZTAB, A_of_z, I_of_z,
                                 RIBW, E_MOD, RHO, G, TOTAL_LIFT_N, SPAN)

VARIANTS = [
    ("lift only            ", "_nograv", dict(spar=False, ribs=False, skin=False)),
    ("lift + solid gravity ", "_solidgrav", dict(spar=True, ribs=True, skin=False)),
    ("full (all gravity)   ", "", dict(spar=True, ribs=True, skin=True)),
]


def beam_tip(t, w, skin_stiff):
    ptab = np.array([ring_props(z, t) for z in ZTAB])
    perim_tab, Iring_tab = ptab[:, 0], ptab[:, 1]
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = TOTAL_LIFT_N / (np.pi * SPAN / 4.0)
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN) ** 2, 0, None))
    if w["spar"]:
        q = q - RHO * G * A_of_z(z)
    if w["skin"]:
        q = q - RHO * G * np.interp(z, ZTAB, perim_tab) * t
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i + 1] + 0.5 * (q[i] + q[i + 1]) * dz
    if w["ribs"]:
        for zr, _, _, P in RIBW:
            V[z < zr - 1e-9] -= P
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i + 1] + 0.5 * (V[i] + V[i + 1]) * dz
    EI = E_MOD * (I_of_z(z) + (np.interp(z, ZTAB, Iring_tab) if skin_stiff else 0.0))
    curv = M / EI
    th = np.zeros(n); v = np.zeros(n)
    for i in range(1, n):
        th[i] = th[i - 1] + 0.5 * (curv[i - 1] + curv[i]) * dz
    for i in range(1, n):
        v[i] = v[i - 1] + 0.5 * (th[i - 1] + th[i]) * dz
    return v[-1]


def main():
    base = sys.argv[1]
    res = []
    print(f"{'variant':22s} {'FEA Uy [mm]':>12s} {'B_hi [mm]':>11s} {'FEA/B_hi-1':>11s} {'B_lo [mm]':>11s}")
    for label, suf, w in VARIANTS:
        job = base + suf
        xyz, solid, shell, root, cload, t_skin = parse_inp(job + ".inp")
        coords = np.vstack([np.zeros((1, 3)), xyz])
        U = read_frd(job + ".frd", "DISP", 3)
        solid_ids = np.unique(solid)
        z = coords[solid_ids, 2]
        x = coords[solid_ids, 0]
        sel = solid_ids[(np.abs(z - SPAN) < 2.0) & (np.abs(x) < 4.5)]
        uy = U[sel, 1]
        if uy.min() == 0.0:
            print("WARNING: zero row among selected nodes")
        fea = uy.mean()
        bhi = beam_tip(t_skin, w, True)
        blo = beam_tip(t_skin, w, False)
        res.append((fea, bhi, blo))
        print(f"{label:22s} {fea:12.6f} {bhi:11.6f} {100*(fea/bhi-1):+10.2f}% {blo:11.6f}   "
              f"(n={len(sel)}, spread {uy.min():.5f}..{uy.max():.5f})")
    (f0, h0, _), (f1, h1, _), (f2, h2, _) = res
    print("\nEffect of each added load on the tip deflection (FEA vs beam B_hi):")
    print(f"  solid gravity: FEA {f1 - f0:+.6f} mm, beam {h1 - h0:+.6f} mm, FEA/beam = {(f1 - f0) / (h1 - h0):.3f}")
    print(f"  skin gravity : FEA {f2 - f1:+.6f} mm, beam {h2 - h1:+.6f} mm, FEA/beam = {(f2 - f1) / (h2 - h1):.3f}")


if __name__ == "__main__":
    main()
