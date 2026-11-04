#!/usr/bin/python3
"""V4: linear through-thickness temperature gradient, free plate (three mid-plane pins).

Usage: /usr/bin/python3 scripts/check_v4.py JOB.dat MESH.inp T_BOTTOM T_TOP T_REF

References (thin-plate theory, linear): curvature kappa = alpha (T_top - T_bottom) / h in both directions,
free plate -> spherical shape, no stress. With the three pinned mid-plane corners at w = 0:
  w(x, y) = kappa/2 [ (L^2 + W^2)/4 - (x - L/2)^2 - (y - W/2)^2 ]
  centre rise w_c = kappa (L^2 + W^2) / 8, edge profile on x = 0 from the formula,
  mid-plane in-plane expansion u_x(L, 0) = alpha (T_mid - T_ref) L.
Tolerances (fixed in docs/specification.md before the run): 2 % on the deflection (centre rise and the
largest error along the x = 0 edge, both relative to the centre rise); von Mises stress below 0.1 % of
E alpha (T_top - T_bottom). The in-plane expansion is a supplementary comparison.
Writes results/v4_gradient.txt and exits with status 1 on failure.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio
from material import MATERIAL

TOL_W = 0.02
TOL_S = 0.001


def main():
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    dat, mesh = sys.argv[1], sys.argv[2]
    tb, tt, tref = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    par = ccxio.read_params(mesh)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mesh)
    blocks = ccxio.parse_dat(dat)
    found = sorted({(b["title"], b["set"]) for b in blocks})

    def last(prefix, name):
        sel = [b for b in blocks if b["title"].lower().startswith(prefix) and b["set"] == name]
        if not sel:
            raise SystemExit("no '%s' block for set %s; blocks found: %s" % (prefix, name, found))
        return sel[-1]

    def disp(name):
        return np.array(last("displacements", name)["rows"])

    E, a = MATERIAL["E"], MATERIAL["alpha"]
    dTg = tt - tb
    kappa = a * dTg / h
    wc_exact = kappa * (L ** 2 + W ** 2) / 8.0
    PA, PB, PC, PD, PCEN = (disp(n)[0] for n in ("P_A", "P_B", "P_C", "P_D", "P_CEN"))
    wc = PCEN[3]
    ec = 100.0 * abs(wc - wc_exact) / abs(wc_exact)

    edge = disp("FRAME_X0")
    y = np.array([nodes[int(r[0])][1] for r in edge])
    w_exact = kappa / 2.0 * ((L ** 2 + W ** 2) / 4.0 - (0.0 - L / 2.0) ** 2 - (y - W / 2.0) ** 2)
    ee = 100.0 * np.abs(edge[:, 3] - w_exact).max() / abs(wc_exact)

    s = np.array(last("stresses", "PANEL")["rows"])
    sxx, syy, szz, sxy, sxz, syz = (s[:, i] for i in range(2, 8))
    vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + sxz ** 2 + syz ** 2))
    ref = E * a * abs(dTg)
    es = 100.0 * vm.max() / ref
    ux_exact = a * ((tb + tt) / 2.0 - tref) * L
    eu = 100.0 * abs(PB[1] - ux_exact) / abs(ux_exact) if ux_exact != 0 else float("nan")
    rf = np.array(last("total force", "NALL")["rows"][0])

    ok = ec <= 100 * TOL_W and ee <= 100 * TOL_W and es <= 100 * TOL_S
    lines = ["V4 linear through-thickness gradient, free plate",
             "job %s, mesh %s" % (os.path.basename(dat), os.path.basename(mesh)),
             "T_bottom %.3f, T_top %.3f, T_ref %.3f degC; L = W = %g m, h = %g m" % (tb, tt, tref, L, h),
             "curvature kappa = alpha dT_g / h = %.6e 1/m; centre rise exact %.6e m (%.3f h)" % (kappa, wc_exact, wc_exact / h),
             "  corner heights (w): P_A %.3e, P_B %.3e, P_C %.3e (pinned at 0), P_D %.3e m" % (PA[3], PB[3], PC[3], PD[3]),
             "  centre rise: FE %.6e m, exact %.6e m, error %.3e %% (tolerance %.1f %%)" % (wc, wc_exact, ec, 100 * TOL_W),
             "  edge x = 0 profile (%d nodes): max |w - exact| = %.3e m = %.3e %% of the centre rise (tolerance %.1f %%)"
             % (len(edge), np.abs(edge[:, 3] - w_exact).max(), ee, 100 * TOL_W),
             "  max von Mises stress over %d integration points: %.4e Pa = %.3e %% of E alpha dT_g = %.4e Pa (tolerance %.1f %%)"
             % (len(vm), vm.max(), es, ref, 100 * TOL_S),
             "  supplementary: mid-plane expansion u_x(P_B) FE %.6e m, exact %.6e m, error %.3e %%" % (PB[1], ux_exact, eu),
             "  total reaction force (fx, fy, fz): %s N (expected 0)" % np.array2string(rf, precision=3),
             "RESULT: %s" % ("PASS" if ok else "FAIL")]
    text = "\n".join(lines) + "\n"
    sys.stdout.write(text)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "results"), exist_ok=True)
    with open(os.path.join(root, "results", "v4_gradient.txt"), "w") as f:
        f.write(text)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
