#!/usr/bin/python3
"""V2: uniform temperature rise, minimally constrained (free expansion).

Usage: /usr/bin/python3 scripts/check_v2.py JOB.dat MESH.inp T_FINAL T_REF

References: u_x(L) = alpha dT L, u_y(W) = alpha dT W, centre = half of that,
u_z(z) = alpha dT (z - h/2) on the centre column, stress = 0, total reaction force = 0.
Tolerances (fixed in docs/specification.md before the run): 1 % on displacements
(relative to the expected value, or to alpha dT h/2 for the thickness expansion),
von Mises stress below 0.1 % of E alpha dT.
Writes results/v2_free_expansion.txt and exits with status 1 on failure.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio
from material import MATERIAL

TOL_U = 0.01
TOL_S = 0.001


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    dat, mesh = sys.argv[1], sys.argv[2]
    t_final, t_ref = float(sys.argv[3]), float(sys.argv[4])
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

    E, nu, a = MATERIAL["E"], MATERIAL["nu"], MATERIAL["alpha"]
    dT = t_final - t_ref
    ex, ey = a * dT * L, a * dT * W
    PA, PB, PC, PD, PCEN = (disp(n)[0] for n in ("P_A", "P_B", "P_C", "P_D", "P_CEN"))

    def rel(u, ue):
        return 100.0 * abs(u - ue) / abs(ue)

    checks = [
        ("P_B  ux", PB[1], ex), ("P_C  uy", PC[2], ey),
        ("P_D  ux", PD[1], ex), ("P_D  uy", PD[2], ey),
        ("P_CEN ux", PCEN[1], 0.5 * ex), ("P_CEN uy", PCEN[2], 0.5 * ey),
    ]
    lines = ["V2 uniform temperature rise, free expansion",
             "job %s, mesh %s" % (os.path.basename(dat), os.path.basename(mesh)),
             "dT = %.1f K, E = %.4g Pa, alpha = %.3g 1/K, L = %g, W = %g, h = %g m" % (dT, E, a, L, W, h),
             "P_A displacement (fixed point): %s" % np.array2string(PA[1:], precision=4),
             "  quantity      FE [m]          exact [m]      error [%]"]
    ok = True
    for name, u, ue in checks:
        e = rel(u, ue)
        ok = ok and e <= 100 * TOL_U
        lines.append("  %-9s %13.6e   %13.6e   %9.3e" % (name, u, ue, e))
    col = np.array(last("displacements", "COL_CEN")["rows"])
    z = np.array([nodes[int(r[0])][2] for r in col])
    uz_exact = a * dT * (z - h / 2)
    scale = a * dT * h / 2
    ez = 100.0 * np.abs(col[:, 3] - uz_exact).max() / scale
    ok = ok and ez <= 100 * TOL_U
    lines.append("  thickness expansion on the centre column: max |uz - exact| = %.3e m = %.3e %% of alpha dT h/2 (exact face value %.3e m)"
                 % (np.abs(col[:, 3] - uz_exact).max(), ez, scale))
    s = np.array(last("stresses", "PANEL")["rows"])
    sxx, syy, szz, sxy, sxz, syz = (s[:, i] for i in range(2, 8))
    vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + sxz ** 2 + syz ** 2))
    ref = E * a * abs(dT)
    es = 100.0 * vm.max() / ref
    ok = ok and es <= 100 * TOL_S
    lines.append("  max von Mises stress over %d integration points: %.4e Pa = %.3e %% of E alpha dT = %.4e Pa (tolerance %.1f %%)"
                 % (len(vm), vm.max(), es, ref, 100 * TOL_S))
    rf = np.array(last("total force", "NALL")["rows"][0])
    lines.append("  total reaction force (fx, fy, fz): %s N (expected 0)" % np.array2string(rf, precision=3))
    lines.append("tolerances: displacements %.1f %%, stress %.1f %% of E alpha dT" % (100 * TOL_U, 100 * TOL_S))
    lines.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    text = "\n".join(lines) + "\n"
    sys.stdout.write(text)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "results"), exist_ok=True)
    with open(os.path.join(root, "results", "v2_free_expansion.txt"), "w") as f:
        f.write(text)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
