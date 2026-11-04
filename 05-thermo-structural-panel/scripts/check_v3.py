#!/usr/bin/python3
"""V3: uniform temperature rise, frame support (in-plane restrained on the mid-plane), linear.

Usage: /usr/bin/python3 scripts/check_v3.py JOB.dat MESH.inp T_FINAL T_REF

References (plane stress in the interior, central region only):
  sigma_x = sigma_y = -E alpha dT / (1 - nu), sigma_z = 0, shear = 0,
  u_z(z) = alpha dT (1 + nu)/(1 - nu) (z - h/2) on the centre column.
Central region: elements whose centroid is within R = 0.1 m of the plate centre in x and y.
Tolerance (fixed in docs/specification.md before the run): 1 %.
The edge reaction force is reported as a supplementary comparison only.
Writes results/v3_restrained.txt and exits with status 1 on failure.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio
from material import MATERIAL

TOL = 0.01
R = 0.1


def main():
    if len(sys.argv) not in (5, 6):
        raise SystemExit(__doc__)
    dat, mesh = sys.argv[1], sys.argv[2]
    t_final, t_ref = float(sys.argv[3]), float(sys.argv[4])
    par = ccxio.read_params(mesh)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mesh)
    eids, conn = ccxio.read_elements(mesh)
    blocks = ccxio.parse_dat(dat)
    found = sorted({(b["title"], b["set"]) for b in blocks})

    def last(prefix, name):
        sel = [b for b in blocks if b["title"].lower().startswith(prefix) and b["set"] == name]
        if not sel:
            raise SystemExit("no '%s' block for set %s; blocks found: %s" % (prefix, name, found))
        return sel[-1]

    E, nu, a = MATERIAL["E"], MATERIAL["nu"], MATERIAL["alpha"]
    dT = t_final - t_ref
    dT_cr = 4 * math.pi ** 2 * h ** 2 / (12 * (1 + nu) * a * L * W)
    sig_ref = -E * a * dT / (1 - nu)

    # central-region element mask from the element centroids (corner nodes only)
    cen = np.array([np.mean([nodes[int(n)] for n in c[:8]], axis=0) for c in conn])
    inside = (np.abs(cen[:, 0] - L / 2) < R) & (np.abs(cen[:, 1] - W / 2) < R)
    sel_ids = set(int(i) for i in eids[inside])
    s = np.array(last("stresses", "PANEL")["rows"])
    keep = np.array([int(e) in sel_ids for e in s[:, 0]])
    sc = s[keep]
    sxx, syy, szz, sxy, sxz, syz = (sc[:, i] for i in range(2, 8))

    lines = ["V3 uniform temperature rise, linear (supports as defined in the job deck)",
             "job %s, mesh %s" % (os.path.basename(dat), os.path.basename(mesh)),
             "dT = %.2f K; closed-form buckling dT_cr = %.2f K (dT/dT_cr = %.2f); h = %g m, L = W = %g m" % (dT, dT_cr, dT / dT_cr, h, L),
             "central region: %d elements, %d integration points" % (inside.sum(), keep.sum()),
             "reference sigma_x = sigma_y = -E alpha dT/(1-nu) = %.6e Pa" % sig_ref,
             "  component   mean [Pa]      min [Pa]       max [Pa]     max |dev| from reference [%]"]
    ok = True
    for name, v in (("sigma_x", sxx), ("sigma_y", syy)):
        dev = 100.0 * np.abs(v - sig_ref).max() / abs(sig_ref)
        ok = ok and dev <= 100 * TOL
        lines.append("  %-9s %13.6e %13.6e %13.6e %10.3e" % (name, v.mean(), v.min(), v.max(), dev))
    for name, v in (("sigma_z", szz), ("sigma_xy", sxy), ("sigma_xz", sxz), ("sigma_yz", syz)):
        dev = 100.0 * np.abs(v).max() / abs(sig_ref)
        ok = ok and dev <= 100 * TOL
        lines.append("  %-9s %13.6e %13.6e %13.6e %10.3e  (|max| as %% of |reference|)" % (name, v.mean(), v.min(), v.max(), dev))

    col = np.array(last("displacements", "COL_CEN")["rows"])
    z = np.array([nodes[int(r[0])][2] for r in col])
    uz_exact = a * dT * (1 + nu) / (1 - nu) * (z - h / 2)
    scale = a * dT * (1 + nu) / (1 - nu) * h / 2
    ez = 100.0 * np.abs(col[:, 3] - uz_exact).max() / scale
    ok = ok and ez <= 100 * TOL
    lines.append("thickness expansion on the centre column: max |uz - exact| = %.3e m = %.3e %% of the exact face value %.4e m"
                 % (np.abs(col[:, 3] - uz_exact).max(), ez, scale))
    lines.append("  uz at the faces: bottom %.4e m, top %.4e m" % (col[np.argmin(z), 3], col[np.argmax(z), 3]))

    f_exact = abs(sig_ref) * h * W
    lines.append("supplementary: edge reaction force, exact (interior stress x h x W) %.2f N" % f_exact)
    prefix = "FACE" if any(b["set"] == "FACE_X0" for b in blocks) else "FRAME"
    for name in (prefix + "_X0", prefix + "_X1", prefix + "_Y0", prefix + "_Y1"):
        rf = last("total force", name)["rows"][0]
        comp = 0 if name.endswith(("_X0", "_X1")) else 1
        lines.append("  %-8s total reaction (fx, fy, fz) = %s N, normal component %.2f N = %.2f %% of exact"
                     % (name, np.array2string(np.array(rf), precision=3), rf[comp], 100.0 * abs(rf[comp]) / f_exact))
    rf = last("total force", "NALL")["rows"][0]
    lines.append("  all nodes total reaction (fx, fy, fz) = %s N (expected 0)" % np.array2string(np.array(rf), precision=3))
    lines.append("tolerance %.1f %% on central-region stresses and thickness expansion" % (100 * TOL))
    lines.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    text = "\n".join(lines) + "\n"
    sys.stdout.write(text)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "results"), exist_ok=True)
    with open(os.path.join(root, "results", sys.argv[5] if len(sys.argv) == 6 else "v3_restrained.txt"), "w") as f:
        f.write(text)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
