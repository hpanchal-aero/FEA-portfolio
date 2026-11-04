#!/usr/bin/python3
"""V7 references: post-buckling deflection of a simply supported square plate (normal in-plane edge displacement restrained,
tangential free, zero edge shear) under a uniform temperature rise, first-mode imperfection.

Usage: /usr/bin/python3 scripts/ref_postbuckling.py [original|corrected]      (default: corrected)
  original  : Amendment 2 reference (single mode, uniform membrane force, Berger-type); writes results/v7_reference.txt
  corrected : Amendment 3 reference (single mode with the full von Karman membrane stress); writes results/v7_reference_corrected.txt

Common definitions: w = W sin(pi x/L) sin(pi y/L), w0 = A sin(pi x/L) sin(pi y/L), equal biaxial uniform membrane compression N,
omega = W/h, a = A/h, tau = dT/dT_cr, dT_cr = 2 pi^2 h^2 / (12 (1 + nu) alpha L^2) (Amendment 1), W = TOTAL centre deflection.

Original (Amendment 2, own derivation, not cross-checked against a publication):
  N = N_cr (1 - A/W);  alpha dT - N (1 - nu)/(E h) = pi^2 (W^2 - A^2)/(8 L^2)
  -> omega^3 - [a^2 + c (tau - 1)] omega - c a = 0,  c = 4/(3 (1 + nu)).
  This first attempt neglects the non-uniform membrane stress induced by the deflection; the first V7 run failed against it
  (see docs/specification.md, Amendment 3).
Corrected (Amendment 3, own derivation, derived AFTER the first V7 result was known; not cross-checked against a publication):
  The compatibility equation with the imperfection gives the particular membrane forces N_x = -kappa cos(2 pi y/L),
  N_y = -kappa cos(2 pi x/L), kappa = E h pi^2 (W^2 - A^2)/(8 L^2), plus the uniform compression N. This field has zero normal
  edge displacement at every edge point and zero edge shear. Galerkin projection on the mode gives
  N = N_cr (1 - A/W) + kappa/2; the zero mean edge displacement gives alpha dT - pi^2 (W^2 - A^2)/(8 L^2) = N (1 - nu)/(E h).
  -> omega^3 - [a^2 + c' (tau - 1)] omega - c' a = 0,  c' = 8/(3 (1 + nu)(3 - nu)),
     N/N_cr = (1 - a/omega) + 3 (1 - nu^2) (omega^2 - a^2)/8   (N = mean edge force per unit width).
Functions reference_wh(), reference_wh_corrected(), reference_n_corrected(), dT_cr() are importable.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from material import MATERIAL

L = 0.5
H = 0.002
TAUS = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5]
AMPS = [0.01, 0.1, 0.5]


def dT_cr(h=H, length=L):
    m = MATERIAL
    return 2 * math.pi ** 2 * h ** 2 / (12 * (1 + m["nu"]) * m["alpha"] * length ** 2)


def c_coeff():
    return 4.0 / (3.0 * (1.0 + MATERIAL["nu"]))


def c_corrected():
    nu = MATERIAL["nu"]
    return 8.0 / (3.0 * (1.0 + nu) * (3.0 - nu))


def _root(a, tau, c):
    r = np.roots([1.0, 0.0, -(a * a + c * (tau - 1.0)), -c * a])
    pos = [z.real for z in r if abs(z.imag) < 1e-9 and z.real > 0]
    if len(pos) != 1:
        raise SystemExit("expected one positive root, found %r for a=%g tau=%g" % (pos, a, tau))
    return pos[0]


def reference_wh(a, tau):
    return _root(a, tau, c_coeff())


def reference_wh_corrected(a, tau):
    return _root(a, tau, c_corrected())


def reference_n_corrected(a, tau):
    nu = MATERIAL["nu"]
    w = reference_wh_corrected(a, tau)
    return (1.0 - a / w) + 3.0 * (1.0 - nu ** 2) * (w * w - a * a) / 8.0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "corrected"
    if mode not in ("original", "corrected"):
        raise SystemExit(__doc__)
    root = os.path.dirname(HERE)
    d = dT_cr()
    if mode == "original":
        out = ["V7 reference: single-mode imperfect-plate curve (immovable edges, uniform temperature rise), h = %.4f m, L = %.3f m" % (H, L),
               "dT_cr (closed form) = %.5f K, c = 4/(3(1+nu)) = %.5f; W = total centre deflection, imperfection included" % (d, c_coeff())]
        for a in AMPS:
            out.append("")
            out.append("imperfection A = %g h = %.4f mm" % (a, a * H * 1e3))
            out.append("  tau     dT [K]     W/h      W [mm]")
            for t in TAUS:
                w = reference_wh(a, t)
                out.append("  %.2f   %8.5f   %.4f   %.4f" % (t, t * d, w, w * H * 1e3))
        path = os.path.join(root, "results", "v7_reference.txt")
    else:
        out = ["V7 corrected reference (Amendment 3): single mode with the full von Karman membrane stress, h = %.4f m, L = %.3f m" % (H, L),
               "dT_cr (closed form) = %.5f K, c' = 8/(3(1+nu)(3-nu)) = %.5f; W = total centre deflection, imperfection included; N = mean edge force per width / N_cr" % (d, c_corrected()),
               "derived AFTER the first V7 result was known: post hoc for A = 0.1 h, blind check for A = 0.01 h and A = 0.5 h"]
        for a in AMPS:
            out.append("")
            out.append("imperfection A = %g h = %.4f mm (%s)" % (a, a * H * 1e3, "post hoc" if abs(a - 0.1) < 1e-12 else "blind"))
            out.append("  tau     dT [K]     W/h      W [mm]    N/N_cr")
            for t in TAUS:
                w = reference_wh_corrected(a, t)
                out.append("  %.2f   %8.5f   %.4f   %.4f   %.4f" % (t, t * d, w, w * H * 1e3, reference_n_corrected(a, t)))
        path = os.path.join(root, "results", "v7_reference_corrected.txt")
    text = "\n".join(out) + "\n"
    print(text, end="")
    with open(path, "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()
