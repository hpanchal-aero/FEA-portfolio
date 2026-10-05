#!/usr/bin/env python3
"""
Stage 4: single-cell Bredt-Batho tip twist estimate (hand-calculation reference, NOT a verification).
Assumptions: closed thin-walled single cell, constant t = 1.0 mm, NACA 0012 section scaled with chord,
ribs and spar web ignored, G = E/(2(1+nu)), net lift (elliptical span, centre of pressure at 25% chord)
acting at the section; shear centre assumed at x_sc (fraction of chord from LE; spar axis = 0.30).
Torque per unit span = q(z) * (x_sc - 0.25) * c(z), nose-up for a lift ahead of the shear centre.
Usage: /usr/bin/python3 scripts/estimate_stage4_twist_bredt.py
"""
import numpy as np

E, NU, T = 73100.0, 0.33, 1.0
G = E / (2.0 * (1.0 + NU))
L, SPAN = 372.78, 1000.0
q0 = L / (SPAN * np.pi / 4.0)

x = np.linspace(0.0, 1.0, 200001)
yt = 0.6 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
a = 2.0 * np.trapz(yt, x)                                   # enclosed area / c^2
p = 2.0 * np.sum(np.hypot(np.diff(x), np.diff(yt)))         # perimeter / c
print(f"section constants (closed-TE NACA 0012): A/c^2 = {a:.5f}, s/c = {p:.4f}, G = {G:.1f} MPa")

z = np.linspace(0.0, SPAN, 200001)
dz = z[1] - z[0]
c = 300.0 - 0.15 * z
q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN) ** 2, 0.0, 1.0))
GJ = G * 4.0 * (a * c ** 2) ** 2 * T / (p * c)
print(f"GJ root = {GJ[0]:.4e} N*mm^2, GJ tip = {GJ[-1]:.4e} N*mm^2")
for xsc in (0.25, 0.30, 0.35, 0.40):
    m = q * (xsc - 0.25) * c
    Tq = (np.cumsum(m[::-1]) * dz)[::-1]                    # torque at z from loads outboard of z
    th = np.cumsum(Tq / GJ) * dz
    row = "  ".join(f"z={zr:4.0f}: {th[int(zr / dz)]:.3e}" for zr in (200, 400, 600, 800, 1000))
    print(f"x_sc = {xsc:.2f} c (arm {xsc - 0.25:+.2f} c): theta [rad] {row}")
