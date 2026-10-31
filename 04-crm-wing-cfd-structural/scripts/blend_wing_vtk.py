#!/usr/bin/python3
"""Blend the wing-patch pressure of two CFD states on an identical mesh to a target lift coefficient.

usage: blend_wing_vtk.py <vtk1> <alpha1_deg> <vtk2> <alpha2_deg> <target_CL> <out.vtk>

Cl is computed from the surface pressure: F = sum (p - pinf) A n (outward from the fluid domain, as in map_cfd_loads.py),
L = Fz cos(alpha) - Fx sin(alpha), Cl = L / (q Aref), q = 10955.9 Pa, Aref = 191.870 m2 (half model).
p_blend = (1 - w) p1 + w p2 with w = (CL_target - Cl1)/(Cl2 - Cl1).  Only the cell array p is blended;
all other arrays and the geometry are those of state 2.  Interpolation is exact for the (linear) structure;
the assumption is linearity of the aerodynamic load in angle of attack between the two states.
"""
import os
import sys

import numpy as np
import pyvista as pv

PINF, Q, AREF = 21662.7, 10955.9, 191.870


def surf(path):
    m = pv.read(path)
    s = m.compute_cell_sizes(length=False, area=True, volume=False)
    s = s.compute_normals(cell_normals=True, point_normals=False, consistent_normals=False,
                          auto_orient_normals=False, split_vertices=False)
    return m, np.asarray(s.cell_data["Area"]), np.asarray(s.cell_data["Normals"]), np.asarray(m.cell_data["p"])


def cl_of(area, nrm, p, alpha):
    F = (((p - PINF) * area)[:, None] * nrm).sum(axis=0)
    if -F[0] * np.sin(np.radians(alpha)) + F[2] * np.cos(np.radians(alpha)) < 0:
        F = -F
    a = np.radians(alpha)
    return (F[2] * np.cos(a) - F[0] * np.sin(a)) / (Q * AREF), F


v1, a1, v2, a2, target, out = sys.argv[1], float(sys.argv[2]), sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
m1, A1, N1, p1 = surf(v1)
m2, A2, N2, p2 = surf(v2)
if m1.n_cells != m2.n_cells or m1.n_points != m2.n_points or not np.allclose(m1.points, m2.points):
    sys.exit("ABORT: the two wing patches are not the same mesh")
cl1, F1 = cl_of(A1, N1, p1, a1)
cl2, F2 = cl_of(A2, N2, p2, a2)
w = (target - cl1) / (cl2 - cl1)
print(f"state 1: {v1}  alpha {a1}  Cl {cl1:.5f}  F {F1[0]:.1f} {F1[1]:.1f} {F1[2]:.1f} N")
print(f"state 2: {v2}  alpha {a2}  Cl {cl2:.5f}  F {F2[0]:.1f} {F2[1]:.1f} {F2[2]:.1f} N")
print(f"target Cl {target}  weight on state 2: w = {w:.5f}   (1 - w = {1 - w:.5f})")
if not (0.0 <= w <= 1.0):
    sys.exit("ABORT: target lift is outside the bracket (extrapolation), not done silently")
pb = (1.0 - w) * p1 + w * p2
aeq = (1.0 - w) * a1 + w * a2
clb, Fb = cl_of(A2, N2, pb, aeq)
print(f"blend:   equivalent alpha {aeq:.4f}   Cl {clb:.5f}   F {Fb[0]:.1f} {Fb[1]:.1f} {Fb[2]:.1f} N")
m2.cell_data["p"] = pb
os.makedirs(os.path.dirname(out), exist_ok=True)
m2.save(out)
print("wrote", out)
