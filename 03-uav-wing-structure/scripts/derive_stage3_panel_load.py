"""
Computes sigma_zz,skin(x, z=200) = -M(z)*y(x,z)/I(z) across the
root-bay panel's loaded edge (z=200mm, the interior rib), using the
SAME tapered-beam M(z) and spar cross-section I(z) already verified in
Stage 1/2 (spar's own bending stiffness -- ribs contribute ~0 per the
Stage 2 finding, so this is not a new quantity, just evaluated at the
panel's actual surface height y(x,z) instead of the flange's fixed
y=H(z)/2).

Prints the stress distribution and the total equivalent compressive
force on this edge, for review BEFORE this is converted into applied
nodal loads in the mesh script.

Usage: python3 scripts/derive_stage3_panel_load.py
"""
import math
import numpy as np

SPAN = 1000.0
ROOT_CHORD, TIP_CHORD = 300.0, 150.0
SPAR_CHORD_FRAC = 0.30
H_ROOT, H_TIP = 28.8, 14.4
WEB_T, CAP_W, CAP_T = 2.0, 8.0, 3.0
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
LIFT = 24.5
THICKNESS = 1.0  # mm, this panel-verify case only
Z_LOAD = 200.0
N_X = 30


def chord_at(z):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * (z / SPAN)


def naca0012_half_thickness(x_over_c, t=0.12):
    return 5 * t * (0.2969 * math.sqrt(x_over_c) - 0.1260 * x_over_c
                     - 0.3516 * x_over_c**2 + 0.2843 * x_over_c**3
                     - 0.1015 * x_over_c**4)


def H_of_z(z):
    return H_ROOT + (H_TIP - H_ROOT) * z / SPAN


def I_of_z(z):
    h = H_of_z(z)
    return (CAP_W * h**3 - (CAP_W - WEB_T) * (h - 2*CAP_T)**3) / 12.0


def A_of_z(z):
    h = H_of_z(z)
    return (h - 2*CAP_T) * WEB_T + 2 * CAP_W * CAP_T


def moment_at(z_query):
    n = 20001
    z = np.linspace(0.0, SPAN, n)
    dz = z[1] - z[0]
    q0 = LIFT / (np.pi * SPAN / 4.0)
    q = q0 * np.sqrt(np.clip(1 - (z / SPAN)**2, 0, None)) - RHO * A_of_z(z) * G
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i+1] + 0.5 * (q[i] + q[i+1]) * dz
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dz
    return float(np.interp(z_query, z, M))


def main():
    M_z = moment_at(Z_LOAD)
    I_z = I_of_z(Z_LOAD)
    c = chord_at(Z_LOAD)
    print(f"z={Z_LOAD}mm: M(z)={M_z:.4f} N*mm, I(z)={I_z:.4f} mm^4, chord={c:.2f}mm")
    print(f"(M(z), I(z) reused unchanged from the verified Stage 1/2 beam reference)\n")

    print(f"{'x/c':>6} {'x(mm)':>8} {'y(mm)':>8} {'sigma_zz(MPa)':>15}")
    total_force = 0.0
    prev_x = None
    xs, sigmas = [], []
    for i in range(N_X):
        frac = i / (N_X - 1)
        x_over_c = SPAR_CHORD_FRAC + frac * (1.0 - SPAR_CHORD_FRAC)
        x = x_over_c * c - SPAR_CHORD_FRAC * c
        y = c * naca0012_half_thickness(x_over_c)
        sigma = -M_z * y / I_z
        xs.append(x); sigmas.append(sigma)
        print(f"{x_over_c:6.3f} {x:8.2f} {y:8.3f} {sigma:15.4f}")

    xs = np.array(xs); sigmas = np.array(sigmas)
    total_force = np.trapz(sigmas, xs) * THICKNESS
    print(f"\nEquivalent total compressive force on this edge "
          f"(integral of sigma*t dx): {total_force:.4f} N")
    print(f"(For comparison, the flat-plate verification's closed-form "
          f"P_cr for a 200x30x1.0mm panel was 9017.4 N -- this is a "
          f"different geometry/loading, so no direct comparison is "
          f"expected, just a sanity-check magnitude.)")


if __name__ == "__main__":
    main()
