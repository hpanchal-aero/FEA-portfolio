"""
Closed-form Euler plate-buckling reference for the Stage 3 *BUCKLE
verification case: 200x30x1.0mm plate, simply supported on all four
edges, uniaxial compression along the 200mm direction.

sigma_cr = k * pi^2 * E / (12*(1-nu^2)) * (t/b)^2
k = (m*b/a + a/(m*b))^2, minimized over positive integer m

Usage: python3 scripts/buckle_verify_reference.py
"""
import numpy as np

E = 73100.0
NU = 0.33
A = 200.0
B = 30.0
T = 1.0


def main():
    D = E * T**3 / (12 * (1 - NU**2))
    best_k, best_m = None, None
    for m in range(1, 15):
        k = (m*B/A + A/(m*B))**2
        if best_k is None or k < best_k:
            best_k, best_m = k, m
    sigma_cr = best_k * np.pi**2 * E / (12*(1-NU**2)) * (T/B)**2
    N_cr = sigma_cr * T
    P_cr = N_cr * B
    print(f"Plate {A}x{B}x{T}mm, simply supported, uniaxial compression")
    print(f"D={D:.4f} N*mm, m={best_m}, k={best_k:.6f}")
    print(f"sigma_cr={sigma_cr:.4f} MPa (NOTE: exceeds ~290-345 MPa Al 2024-T3 "
          f"yield -- this plate would yield before elastic buckling; noted, "
          f"not a concern for this eigenvalue-solver verification exercise)")
    print(f"N_cr={N_cr:.4f} N/mm, P_cr={P_cr:.4f} N (reference load for *BUCKLE)")


if __name__ == "__main__":
    main()
