"""
Stage 4a (right-angle frame, sharp corner, no gusset) — closed-form
analytical reference, computed BEFORE any FEA (magnitudes locked from
the original derivation). Corrected AFTER first FEA comparison: the
original version had the axial-force sign and the
tension/compression face labeling backwards. Magnitudes were always
correct (statics moment/force values unaffected); only the sign
convention was wrong. Correction derived and cross-checked two
independent ways before accepting:

  1. Free-body equilibrium at a cut in Flange A: the internal axial
     force required to balance the tip load (applied toward the
     fixed support) is COMPRESSIVE, not tensile.
  2. Bend-closing physical intuition: the applied load closes the
     frame's interior angle slightly, putting the CONCAVE (inside)
     surface in compression and the CONVEX (outside) surface in
     tension - standard behavior for a bent member, analogous to an
     elbow.

Both checks agree with each other and with the FEA result once the
correct stress component is extracted per flange (see
diagnose_frame_results.py, which had a separate, independent bug:
using sigma_xx for Flange B instead of sigma_zz).

Geometry:
  Flange A: x in [0,60mm], z=0 (fixed at x=0)
  Flange B: z in [0,60mm], x=60 (free tip at z=60)
  Cross-section (both flanges): 40mm wide (y), 4mm thick (bending
  direction), constant. Same material/section as Stages 1-3.

Load: F = 235.44 N at the Flange B tip, in the -x direction.

Reentrant (concave) corner is at (x=58, z=2) - i.e. Flange A's z=+2
face and Flange B's x=58 face are the CONCAVE (inside-of-bend) side.
Flange A's z=-2 face and Flange B's x=62 face are the CONVEX
(outside-of-bend) side.
"""

F = 235.44       # N
L = 60.0         # mm, length of each flange
E = 71700.0      # MPa
WIDTH = 40.0     # mm
THICK = 4.0      # mm

A = WIDTH * THICK
I = WIDTH * THICK**3 / 12.0
c = THICK / 2.0

print("=== Cross-section properties (same as Stages 1-3) ===")
print(f"A = {A} mm^2")
print(f"I = {I:.6e} mm^4")
print(f"c = {c} mm")

print("\n=== Flange B (transverse bending, corner at s=L) ===")
print("Concave (inside of bend) face: x=58   Convex (outside): x=62")
print(f"{'s (mm, from tip)':>18} {'M (N.mm)':>12} {'concave sigma_zz (MPa)':>24} {'convex sigma_zz (MPa)':>24}")
for s in [0, 15, 30, 45, 60]:
    M = F * s
    sigma_mag = M * c / I
    print(f"{s:>18.1f} {M:>12.2f} {-sigma_mag:>24.4f} {sigma_mag:>24.4f}")

print("\n=== Flange A (constant bending + constant axial, corrected signs) ===")
print("Concave (inside of bend) face: z=+2   Convex (outside): z=-2")
M_A = F * L
N_A = -F   # CORRECTED: compressive, not tensile
sigma_bend_mag = M_A * c / I
sigma_axial = N_A / A
print(f"Constant M = {M_A:.2f} N.mm along entire flange")
print(f"Constant N (axial) = {N_A:.2f} N along entire flange (compressive)")
print(f"sigma_bending magnitude (surface) = {sigma_bend_mag:.4f} MPa")
print(f"sigma_axial = {sigma_axial:.4f} MPa")
sigma_concave = -sigma_bend_mag + sigma_axial   # inside: bending compression + axial compression
sigma_convex = sigma_bend_mag + sigma_axial     # outside: bending tension + axial compression
print(f"sigma_total, concave (z=+2, inside):  {sigma_concave:.4f} MPa (compression)")
print(f"sigma_total, convex  (z=-2, outside): {sigma_convex:.4f} MPa (tension)")

# ---- Tip deflection - unaffected by the sign correction (a scalar
# magnitude derived from strain energy, not sensitive to the
# tension/compression face-labeling error) ----
delta_B = F * L**3 / (3 * E * I)
delta_A_bending = F * L**3 / (E * I)
delta_A_axial = F * L / (E * A)
delta_total = delta_B + delta_A_bending + delta_A_axial

print("\n=== Tip deflection (Castigliano / unit-load method, unchanged) ===")
print(f"Flange B bending contribution:  {delta_B:.4f} mm")
print(f"Flange A bending contribution:  {delta_A_bending:.4f} mm")
print(f"Flange A axial contribution:    {delta_A_axial:.4f} mm")
print(f"TOTAL tip deflection magnitude: {delta_total:.4f} mm (direction: -x, matching applied load)")
