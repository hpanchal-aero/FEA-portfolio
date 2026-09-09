"""
Stage 4a (right-angle frame) - full result extraction: reaction
equilibrium, tip displacement, and stress at the stations locked in
frame_analytical_reference.py, compared against the CORRECTED
closed-form targets (see that script's docstring for the sign-error
history).

Fixes two bugs from the first pass:
  1. Flange B's bending stress component is sigma_zz, not sigma_xx
     (its beam axis is z) - corrected here.
  2. Tension/compression face labels and the axial-force sign were
     backwards in the original analytical reference - corrected
     there; this script now compares against the corrected targets.
"""

import numpy as np

FRD_FILE = "analysis_frame.frd"
NSET_FILE = "nsets_frame.inp"

F = 235.44
L = 60.0

# Locked, CORRECTED analytical targets (magnitude progression
# unchanged from the original derivation; concave/convex sign
# assignment corrected)
TARGETS_B_MAG = {0: 0.0, 15: 33.1088, 30: 66.2175, 45: 99.3262, 60: 132.4350}
# Flange B concave (x=58) = -magnitude, convex (x=62) = +magnitude

TARGET_A_CONCAVE = -133.9065   # z=+2, compression
TARGET_A_CONVEX = 130.9635     # z=-2, tension
TARGET_TIP_DEFL = 4.4342


def read_nset(path, name):
    node_ids = []
    with open(path) as f:
        lines = f.readlines()
    capture = False
    for line in lines:
        s = line.strip()
        if s.upper().startswith("*NSET"):
            capture = (f"NSET={name}".upper() in s.upper().replace(" ", ""))
            continue
        if s.startswith("*"):
            capture = False
            continue
        if capture:
            parts = [p.strip() for p in s.split(",") if p.strip()]
            node_ids.extend(int(p) for p in parts)
    return set(node_ids)


def parse_frd(path):
    coords, disp, stress, rf = {}, {}, {}, {}
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("    2C"):
            i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                coords[nid] = (float(rec[13:25]), float(rec[25:37]), float(rec[37:49]))
                i += 1
            continue
        if "DISP" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                disp[nid] = (float(rec[13:25]), float(rec[25:37]), float(rec[37:49]))
                i += 1
            continue
        if "STRESS" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                vals = [float(rec[13 + 12*k:25 + 12*k]) for k in range(6)]
                stress[nid] = tuple(vals)
                i += 1
            continue
        if "FORC" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                rf[nid] = (float(rec[13:25]), float(rec[25:37]), float(rec[37:49]))
                i += 1
            continue
        i += 1
    return coords, disp, stress, rf


print(f"Reading {FRD_FILE}")
coords, disp, stress, rf = parse_frd(FRD_FILE)
print(f"Parsed {len(coords)} coords, {len(disp)} disp, {len(stress)} stress, {len(rf)} RF")

root_nodes = read_nset(NSET_FILE, "root")
tip_nodes = read_nset(NSET_FILE, "tip")

# ---- 1. Equilibrium check (unaffected by the fix, reported again for completeness) ----
root_rf = [rf[n] for n in root_nodes if n in rf]
sum_fx = sum(r[0] for r in root_rf)
sum_fy = sum(r[1] for r in root_rf)
sum_fz = sum(r[2] for r in root_rf)
print("\n=== EQUILIBRIUM CHECK ===")
print(f"Root Sum Fx: {sum_fx:.4f} N  (expect ~+235.44)")
print(f"Root Sum Fy: {sum_fy:.4f} N  (expect ~0)")
print(f"Root Sum Fz: {sum_fz:.4f} N  (expect ~0)")
err = 100 * (abs(sum_fx) - F) / F
print(f"Relative equilibrium error (Fx): {err:+.4f}%")

# ---- 2. Tip displacement ----
tip_ux = [disp[n][0] for n in tip_nodes if n in disp]
mean_ux = np.mean(tip_ux)
print("\n=== TIP DISPLACEMENT ===")
print(f"Mean Ux at tip: {mean_ux:.4f} mm")
print(f"Analytical target: {-TARGET_TIP_DEFL:.4f} mm (magnitude {TARGET_TIP_DEFL})")
pct = 100 * (mean_ux - (-TARGET_TIP_DEFL)) / TARGET_TIP_DEFL
print(f"% diff: {pct:+.2f}%")

# ---- 3. Flange B stress (sigma_zz, index 2) at both concave/convex faces ----
print("\n=== FLANGE B STRESS (sigma_zz) ===")
print(f"{'s(mm)':>6} {'z':>6} {'face':>8} {'theory(MPa)':>12} {'FEA(MPa)':>12} {'%diff':>8}")
for s, mag in TARGETS_B_MAG.items():
    if s == 0:
        continue
    z_target = L - s
    for face_x, face_name, sign in [(58.0, "concave", -1), (62.0, "convex", +1)]:
        theory = sign * mag
        candidates = []
        for nid, (x, y, z) in coords.items():
            if abs(z - z_target) < 0.75 and abs(y - 20.0) < 1.5 and abs(x - face_x) < 0.5:
                if nid in stress:
                    candidates.append(stress[nid][2])  # sigma_zz, index 2
        if candidates:
            szz = np.mean(candidates)
            pct = 100 * (szz - theory) / theory if theory != 0 else float('nan')
            print(f"{s:>6.0f} {z_target:>6.1f} {face_name:>8} {theory:>12.4f} {szz:>12.4f} {pct:>7.2f}%  (n={len(candidates)})")
        else:
            print(f"{s:>6.0f} {z_target:>6.1f} {face_name:>8} {theory:>12.4f} {'NO MATCH':>12}")

# ---- 4. Flange A stress at mid-span (x=30) ----
print("\n=== FLANGE A STRESS (sigma_xx) at x=30 (mid-span) ===")
concave = []  # z = +2
convex = []   # z = -2
for nid, (x, y, z) in coords.items():
    if abs(x - 30.0) < 0.5 and abs(y - 20.0) < 1.0:
        if nid not in stress:
            continue
        if abs(z - 2.0) < 0.3:
            concave.append(stress[nid][0])
        elif abs(z - (-2.0)) < 0.3:
            convex.append(stress[nid][0])

if concave:
    mean_c = np.mean(concave)
    pct_c = 100 * (mean_c - TARGET_A_CONCAVE) / TARGET_A_CONCAVE
    print(f"Concave (z=+2) sigma_xx: {mean_c:.4f} MPa vs target {TARGET_A_CONCAVE:.4f} MPa ({pct_c:+.2f}%), n={len(concave)}")
else:
    print("Concave: NO MATCH")

if convex:
    mean_v = np.mean(convex)
    pct_v = 100 * (mean_v - TARGET_A_CONVEX) / TARGET_A_CONVEX
    print(f"Convex (z=-2) sigma_xx: {mean_v:.4f} MPa vs target {TARGET_A_CONVEX:.4f} MPa ({pct_v:+.2f}%), n={len(convex)}")
else:
    print("Convex: NO MATCH")
