"""
Project 01 - Aerospace Mounting Bracket
Stage 3 (fillet study): extract peak stress at the two fillet surfaces
from the CalculiX .frd output, and compare against the locked Peterson
reference (Kt=1.429, sigma_nom=66.22 MPa, predicted peak=94.6 MPa).

Reuses the .frd-parsing approach validated in Stage 2's
diagnose_hole_stress.py.

Run from mesh/fillet_study/.
"""

import re
import numpy as np

FRD_FILE = "analysis_fillet.frd"
NSET_FILE = "nsets_fillet.inp"

SIGMA_NOM = 66.22   # MPa, locked reference (6M/(t*d^2) at x=30, beam theory)
KT_PETERSON = 1.429
PRED_PEAK = KT_PETERSON * SIGMA_NOM

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

def parse_frd_nodes_and_stress(path):
    """Parse node coordinates (2C block) and nodal stress (STRESS block,
    last occurrence = the real result, not intermediate substep)."""
    coords = {}
    stress = {}  # node_id -> (sxx, syy, szz, sxy, syz, szx)

    with open(path) as f:
        lines = f.readlines()

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("    2C"):
            i += 1
            while i < n and lines[i].startswith(" -1"):
                # Fixed-width CalculiX .frd node record: A3,I10,3E12.5
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                x = float(rec[13:25])
                y = float(rec[25:37])
                z = float(rec[37:49])
                coords[nid] = (x, y, z)
                i += 1
            continue
        if "STRESS" in line and line.startswith(" -4"):
            i += 1
            # skip the -5 component-definition lines
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                # Fixed-width CalculiX .frd stress record: A3,I10,6E12.5
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                vals = [float(rec[13 + 12*k:25 + 12*k]) for k in range(6)]
                stress[nid] = tuple(vals)
                i += 1
            continue
        i += 1

    return coords, stress

print(f"Reading NSETs from {NSET_FILE}")
top_nodes = read_nset(NSET_FILE, "fillet_top")
bottom_nodes = read_nset(NSET_FILE, "fillet_bottom")
print(f"fillet_top: {len(top_nodes)} nodes")
print(f"fillet_bottom: {len(bottom_nodes)} nodes")

print(f"\nReading {FRD_FILE}")
coords, stress = parse_frd_nodes_and_stress(FRD_FILE)
print(f"Parsed {len(coords)} node coordinates, {len(stress)} nodal stress results")

def report(label, node_set):
    rows = []
    for nid in node_set:
        if nid not in stress or nid not in coords:
            continue
        sxx, syy, szz, sxy, syz, szx = stress[nid]
        x, y, z = coords[nid]
        rows.append((nid, x, y, z, sxx))

    rows.sort(key=lambda r: abs(r[4]), reverse=True)

    print(f"\n=== {label}: top 10 by |sigma_xx| ===")
    print(f"{'node':>8} {'x':>8} {'y':>8} {'z':>8} {'sigma_xx (MPa)':>15}")
    for nid, x, y, z, sxx in rows[:10]:
        print(f"{nid:>8} {x:>8.3f} {y:>8.3f} {z:>8.3f} {sxx:>15.4f}")

    if rows:
        peak = rows[0]
        return peak
    return None

peak_top = report("fillet_top", top_nodes)
peak_bottom = report("fillet_bottom", bottom_nodes)

print("\n=== SUMMARY vs. LOCKED PETERSON REFERENCE ===")
print(f"sigma_nom (locked, 6M/(t*d^2) at x=30):  {SIGMA_NOM:.2f} MPa")
print(f"Kt (Peterson, r/d=0.25, D/d=1.5, h/r=1.0): {KT_PETERSON:.3f}")
print(f"Predicted peak stress:                     {PRED_PEAK:.2f} MPa")

for label, peak in [("Top fillet", peak_top), ("Bottom fillet", peak_bottom)]:
    if peak is None:
        print(f"{label}: no data found")
        continue
    nid, x, y, z, sxx = peak
    kt_observed = abs(sxx) / SIGMA_NOM
    pct_diff = 100.0 * (abs(sxx) - PRED_PEAK) / PRED_PEAK
    print(f"\n{label}: node {nid} at ({x:.3f}, {y:.3f}, {z:.3f})")
    print(f"  sigma_xx = {sxx:.4f} MPa")
    print(f"  Observed Kt = {kt_observed:.4f}")
    print(f"  vs. Peterson predicted {PRED_PEAK:.2f} MPa: {pct_diff:+.2f}%")
