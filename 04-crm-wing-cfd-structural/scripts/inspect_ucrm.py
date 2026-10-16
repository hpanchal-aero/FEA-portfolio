#!/usr/bin/env python3
"""Read-only inspection of the uCRM-9 dataset: BDF card inventory, material/property
cards, run decks, load files, and the two reference spreadsheets."""
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "data" / "ucrm9" / "uCRM9"


def card_counts(path):
    c = Counter()
    for line in path.read_text(errors="replace").splitlines():
        if not line or line[0] in "$+ \t":
            continue
        key = re.split(r"[ ,*\t]", line, maxsplit=1)[0].upper()
        if key and key[0].isalpha():
            c[key] += 1
    return c


def grep_cards(path, keys, limit):
    out = []
    for line in path.read_text(errors="replace").splitlines():
        if not line or line[0] in "$+ \t":
            continue
        key = re.split(r"[ ,*\t]", line, maxsplit=1)[0].upper()
        if key in keys:
            out.append(line.rstrip())
            if len(out) >= limit:
                break
    return out


def head(path, n):
    return path.read_text(errors="replace").splitlines()[:n]


print("=== ROOT ===", ROOT, "exists:", ROOT.is_dir())

for lvl in ("coarse", "medium", "fine"):
    p = ROOT / "Mesh" / lvl.capitalize() / "FEM" / f"uCRM-9_wingbox_{lvl}.bdf"
    print(f"\n=== card inventory: {p.name} ===")
    for k, v in sorted(card_counts(p).items(), key=lambda kv: -kv[1]):
        print(f"  {k:10s} {v}")

coarse = ROOT / "Mesh" / "Coarse" / "FEM"
wb = coarse / "uCRM-9_wingbox_coarse.bdf"

print("\n=== MAT / PSHELL / PCOMP / PBAR / SPC cards (coarse, first 12 each) ===")
for keys in (("MAT1", "MAT8", "MAT2"), ("PSHELL",), ("PCOMP",), ("PBAR", "PBARL", "PBEAM", "PROD"), ("SPC", "SPC1", "SPCADD")):
    print(" ", keys)
    for l in grep_cards(wb, set(keys), 12):
        print("   ", l)

print("\n=== first 25 lines of coarse wingbox BDF ===")
for l in head(wb, 25):
    print("  ", l)

for name in ("run_uCRM-9_static_cruise_coarse.bdf", "run_uCRM-9_modal_coarse.bdf"):
    print(f"\n=== {name} (full) ===")
    for l in head(coarse / name, 60):
        print("  ", l)

for name, n in (("uCRM-9_aero_forces_1.0g_coarse.bdf", 12), ("uCRM-9_engine_mass_coarse.bdf", 10), ("uCRM-9_LETE_masses_coarse.bdf", 8)):
    print(f"\n=== {name} (first {n} lines) ===")
    for l in head(coarse / name, n):
        print("  ", l)
    print("   card counts:", dict(card_counts(coarse / name)))

print("\n=== uCRM-9_wingbox_panel_sizing.txt (first 20 lines) ===")
for l in head(ROOT / "Geometry" / "uCRM-9_wingbox_panel_sizing.txt", 20):
    print("  ", l)

print("\n=== spreadsheets ===")
try:
    import openpyxl
except ImportError:
    sys.exit("openpyxl not importable: tell me and I will use a zipfile/XML fallback")

for name in ("uCRM-9_Design_Specifications.xlsx", "uCRM-9_Results.xlsx"):
    wbk = openpyxl.load_workbook(ROOT / "Data" / name, data_only=True)
    for ws in wbk.worksheets:
        print(f"\n--- {name} :: sheet '{ws.title}' ({ws.max_row} rows x {ws.max_column} cols) ---")
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            if i >= 70:
                print("   ... (truncated)")
                break
            if any(v is not None for v in row):
                print("  ", [v for v in row])
