#!/usr/bin/env python3
"""Read-only audit of the uCRM-9 NASTRAN inputs for anything the CalculiX converter ignores.

  1. card inventory, executive/case-control lines, LOAD/GRAV/SPC1/SPCADD/MPC/PARAM cards of every file
  2. wingbox CQUADR: non-zero ZOFFS, non-zero continuation fields (TFLAG, T1..T4); first raw cards
  3. wingbox PSHELL: 12I/T3, TS/T, Z1, Z2, NSM, MID2, MID4 usage
  4. SPC cards grouped by set id (the converter merges all of them)
  5. aerodynamic FORCE cards grouped by (set id, coordinate system)
Usage: /usr/bin/python3 scripts/check_bdf_setup.py [coarse|medium|fine]
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

level = sys.argv[1] if len(sys.argv) > 1 else "coarse"
if level not in ("coarse", "medium", "fine"):
    sys.exit("level must be coarse, medium or fine")
fem = HERE.parent / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"
if not fem.is_dir():
    sys.exit(f"ABORT: {fem} not found")


def blank(s):
    return s.strip() == ""


def nonzero(s):
    if blank(s):
        return False
    try:
        return num(s) != 0.0
    except Exception:
        return True


def card_name(ln):
    if not ln or ln[0] in "$+* \t":
        return None
    head = ln.split(",")[0] if "," in ln[:9] else ln[:8]
    return head.strip().rstrip("*").upper()


print(f"== files in {fem}")
for p in sorted(x for x in fem.iterdir() if x.is_file()):
    raw = p.read_bytes()
    if b"\x00" in raw[:4096]:
        print(f"\n{p.name}: binary, skipped")
        continue
    lines = raw.decode("utf-8", errors="replace").splitlines()
    bulk = next((k for k, ln in enumerate(lines) if re.match(r"\s*BEGIN\s+BULK", ln, re.I)), None)
    inv = Counter()
    for k, ln in enumerate(lines):
        if bulk is not None and k <= bulk:
            continue
        nm = card_name(ln)
        if nm:
            inv[nm] += 1
    print(f"\n{p.name}: {len(lines)} lines")
    print("  cards: " + ", ".join(f"{n}={c}" for n, c in inv.most_common(18)))
    if bulk is not None:
        ctl = [ln for ln in lines[:bulk + 1] if ln.strip() and not ln.startswith("$")]
        print(f"  executive/case control ({len(ctl)} lines up to BEGIN BULK):")
        for ln in ctl[:60]:
            print("    " + ln[:100])
    shown = Counter()
    for k, ln in enumerate(lines):
        if bulk is not None and k <= bulk:
            continue
        nm = card_name(ln)
        if nm in ("LOAD", "GRAV", "SPCADD", "SPC1", "MPC", "MPCADD", "PARAM", "EIGRL", "TEMP") and shown[nm] < 6:
            shown[nm] += 1
            print("  " + ln[:100])

wb = fem / f"uCRM-9_wingbox_{level}.bdf"
cards = read_cards(wb)
print(f"\n== wingbox {wb.name}")
nq = zoff = cont = 0
for name, f in cards:
    if name == "CQUADR":
        nq += 1
        zoff += nonzero(f[7])
        cont += any(nonzero(x) for x in f[8:])
print(f"CQUADR {nq}: non-zero ZOFFS {zoff}; non-zero continuation fields (TFLAG/T1-T4) {cont}")
text = wb.read_text(errors="replace").splitlines()
shown = 0
for k, ln in enumerate(text):
    if card_name(ln) == "CQUADR" and shown < 2:
        print("  raw:", repr(ln))
        print("      ", repr(text[k + 1]) if k + 1 < len(text) else "(end of file)")
        shown += 1

ps = [f for name, f in cards if name == "PSHELL"]


def top(i):
    return Counter(f[i].strip() or "(blank)" for f in ps).most_common(5)


print(f"PSHELL {len(ps)}: 12I/T3 {top(4)}")
print(f"  TS/T {top(6)}")
print(f"  Z1 {top(8)}  Z2 {top(9)}")
print(f"  NSM non-zero {sum(nonzero(f[7]) for f in ps)}   MID2 set {sum(not blank(f[3]) for f in ps)}   "
      f"MID4 set {sum(not blank(f[10]) for f in ps)}")

sets = defaultdict(lambda: [0, Counter()])
for name, f in cards:
    if name == "SPC":
        i = 1
        while i + 1 < len(f):
            if not blank(f[i]) and ival(f[i]) > 0:
                sets[f[0].strip()][0] += 1
                sets[f[0].strip()][1][f[i + 1].strip()] += 1
            i += 3
print("SPC sets:")
for sid, (n, c) in sorted(sets.items()):
    print(f"  SID {sid}: {n} constrained entries; components {c.most_common(6)}")

aero = fem / f"uCRM-9_aero_forces_1.0g_{level}.bdf"
rows = [(f[0].strip(), f[2].strip() or "0", num(f[3])) for n, f in read_cards(aero) if n == "FORCE"]
print(f"\n== {aero.name}: {len(rows)} FORCE cards; (set id, CID) -> count {dict(Counter((s, c) for s, c, _ in rows))}")
if rows:
    print(f"   scale F: min {min(r[2] for r in rows):.6g}  max {max(r[2] for r in rows):.6g}")
