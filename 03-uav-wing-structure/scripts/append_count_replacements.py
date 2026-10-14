#!/usr/bin/env python3
"""
Stage 6B: replacements L152-L155 for L141, L142, L147, L148 (rejected by the sliver gate: free piece 4.2-4.6 mm^2
< 5 mm^2 for a rib at 688-690 mm). The offending rib is moved to 680 mm (L070 passed there); other ribs unchanged.
The rejected layouts stay in the table. Requires exactly 129 data rows; nothing else is modified.
Usage: /usr/bin/python3 scripts/append_count_replacements.py
"""
import csv

PATH = "results/stage6_layouts.csv"
K_MASS = 0.12104 * 0.00278
NEW = [("L141", "7 ribs L006 + outboard (replacement)", [60, 130, 200, 400, 680]),
       ("L142", "8 ribs uniform (replacement, rib 690 -> 680)", [142.9, 285.7, 428.6, 571.4, 680, 857.1]),
       ("L147", "9 ribs power p=1.3 (replacement, rib 688 -> 680)", [67, 164.9, 279.4, 406.1, 542.8, 680, 840.6]),
       ("L148", "9 ribs L006 + outboard (replacement)", [60, 130, 200, 400, 550, 680, 850])]


def main():
    raw = open(PATH, newline="").read()
    rows = list(csv.DictReader(open(PATH, newline="")))
    if len(rows) != 129:
        raise SystemExit("expected 129 rows, found %d: nothing written" % len(rows))
    byid = {r["id"]: r for r in rows}
    fields = list(rows[0].keys())
    out = []
    for i, (src, fam, z) in enumerate(NEW):
        pts = [0.0] + z + [1000.0]
        if any(b - a < 40.0 - 1e-9 for a, b in zip(pts, pts[1:])):
            raise SystemExit("gap violation in " + fam)
        if src not in byid:
            raise SystemExit("source layout %s missing" % src)
        out.append(dict(id="L%03d" % (152 + i), family=fam, z1="", z2="", z3="", z4="",
                        stations=",".join("%g" % v for v in pts),
                        est_rib_mass_g=round(K_MASS * sum((300.0 - 0.15 * v) ** 2 for v in pts), 1),
                        note="replaces %s" % src))
    with open(PATH, "a", newline="") as f:
        if not raw.endswith("\n"):
            f.write("\r\n")
        w = csv.DictWriter(f, fieldnames=fields)
        for r in out:
            w.writerow(r)
    assert open(PATH, newline="").read().startswith(raw)
    for r in out:
        print(" %s  %-52s %-52s %6.1f g  %s" % (r["id"], r["family"], r["stations"], r["est_rib_mass_g"], r["note"]))
    print("appended 4 layouts; existing rows unchanged")


if __name__ == "__main__":
    main()
