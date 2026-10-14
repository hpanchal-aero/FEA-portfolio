#!/usr/bin/env python3
"""
Stage 6B (rib-count study): layouts with 2-5 ribs (L100-L137) and 7-10 ribs (L138...) in total.
Root rib (z=0) and tip rib (z=1000) fixed; rib thickness 1.5 mm. Heuristic search, not a global optimum.
Rules: min gap 40 mm; no interior rib centre in (690, 745) mm (slivers appeared at ribs 695-722; 750 passed twice);
ribs falling in the band are nudged to the nearest edge and labelled; layouts still invalid are skipped and reported.
Safe to rerun: with 77 rows it appends everything; with the 38 low-count rows already present it verifies them and
appends only the rest; any other row count aborts without writing.
Usage: /usr/bin/python3 scripts/append_count_layouts.py
"""
import csv

PATH = "results/stage6_layouts.csv"
SPAN, MIN_GAP, BAND = 1000.0, 40.0, (690.0, 745.0)
K_MASS = 0.12104 * 0.00278

LOW = [
    ("2 ribs: root+tip only", [[]]),
    ("3 ribs", [[p] for p in (100, 150, 200, 250, 300, 400, 500, 600, 800)]),
    ("4 ribs uniform", [[333.3, 666.7]]),
    ("4 ribs", [[60, 130], [80, 160], [100, 200], [100, 300], [150, 300], [150, 400], [200, 400],
                [100, 400], [60, 200], [200, 500], [300, 600], [400, 600], [120, 250], [50, 100]]),
    ("5 ribs uniform", [[250, 500, 750]]),
    ("5 ribs", [[60, 130, 200], [80, 160, 240], [100, 200, 300], [60, 130, 400], [100, 250, 400],
                [150, 300, 450], [100, 200, 400], [200, 400, 600], [150, 300, 600], [100, 300, 500],
                [80, 200, 350], [50, 100, 150]]),
]


def valid(z):
    pts = [0.0] + list(z) + [SPAN]
    if any(b - a < MIN_GAP - 1e-9 for a, b in zip(pts, pts[1:])):
        return False
    return not any(BAND[0] < v < BAND[1] for v in z)


def mass(z):
    return K_MASS * sum((300.0 - 0.15 * v) ** 2 for v in [0.0] + list(z) + [SPAN])


def snap(z):
    out, nudged = [], False
    for v in z:
        if BAND[0] < v < BAND[1]:
            v = BAND[0] if v - BAND[0] <= BAND[1] - v else BAND[1]
            nudged = True
        out.append(round(v, 1))
    return sorted(out), nudged


def build():
    items, skipped = [], []
    for fam, lst in LOW:
        for z in lst:
            z = [round(float(v), 1) for v in z]
            if not valid(z):
                raise SystemExit("invalid low-count layout %s %s" % (fam, z))
            items.append((fam, z))
    for n in (7, 8, 9, 10):
        k = n - 2
        cands = [("uniform", [SPAN * i / (n - 1) for i in range(1, n - 1)]),
                 ("power p=1.3", [SPAN * (i / (n - 1)) ** 1.3 for i in range(1, n - 1)]),
                 ("power p=1.6", [SPAN * (i / (n - 1)) ** 1.6 for i in range(1, n - 1)]),
                 ("L006 + outboard", [60.0, 130.0, 200.0, 400.0] + [400.0 + 600.0 * j / (k - 4 + 1) for j in range(1, k - 4 + 1)])]
        for name, z in cands:
            z, nudged = snap(z)
            label = "%d ribs %s%s" % (n, name, " (nudged out of band)" if nudged else "")
            if valid(z):
                items.append((label, z))
            else:
                skipped.append((label, z))
    return items, skipped


def main():
    raw = open(PATH, newline="").read()
    rows = list(csv.DictReader(open(PATH, newline="")))
    fields = list(rows[0].keys())
    items, skipped = build()
    n_low = sum(len(l) for _, l in LOW)
    if len(rows) == 77:
        have = {}
    elif len(rows) == 77 + n_low:
        have = {r["id"]: r for r in rows[77:]}
    else:
        raise SystemExit("unexpected row count %d (expected 77 or %d): nothing written" % (len(rows), 77 + n_low))
    new = []
    for i, (fam, z) in enumerate(items):
        lid = "L%03d" % (100 + i)
        st = ",".join("%g" % v for v in [0.0] + z + [SPAN])
        if lid in have:
            if have[lid]["stations"] != st:
                raise SystemExit("existing %s differs from the plan (%s vs %s): nothing written" % (lid, have[lid]["stations"], st))
            continue
        row = dict(id=lid, family=fam, z1="", z2="", z3="", z4="", stations=st,
                   est_rib_mass_g=round(mass(z), 1), note="rib-count study")
        if len(z) <= 4:
            for j, v in enumerate(z, 1):
                row["z%d" % j] = v
        new.append(row)
    with open(PATH, "a", newline="") as f:
        if not raw.endswith("\n"):
            f.write("\r\n")
        w = csv.DictWriter(f, fieldnames=fields)
        for r in new:
            w.writerow(r)
    assert open(PATH, newline="").read().startswith(raw)
    for r in new:
        print(" %s  %-34s %-60s %6.1f g" % (r["id"], r["family"], r["stations"], r["est_rib_mass_g"]))
    for label, z in skipped:
        print(" SKIPPED (invalid after nudging): %s %s" % (label, z))
    print("appended %d layouts; existing rows unchanged; skipped %d" % (len(new), len(skipped)))


if __name__ == "__main__":
    main()
