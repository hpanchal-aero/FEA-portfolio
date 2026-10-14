#!/usr/bin/env python3
"""
Stage 6: append 8 replacement layouts (L069-L076) for the 8 layouts rejected by the sliver gate
(L004, L011, L019, L021, L026, L029, L030, L035). The rejected layouts stay in the table and records.
Exclusion band for interior rib centres widened to (685, 755) mm: the failing free pieces appeared for ribs
at about 695-722 mm. Same family as the rejected layout; random ones use a new seed (SEED + 2) and the
same 80 mm separation rule against all existing layouts. The first 69 rows are not modified.
Usage: /usr/bin/python3 scripts/append_replacement_layouts.py
"""
import csv
import random

SPAN, MIN_GAP, EXCL, SEED = 1000.0, 40.0, (685.0, 755.0), 20261011 + 2
K_MASS = 0.12104 * 0.00278
PATH = "results/stage6_layouts.csv"
REPL = [("L004", "power-law p=1.6"), ("L011", "handcrafted"), ("L019", "random"), ("L021", "random"),
        ("L026", "random"), ("L029", "random"), ("L030", "random"), ("L035", "random")]


def valid(z):
    pts = [0.0] + list(z) + [SPAN]
    if any(b - a < MIN_GAP - 1e-9 for a, b in zip(pts, pts[1:])):
        return False
    return not any(EXCL[0] < v < EXCL[1] for v in z)


def mass(z):
    return K_MASS * sum((300.0 - 0.15 * v) ** 2 for v in [0.0] + list(z) + [SPAN])


def main():
    raw = open(PATH).read()
    rows = list(csv.DictReader(open(PATH)))
    if len(rows) != 69:
        raise SystemExit("expected 69 rows, found %d: already appended? nothing written" % len(rows))
    fields = list(rows[0].keys())
    existing = [[float(r[k]) for k in ("z1", "z2", "z3", "z4")] for r in rows]
    byid = {r["id"]: r for r in rows}
    rng = random.Random(SEED)
    new = []

    def add(src, family, z, note):
        z = [round(float(v), 1) for v in z]
        if not valid(z):
            raise SystemExit("invalid replacement for %s: %s" % (src, z))
        existing.append(z)
        i = len(rows) + len(new)
        new.append(dict(id="L%03d" % i, family=family, z1=z[0], z2=z[1], z3=z[2], z4=z[3],
                        stations=",".join("%g" % v for v in [0.0] + z + [SPAN]),
                        est_rib_mass_g=round(mass(z), 1), note=note))

    for src, fam in REPL:
        old = byid[src]
        if fam == "power-law p=1.6":
            z = [SPAN * (i / 5.0) ** 1.6 for i in range(1, 5)]
            z = [685.0 if EXCL[0] < v < EXCL[1] else v for v in z]
            add(src, "power-law p=1.6 (replacement)", sorted(z), "replaces %s; rib moved to 685 (band edge)" % src)
        elif fam == "handcrafted":
            add(src, "handcrafted (replacement)", [150, 300, 680, 850], "replaces %s (rib 700 -> 680)" % src)
        else:
            for _ in range(200000):
                z = sorted(round(rng.uniform(40.0, 960.0), 1) for _ in range(4))
                if valid(z) and all(sum(abs(a - b) for a, b in zip(z, e)) >= 80.0 for e in existing):
                    break
            else:
                raise SystemExit("could not place a replacement for %s" % src)
            add(src, "random (replacement)", z, "replaces %s" % src)
    with open(PATH, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        for r in new:
            w.writerow(r)
    assert open(PATH).read().startswith(raw)
    for r in new:
        print(" %s  %-30s %-34s est rib mass %6.1f g  %s" % (r["id"], r["family"], r["stations"], r["est_rib_mass_g"], r["note"]))
    print("appended %d layouts; first 69 rows unchanged" % len(new))


if __name__ == "__main__":
    main()
