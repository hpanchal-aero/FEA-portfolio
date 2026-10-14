#!/usr/bin/env python3
"""
Stage 6 layout set: 4 interior rib stations between fixed ribs at z = 0 and z = 1000, 6 ribs in total.
Constraints: minimum gap 40 mm between neighbouring ribs (including the root and tip ribs); no interior rib centre
in the exclusion zone (715, 750) mm around the apex of the inherited wedge piece (z = 742.28), where a rib edge
just below the apex leaves a sliver remnant (Stage 5, 20 mm). The driver's imprint gate is the final check.
Layouts: L000 baseline; L001-L005 power-law grading; L006-L013 handcrafted; L014-L016 baseline shifted by +1, +2, -2 mm
(mesh-discretization-scatter estimate); L017-L056 random (fixed seed, mutually at least 80 mm apart in the sum of |dz|);
L057-L062 root-biased random (all 4 interior ribs in z = 40-400 mm); L063-L068 root-biased random (3 interior ribs in
z = 40-400 mm and 1 in z = 400-960 mm). The root-biased family has its own seed so L000-L056 are unchanged.
Writes results/stage6_layouts.csv. Estimated rib mass = 0.12104 * sum(chord^2) * 0.00278 g (chord = 300 - 0.15 z),
calibrated on the 6-rib baseline (107.5 g); the exact mass comes from the geometry builder.
Usage: /usr/bin/python3 scripts/make_layouts.py
"""
import csv
import os
import random

SPAN = 1000.0
MIN_GAP = 40.0
EXCL = (715.0, 750.0)
SEED = 20261011
SEED_ROOT = SEED + 1
K_MASS = 0.12104 * 0.00278
N_RANDOM_END = 57        # L017-L056 random
N_TOTAL = 69             # + 12 root-biased


def valid(z):
    pts = [0.0] + list(z) + [SPAN]
    if any(b - a < MIN_GAP - 1e-9 for a, b in zip(pts, pts[1:])):
        return False
    return not any(EXCL[0] < v < EXCL[1] for v in z)


def snap(z):
    out, snapped = [], False
    for v in z:
        if EXCL[0] < v < EXCL[1]:
            v = EXCL[0] if v - EXCL[0] <= EXCL[1] - v else EXCL[1]
            snapped = True
        out.append(v)
    return sorted(out), snapped


def est_rib_mass(z):
    pts = [0.0] + list(z) + [SPAN]
    return K_MASS * sum((300.0 - 0.15 * v) ** 2 for v in pts)


def main():
    rows = []

    def add(family, z, note=""):
        z = [round(float(v), 1) for v in z]
        if not valid(z):
            raise SystemExit(f"invalid layout for {family}: {z}")
        rows.append(dict(id="L%03d" % len(rows), family=family, z1=z[0], z2=z[1], z3=z[2], z4=z[3],
                         stations=",".join("%g" % v for v in [0.0] + z + [SPAN]),
                         est_rib_mass_g=round(est_rib_mass(z), 1), note=note))

    base = [200.0, 400.0, 600.0, 800.0]
    add("uniform baseline", base)
    for p in (0.6, 0.8, 1.3, 1.6, 2.0):
        z, snapped = snap([SPAN * (i / 5.0) ** p for i in range(1, 5)])
        add("power-law p=%g" % p, z, "snapped out of the exclusion zone" if snapped else "")
    for z in ([60, 130, 200, 400], [80, 160, 240, 420], [50, 100, 150, 300], [100, 200, 300, 600],
              [100, 220, 780, 900], [150, 300, 700, 850], [300, 420, 540, 660], [400, 600, 800, 900]):
        add("handcrafted", z)
    for d in (1.0, 2.0, -2.0):
        add("baseline shifted %+g mm (mesh-scatter est.)" % d, [v + d for v in base])

    rng = random.Random(SEED)
    existing = [[r["z1"], r["z2"], r["z3"], r["z4"]] for r in rows]
    tries = 0
    while len(rows) < N_RANDOM_END:
        tries += 1
        if tries > 200000:
            raise SystemExit("could not place the random layouts")
        z = sorted(round(rng.uniform(40.0, 960.0), 1) for _ in range(4))
        if not valid(z):
            continue
        if any(sum(abs(a - b) for a, b in zip(z, e)) < 80.0 for e in existing):
            continue
        existing.append(z)
        add("random", z)

    rng2 = random.Random(SEED_ROOT)
    tries2 = 0
    while len(rows) < N_TOTAL:
        tries2 += 1
        if tries2 > 200000:
            raise SystemExit("could not place the root-biased layouts")
        four_below = len(rows) < N_RANDOM_END + 6
        if four_below:
            z = sorted(round(rng2.uniform(40.0, 400.0), 1) for _ in range(4))
            fam = "root-biased random (4 below z=400)"
        else:
            z = sorted([round(rng2.uniform(40.0, 400.0), 1) for _ in range(3)]
                       + [round(rng2.uniform(400.0, 960.0), 1)])
            fam = "root-biased random (3 below z=400)"
        if not valid(z):
            continue
        if any(sum(abs(a - b) for a, b in zip(z, e)) < 80.0 for e in existing):
            continue
        existing.append(z)
        add(fam, z)

    os.makedirs("results", exist_ok=True)
    with open("results/stage6_layouts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} layouts written to results/stage6_layouts.csv "
          f"(random placement took {tries} + {tries2} tries)\n")
    print(" id    family                                      stations                          est rib mass [g]")
    for r in rows[N_RANDOM_END - 40 + 0:] if False else rows:
        print(f" {r['id']}  {r['family']:44s}{r['stations']:34s}{r['est_rib_mass_g']:8.1f}  {r['note']}")
    masses = [r["est_rib_mass_g"] for r in rows]
    n3 = sum(1 for r in rows if sum(1 for k in ("z1", "z2", "z3", "z4") if r[k] < 400.0) >= 3)
    print(f"\nestimated rib mass: min {min(masses):.1f} g, max {max(masses):.1f} g, baseline {masses[0]:.1f} g")
    print(f"layouts with >= 3 interior ribs below z = 400: {n3} of {len(rows)}")


if __name__ == "__main__":
    main()
