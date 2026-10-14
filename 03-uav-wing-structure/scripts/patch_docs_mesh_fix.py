#!/usr/bin/env python3
"""
Guarded fix of three README sentences in section 6 (Stage 6 mesh sensitivity).
Aborts unless each anchor occurs exactly once and the dev-log source lines match.
Backup: /tmp/README.md.bak_docsfix
Run: /usr/bin/python3 scripts/patch_docs_mesh_fix.py
"""
import csv
import shutil

# --- verify the dev-log source of the Stage 3 numbers ---
L = open("development-log.md").read().split("\n")
if "buckling factors at" not in L[736] or "t = 0.5: 0.9122 -> 0.9325 (+2.2%)" not in L[738]:
    raise SystemExit("ABORT: development-log.md lines 737/739 do not match the expected Stage 3 text")

# --- numbers from the CSV ---
D = {}
for r in csv.DictReader(open("results/stage6d_mesh_sensitivity.csv")):
    D[(r["layout"], float(r["skin_t"]), float(r["skin_mesh"]))] = float(r["buckling"])


def gap(skin, k):
    return 100 * (D[("L006", skin, k)] / D[("L126", skin, k)] - 1)


g1 = [gap(1.0, k) for k in (10, 8, 6)]
g06 = [gap(0.6, k) for k in (10, 8, 6)]

R = open("README.md").read()


def edit(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit("ABORT: anchor '%s' occurs %d times (need 1)" % (label, n))
    return text.replace(old, new)


R = edit(R,
         "layout rankings (L006 exceeds L126 by 1.43%, 1.41% and 1.47% at 10, 8 and 6 mm, i.e. a constant gap)",
         "layout rankings (L006 exceeds L126 in buckling by %.2f%%, %.2f%% and %.2f%% at 10, 8 and 6 mm skin mesh with 1.0 mm skin, and by %.2f%%, %.2f%% and %.2f%% with 0.6 mm skin, i.e. an essentially constant gap)" % (tuple(g1) + tuple(g06)),
         "gap")
R = edit(R,
         "lowers buckling by -0.83% to -0.51% at 1.0 mm skin and by -1.50% to -1.48% at 0.6 mm",
         "lowers buckling by 0.51% to 0.83% at 1.0 mm skin and by 1.48% to 1.50% at 0.6 mm",
         "double negative")
R = edit(R,
         "**Earlier evidence is consistent in size but not in sign:** the Stage 3 skin-mesh check (8 → 6 mm, development log) gave buckling −0.5% at t = 1.0 mm and +2.2% at t = 0.5 mm,",
         "**The sign of the skin-mesh effect is not consistent between cases:** the Stage 3 check on its baseline layout (8 → 6 mm, development log) gave buckling −0.5% at t = 1.0 mm and +2.2% at t = 0.5 mm (0.9122 → 0.9325), whereas the Stage 6 layouts at 0.6 mm give about −1.5%,",
         "stage3 sentence")

shutil.copy("README.md", "/tmp/README.md.bak_docsfix")
open("README.md", "w").write(R)
print("FIX_OK  gaps 1.0mm: %s   0.6mm: %s" % (["%.2f" % x for x in g1], ["%.2f" % x for x in g06]))

