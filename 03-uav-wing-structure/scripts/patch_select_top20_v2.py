#!/usr/bin/env python3
"""
Guarded patch (v2) for select_top20.py and select_top20_final.py.
- adds --first-pass (layouts with id number < 100, the set the tier 2 selection was made from)
- output: --first-pass -> results/stage6_fronts.csv ; otherwise results/stage6_fronts_all.csv ; tolerance run adds _tol
Aborts unless each anchor occurs exactly once. Backups in /tmp (*.bak_fp). Compiles both files.
Run: /usr/bin/python3 scripts/patch_select_top20_v2.py
"""
import py_compile
import shutil


def edit(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit("ABORT: anchor '%s' occurs %d times (need 1)" % (label, n))
    return text.replace(old, new)


P1 = "scripts/select_top20.py"
P2 = "scripts/select_top20_final.py"
S = open(P1).read()
F = open(P2).read()

S = edit(S, '    ap.add_argument("--tol-tip", type=float, default=0.0)\n',
         '    ap.add_argument("--tol-tip", type=float, default=0.0)\n'
         '    ap.add_argument("--first-pass", action="store_true",\n'
         '                    help="only layouts with id number < 100 (the set the tier 2 selection was made from)")\n', "args")
S = edit(S, '        if r["status"] == "complete":',
         '        if r["status"] == "complete" and (not a.first_pass or int(r["id"][1:]) < 100):', "filter")
S = edit(S, '    out_path = "results/stage6_fronts.csv" if (a.tol_buckling == 0.0 and a.tol_tip == 0.0) else "results/stage6_fronts_tol.csv"',
         '    out_path = "results/stage6_fronts" + ("" if a.first_pass else "_all") + \\\n'
         '               ("" if (a.tol_buckling == 0.0 and a.tol_tip == 0.0) else "_tol") + ".csv"', "out_path")
S = edit(S, "Writes results/stage6_fronts.csv.",
         "Writes results/stage6_fronts.csv with --first-pass (layouts L000-L099, the set the tier 2 selection was made from), results/stage6_fronts_all.csv without it; a tolerance run adds the suffix _tol.", "docstring")
F = edit(F, "(from select_top20.py run with default tolerances)", "(from select_top20.py --first-pass with default tolerances)", "final docstring")

shutil.copy(P1, "/tmp/select_top20.py.bak_fp")
shutil.copy(P2, "/tmp/select_top20_final.py.bak_fp")
open(P1, "w").write(S)
open(P2, "w").write(F)
py_compile.compile(P1, doraise=True)
py_compile.compile(P2, doraise=True)
print("PATCH_OK")	
