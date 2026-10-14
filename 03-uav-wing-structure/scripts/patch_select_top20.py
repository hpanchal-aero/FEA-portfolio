#!/usr/bin/env python3
"""
Guarded patch: select_top20.py writes the tolerance variant to results/stage6_fronts_tol.csv,
leaving results/stage6_fronts.csv for the no-tolerance run only. Prints the path written.
Aborts unless each anchor occurs exactly once. Backup: /tmp/select_top20.py.bak_tol
Run: /usr/bin/python3 scripts/patch_select_top20.py
"""
import py_compile
import shutil

P = "scripts/select_top20.py"
S = open(P).read()


def edit(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit("ABORT: anchor '%s' occurs %d times (need 1)" % (label, n))
    return text.replace(old, new)


S = edit(S, '    with open("results/stage6_fronts.csv", "w", newline="") as f:',
         '    out_path = "results/stage6_fronts.csv" if (a.tol_buckling == 0.0 and a.tol_tip == 0.0) else "results/stage6_fronts_tol.csv"\n'
         '    with open(out_path, "w", newline="") as f:', "open line")
S = edit(S, "% (len(rows), a.tol_buckling, a.tol_tip))",
         "% (len(rows), a.tol_buckling, a.tol_tip))\n    print(\"wrote \" + out_path)", "print line")
shutil.copy(P, "/tmp/select_top20.py.bak_tol")
open(P, "w").write(S)
py_compile.compile(P, doraise=True)
print("PATCH_OK")
