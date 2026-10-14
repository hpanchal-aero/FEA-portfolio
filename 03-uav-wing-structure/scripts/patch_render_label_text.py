#!/usr/bin/env python3
"""
Guarded patch for scripts/render_stage6_3d.py: the buckling colour-bar label without '|' characters
(the plotting font does not draw them). Aborts unless the anchor occurs exactly once.
Backup /tmp/render_stage6_3d.py.bak_label. Compiles the result.
Run: /usr/bin/python3 scripts/patch_render_label_text.py
"""
import py_compile
import shutil

P = "scripts/render_stage6_3d.py"
S = open(P).read()
old = 'bar_label(p, "|U| / max|U| (buckling mode 1)")'
if S.count(old) != 1:
    raise SystemExit("ABORT: anchor occurs %d times (need 1)" % S.count(old))
S = S.replace(old, 'bar_label(p, "displacement magnitude / its maximum (buckling mode 1)")')
shutil.copy(P, "/tmp/render_stage6_3d.py.bak_label")
open(P, "w").write(S)
py_compile.compile(P, doraise=True)
print("PATCH_OK")
