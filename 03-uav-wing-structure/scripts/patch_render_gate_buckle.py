#!/usr/bin/env python3
"""
Guarded patch for scripts/render_stage6_3d.py: buckling-mode stills only with RENDER_BUCKLE=1, with the reason in the docstring.
Aborts unless each anchor occurs exactly once. Backup /tmp/render_stage6_3d.py.bak_gate. Compiles the result.
Run: /usr/bin/python3 scripts/patch_render_gate_buckle.py
"""
import py_compile
import shutil

P = "scripts/render_stage6_3d.py"
S = open(P).read()


def edit(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit("ABORT: anchor '%s' occurs %d times (need 1)" % (label, n))
    return text.replace(old, new)


S = edit(S, 'maximum modal displacement (L006_buckling_mode1_zoom.png).',
         'maximum modal displacement (L006_buckling_mode1_zoom.png), only with RENDER_BUCKLE=1. These two stills are\n'
         'NOT published: their peak location (TE rib corner, z = 201 mm) disagrees with modal_span_distribution.py\n'
         '(x = 85, z = 260 mm), and the cause (probably the skin vertex-to-node mapping in render_buckle) is unresolved.', "docstring")
S = edit(S, '    if "L006" in cases:\n        render_buckle()\n',
         '    if "L006" in cases and os.environ.get("RENDER_BUCKLE") == "1":\n        render_buckle()\n', "main gate")
shutil.copy(P, "/tmp/render_stage6_3d.py.bak_gate")
open(P, "w").write(S)
py_compile.compile(P, doraise=True)
print("PATCH_OK")
