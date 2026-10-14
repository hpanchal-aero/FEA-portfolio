#!/usr/bin/env python3
"""
Guarded patch for scripts/render_stage6_3d.py:
 - caption note on trailing-edge rib-tip peaks above the colour scale
 - print solid-vertex max/p99 excluding the TE zone (dTE < 3 mm) and root (z <= 10 mm)
Aborts unless each anchor occurs exactly once. Backup /tmp/render_stage6_3d.py.bak_cap. Compiles the result.
Run: /usr/bin/python3 scripts/patch_render_caption.py
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


S = edit(S, '        p.add_text("%s\\n%gx deformation, %d%% of limit load; placeholder load; skin stress not shown"\n',
         '        p.add_text("%s\\n%gx deformation, %d%% of limit load; placeholder load; skin stress not shown\\n'
         'trailing-edge rib-tip peaks exceed the scale (mesh-dependent singularity, excluded from design stress)"\n', "caption")
S = edit(S, "    sk = np.where(is_solid, vm, np.nan)\n",
         "    dte = 0.7 * (300.0 - 0.15 * P[:, 2]) - P[:, 0]\n"
         "    ok = is_solid & (dte >= 3.0) & (P[:, 2] > 10.0)\n"
         "    print(\"[%s] solid-vertex stress excluding the TE zone (dTE < 3 mm) and root (z <= 10 mm): %d vertices, max %.1f MPa, p99 %.1f MPa\"\n"
         "          % (lid, ok.sum(), vm[ok].max(), np.percentile(vm[ok], 99)))\n"
         "    sk = np.where(is_solid, vm, np.nan)\n", "stress print")
shutil.copy(P, "/tmp/render_stage6_3d.py.bak_cap")
open(P, "w").write(S)
py_compile.compile(P, doraise=True)
print("PATCH_OK")
EOF_MARKER = None
