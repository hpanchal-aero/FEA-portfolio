#!/usr/bin/env python3
"""
Guarded patch for scripts/render_stage6_3d.py: colour-bar titles drawn as plain text below the bar
(the scalar-bar widget's own title overlaps its tick labels). Aborts unless each anchor occurs exactly once.
Backup /tmp/render_stage6_3d.py.bak_barlabel. Compiles the result.
Run: /usr/bin/python3 scripts/patch_render_bar_label.py
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


S = edit(S, "def load_case(lid):",
         "def bar_label(p, text):\n"
         "    p.add_text(text, position=(int(0.25 * p.window_size[0]), 6), font_size=10, color=\"black\")\n\n\n"
         "def load_case(lid):", "bar_label definition")
S = edit(S, 'scalar_bar_args=dict(SB, title="von Mises, solid parts [MPa] (capped at %.0f, same scale for all cases)" % cap))',
         'scalar_bar_args=dict(SB, title=" "))\n'
         '            bar_label(p, "von Mises, solid parts [MPa] (capped at %.0f, same scale for all cases)" % cap)', "stress bar")
S = edit(S, 'scalar_bar_args=dict(SB, title="Uy [mm] (limit load)"))',
         'scalar_bar_args=dict(SB, title=" "))\n'
         '            bar_label(p, "Uy [mm] (limit load)")', "uy bar")
S = edit(S, 'scalar_bar_args=dict(SB, title="|U| / max|U| (buckling mode 1)"))',
         'scalar_bar_args=dict(SB, title=" "))\n'
         '        bar_label(p, "|U| / max|U| (buckling mode 1)")', "buckling bar")
shutil.copy(P, "/tmp/render_stage6_3d.py.bak_barlabel")
open(P, "w").write(S)
py_compile.compile(P, doraise=True)
print("PATCH_OK")
