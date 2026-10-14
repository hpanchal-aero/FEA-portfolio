#!/usr/bin/env python3
"""
Guarded patch for scripts/render_stage6_3d.py:
 - normalise the buckling mode by the maximum displacement NORM (colour bar then tops out at 1.00)
 - narrower view angle (20 deg) for the full-wing stills, ramp and buckling overview
 - colour-bar layout (height, fonts) shared by all three bars
Aborts unless each anchor occurs exactly once. Backup /tmp/render_stage6_3d.py.bak_layout. Compiles the result.
Run: /usr/bin/python3 scripts/patch_render_layout.py
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


S = edit(S, "ZOOM_OFFSET = np.array([-120.0, 380.0, 260.0])",
         "SB = dict(width=0.5, height=0.08, position_x=0.25, position_y=0.03, title_font_size=11, label_font_size=10)\n"
         "ZOOM_OFFSET = np.array([-120.0, 380.0, 260.0])", "SB definition")
S = edit(S, 'scalar_bar_args=dict(title="von Mises, solid parts [MPa] (capped at %.0f, same scale for all cases)" % cap))',
         'scalar_bar_args=dict(SB, title="von Mises, solid parts [MPa] (capped at %.0f, same scale for all cases)" % cap))', "stress bar")
S = edit(S, 'scalar_bar_args=dict(title="Uy [mm] (limit load)"))',
         'scalar_bar_args=dict(SB, title="Uy [mm] (limit load)"))', "uy bar")
S = edit(S, 'scalar_bar_args=dict(title="|U| / max (buckling mode 1)"))',
         'scalar_bar_args=dict(SB, title="|U| / max|U| (buckling mode 1)"))', "buckling bar")
S = edit(S, "        p.camera_position = CAM\n",
         "        p.camera_position = CAM\n        p.camera.view_angle = 20.0\n", "scene camera")
S = edit(S, "        p.camera_position = cam\n",
         "        p.camera_position = cam\n        if cam is CAM:\n            p.camera.view_angle = 20.0\n", "buckling camera")
S = edit(S, "    amp = np.abs(u).max()", "    amp = np.linalg.norm(u, axis=1).max()", "amp")
shutil.copy(P, "/tmp/render_stage6_3d.py.bak_layout")
open(P, "w").write(S)
py_compile.compile(P, doraise=True)
print("PATCH_OK")
