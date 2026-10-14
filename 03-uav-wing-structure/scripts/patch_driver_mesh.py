#!/usr/bin/env python3
"""Guarded patch: add --mesh SOLID,SKIN and --stage conv to run_layout_study.py (default 4,8 mm behaviour unchanged)."""
import py_compile
import shutil
import sys

P = "scripts/run_layout_study.py"
src = open(P).read()

CONV_GATE = '''            if a.stage == "conv":
                bp = rec_path(c.tag.lstrip("_"))
                if not os.path.exists(bp):
                    LOG.say("[%s] conv stage skipped: base record %s is missing" % (c.key, bp))
                    continue
                base = json.load(open(bp))
                if not base.get("complete") or "imprint" not in base.get("steps", {}):
                    LOG.say("[%s] conv stage skipped: base record is not complete" % c.key)
                    continue
                c.rec["steps"].setdefault("imprint", base["steps"]["imprint"])
                miss = [k for k in ("step", "json", "vol") if not os.path.exists(c.N[k])]
                if miss:
                    LOG.say("[%s] conv stage skipped: missing reusable files %s" % (c.key, miss))
                    continue
'''

EDITS = [
    # helper + globals
    ('def skin_p():\n    return ("%g" % SKIN_T).replace(".", "p")\n',
     'def skin_p():\n    return ("%g" % SKIN_T).replace(".", "p")\n\n\n'
     'MESH_S, MESH_K = 4.0, 8.0     # solid / skin maxh [mm]; changed by --mesh (stage "conv" only)\n'
     'MESH_SUFFIX = ""              # "" at 4,8 mm, "_s4_k6" otherwise (record keys, summary names)\n\n\n'
     'def mesh_p():\n    return ("s%g_k%g" % (MESH_S, MESH_K)).replace(".", "p")\n'),
    # names
    ('    n["raw"] = "mesh/stage3_skin/raw_s4_k8%s_closeedgefacnone.npz" % tag\n',
     '    n["raw"] = "mesh/stage3_skin/raw_%s%s_closeedgefacnone.npz" % (mesh_p(), tag)\n'),
    ('    n["static_job"] = "simulation/stage3_skin_s4_k8%s_closeedgefacnone_t%s_n3p8" % (tag, skin_p())\n',
     '    n["static_job"] = "simulation/stage3_skin_%s%s_closeedgefacnone_t%s_n3p8" % (mesh_p(), tag, skin_p())\n'),
    ('    n["job"] = "simulation/stage4_s4_k8%s_closeedgefacnone_t%s_n3p8" % (tag, skin_p())\n',
     '    n["job"] = "simulation/stage4_%s%s_closeedgefacnone_t%s_n3p8" % (mesh_p(), tag, skin_p())\n'),
    ('    n["viewer"] = "mesh/stage3_skin/stage3_skin_s4_k8%s_closeedgefacnone_t%s_n3p8_viewer.npz" % (tag, skin_p())\n',
     '    n["viewer"] = "mesh/stage3_skin/stage3_skin_%s%s_closeedgefacnone_t%s_n3p8_viewer.npz" % (mesh_p(), tag, skin_p())\n'),
    # mesh step sizes
    ('"scripts/mesh_stage3_skin_raw.py", "4.0", "8.0"]',
     '"scripts/mesh_stage3_skin_raw.py", "%g" % MESH_S, "%g" % MESH_K]'),
    # record key
    ('        self.key = self.tag.lstrip("_") + SKIN_SUFFIX\n',
     '        self.key = self.tag.lstrip("_") + SKIN_SUFFIX + MESH_SUFFIX\n'),
    # step lists and completion flag
    ('SKIN_STEPS = PASS1[3:] + PASS2',
     'CONV_STEPS = PASS1[2:] + PASS2      # mesh, inp, load, then all of pass 2 (geometry and imprint are reused)\nSKIN_STEPS = PASS1[3:] + PASS2'),
    ('    if steps is PASS2 or steps is SKIN_STEPS:\n',
     '    if steps is PASS2 or steps is SKIN_STEPS or steps is CONV_STEPS:\n'),
    # summary
    ('            key = tag.lstrip("_") + SKIN_SUFFIX\n',
     '            key = tag.lstrip("_") + SKIN_SUFFIX + MESH_SUFFIX\n'),
    ('    if SKIN_T != 1.0:\n        out = out.replace(".csv", "_sk%s.csv" % skin_p())\n',
     '    if SKIN_T != 1.0:\n        out = out.replace(".csv", "_sk%s.csv" % skin_p())\n'
     '    if MESH_SUFFIX:\n        out = out.replace(".csv", "%s.csv" % MESH_SUFFIX)\n'),
    # main
    ('choices=["check", "pass1", "pass2", "summary", "skin"])',
     'choices=["check", "pass1", "pass2", "summary", "skin", "conv"])'),
    ('    ap.add_argument("--skin-t", type=float, default=1.0)\n',
     '    ap.add_argument("--skin-t", type=float, default=1.0)\n    ap.add_argument("--mesh", default="4,8", help="solid,skin maxh [mm]; only for --stage conv (or summary)")\n'),
    ('    global LOG, SKIN_T, SKIN_SUFFIX\n',
     '    global LOG, SKIN_T, SKIN_SUFFIX, MESH_S, MESH_K, MESH_SUFFIX\n'),
    ('    SKIN_SUFFIX = "" if SKIN_T == 1.0 else "_sk" + skin_p()\n',
     '    SKIN_SUFFIX = "" if SKIN_T == 1.0 else "_sk" + skin_p()\n'
     '    try:\n        MESH_S, MESH_K = [float(x) for x in a.mesh.split(",")]\n'
     '    except ValueError:\n        sys.exit("ERROR: --mesh must look like 4,6")\n'
     '    if (MESH_S, MESH_K) != (4.0, 8.0) and a.stage in ("pass1", "pass2", "skin"):\n'
     '        sys.exit("ERROR: --mesh other than 4,8 is only for --stage conv (or summary).")\n'
     '    if a.stage == "conv" and ((MESH_S, MESH_K) == (4.0, 8.0) or SKIN_T != 1.0):\n'
     '        sys.exit("ERROR: --stage conv needs --mesh other than 4,8 and --skin-t 1.0.")\n'
     '    MESH_SUFFIX = "" if (MESH_S, MESH_K) == (4.0, 8.0) else "_" + mesh_p()\n'),
    ('    steps = PASS1 if a.stage == "pass1" else (SKIN_STEPS if a.stage == "skin" else PASS2)\n',
     '    steps = (PASS1 if a.stage == "pass1" else SKIN_STEPS if a.stage == "skin"\n'
     '             else CONV_STEPS if a.stage == "conv" else PASS2)\n'),
    ('            if a.stage == "pass2" and not step_done(c, "load"):\n',
     CONV_GATE + '            if a.stage == "pass2" and not step_done(c, "load"):\n'),
]

bad = []
for old, new in EDITS:
    n = src.count(old)
    if n != 1:
        bad.append((n, old[:80]))
if bad:
    for n, o in bad:
        print("ANCHOR COUNT %d (need 1): %r" % (n, o))
    sys.exit("ABORT: nothing modified.")
shutil.copy(P, P + ".bak_mesh")
for old, new in EDITS:
    src = src.replace(old, new)
open(P, "w").write(src)
py_compile.compile(P, doraise=True)
print("PATCH_OK  (%d edits, backup %s.bak_mesh)" % (len(EDITS), P))
