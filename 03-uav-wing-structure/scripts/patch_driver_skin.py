#!/usr/bin/env python3
"""Guarded patch: add --skin-t and --stage skin to run_layout_study.py (default 1.0 mm behaviour unchanged)."""
import py_compile
import shutil
import sys

P = "scripts/run_layout_study.py"
src = open(P).read()

EDITS = [
    # globals + helper
    ('DEFAULT_RIB_T = 1.5\n',
     'DEFAULT_RIB_T = 1.5\nSKIN_T = 1.0          # skin thickness [mm]; changed by --skin-t (stage "skin" only)\nSKIN_SUFFIX = ""     # "" at 1.0 mm, "_sk0p8" otherwise (record keys, summary names)\n\n\n'
     'def skin_p():\n    return ("%g" % SKIN_T).replace(".", "p")\n'),
    # job names
    ('    n["static_job"] = "simulation/stage3_skin_s4_k8%s_closeedgefacnone_t1_n3p8" % tag\n',
     '    n["static_job"] = "simulation/stage3_skin_s4_k8%s_closeedgefacnone_t%s_n3p8" % (tag, skin_p())\n'),
    ('    n["job"] = "simulation/stage4_s4_k8%s_closeedgefacnone_t1_n3p8" % tag\n',
     '    n["job"] = "simulation/stage4_s4_k8%s_closeedgefacnone_t%s_n3p8" % (tag, skin_p())\n'),
    ('    n["viewer"] = "mesh/stage3_skin/stage3_skin_s4_k8%s_closeedgefacnone_t1_n3p8_viewer.npz" % tag\n',
     '    n["viewer"] = "mesh/stage3_skin/stage3_skin_s4_k8%s_closeedgefacnone_t%s_n3p8_viewer.npz" % (tag, skin_p())\n'),
    # record key
    ('        self.key = self.tag.lstrip("_")\n',
     '        self.key = self.tag.lstrip("_") + SKIN_SUFFIX\n'),
    # record matching / defaults
    ('        if d.get("stations") == c.L["stations"] and d.get("rib_t") == c.rib_t:\n',
     '        if d.get("stations") == c.L["stations"] and d.get("rib_t") == c.rib_t and d.get("skin_t", 1.0) == SKIN_T:\n'),
    ('tag=c.tag, steps={}, metrics={}, failed_step=None, failed_reason=None, complete=False)',
     'tag=c.tag, skin_t=SKIN_T, steps={}, metrics={}, failed_step=None, failed_reason=None, complete=False)'),
    # input builder thickness
    ('c.N["raw"], "1.0", "372.78", "3.8"]',
     'c.N["raw"], "%g" % SKIN_T, "372.78", "3.8"]'),
    # skin step list
    ('("buckle_solve", step_buckle_solve), ("buckle_extract", step_buckle_extract), ("cleanup", step_cleanup)]\n',
     '("buckle_solve", step_buckle_solve), ("buckle_extract", step_buckle_extract), ("cleanup", step_cleanup)]\n'
     'SKIN_STEPS = PASS1[3:] + PASS2      # inp, load, then all of pass 2 (mesh/imprint/geometry are reused)\n'),
    ('    if steps is PASS2:\n        c.rec["complete"] = True\n',
     '    if steps is PASS2 or steps is SKIN_STEPS:\n        c.rec["complete"] = True\n'),
    # summary
    ('            key = tag.lstrip("_")\n',
     '            key = tag.lstrip("_") + SKIN_SUFFIX\n'),
    ('    counts = {}\n    with open(out, "w", newline="") as f:\n',
     '    if SKIN_T != 1.0:\n        out = out.replace(".csv", "_sk%s.csv" % skin_p())\n    counts = {}\n    with open(out, "w", newline="") as f:\n'),
    # main: option, stage, globals, guards
    ('choices=["check", "pass1", "pass2", "summary"])',
     'choices=["check", "pass1", "pass2", "summary", "skin"])'),
    ('    ap.add_argument("--rib-t", type=float, default=DEFAULT_RIB_T)\n',
     '    ap.add_argument("--rib-t", type=float, default=DEFAULT_RIB_T)\n    ap.add_argument("--skin-t", type=float, default=1.0)\n'),
    ('def main():\n    global LOG\n',
     'def main():\n    global LOG, SKIN_T, SKIN_SUFFIX\n'),
    ('    a = ap.parse_args()\n    os.makedirs(RES, exist_ok=True)\n',
     '    a = ap.parse_args()\n'
     '    if a.stage == "skin" and a.skin_t == 1.0:\n'
     '        sys.exit("ERROR: --stage skin needs --skin-t other than 1.0 (1.0 mm is the existing pass 1/pass 2 study).")\n'
     '    if a.stage in ("pass1", "pass2") and a.skin_t != 1.0:\n'
     '        sys.exit("ERROR: --skin-t other than 1.0 is only for --stage skin (or summary).")\n'
     '    if a.skin_t <= 0:\n        sys.exit("ERROR: --skin-t must be positive.")\n'
     '    SKIN_T = a.skin_t\n'
     '    SKIN_SUFFIX = "" if SKIN_T == 1.0 else "_sk" + skin_p()\n'
     '    os.makedirs(RES, exist_ok=True)\n'),
    ('    steps = PASS1 if a.stage == "pass1" else PASS2\n',
     '    steps = PASS1 if a.stage == "pass1" else (SKIN_STEPS if a.stage == "skin" else PASS2)\n'),
    ('            c.rec = load_rec(c)\n            if a.stage == "pass2" and not step_done(c, "load"):\n',
     '            c.rec = load_rec(c)\n'
     '            if a.stage == "skin":\n'
     '                miss = [k for k in ("raw", "json", "vol") if not os.path.exists(c.N[k])]\n'
     '                if miss:\n'
     '                    LOG.say("[%s] skin stage skipped: missing reusable files %s" % (c.key, miss))\n'
     '                    continue\n'
     '            if a.stage == "pass2" and not step_done(c, "load"):\n'),
    ('    LOG.say("=== driver start: stage=%s layouts=%d rib_t=%g retry_failed=%s ===" % (a.stage, len(sel), a.rib_t, a.retry_failed))\n',
     '    LOG.say("=== driver start: stage=%s layouts=%d rib_t=%g skin_t=%g retry_failed=%s ===" % (a.stage, len(sel), a.rib_t, SKIN_T, a.retry_failed))\n'),
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
shutil.copy(P, P + ".bak_skin")
for old, new in EDITS:
    src = src.replace(old, new)
open(P, "w").write(src)
py_compile.compile(P, doraise=True)
print("PATCH_OK  (%d edits, backup %s.bak_skin)" % (len(EDITS), P))
