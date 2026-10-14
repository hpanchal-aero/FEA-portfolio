#!/usr/bin/env python3
"""Guarded patch: store skin p99 / p99.9 / design stress / solid-only max in run_layout_study.py records and CSV."""
import py_compile
import shutil
import sys

P = "scripts/run_layout_study.py"
src = open(P).read()

BLOCK = '''    sr = run([SYS_PY, "scripts/stress_metrics_cli.py", c.N["job"]], "stress-metrics", c, TIMEOUT["extract"])
    sm = None
    for line in (sr["text"] or "").splitlines():
        if line.startswith("STRESS_JSON "):
            sm = json.loads(line[len("STRESS_JSON "):])
    if sr["rc"] != 0 or sr["timed_out"] or sm is None:
        return False, {"wall_s": round(r["wall"], 1)}, "stress metrics failed rc=%s; see %s" % (sr["rc"], sr["log"])
    sm = {k: sm[k] for k in ("skin_p99", "skin_p99_9", "skin_max_dTE_ge3", "solid_only_max")}
'''

EDITS = [
    ('    c.rec["metrics"].update(tip_mm=d["uy_spar_axis_mm"], mass_kg=mass,',
     BLOCK + '    c.rec["metrics"].update(tip_mm=d["uy_spar_axis_mm"], mass_kg=mass,'),
    ('bracket_lo_mm=d["bracket_lo_mm"], bracket_hi_mm=d["bracket_hi_mm"])\n',
     'bracket_lo_mm=d["bracket_lo_mm"], bracket_hi_mm=d["bracket_hi_mm"])\n    c.rec["metrics"].update(sm)\n'),
    ('"b1", "b2", "b3", "b4", "b5", "b6", "vm_skin_max_MPa", "vm_solid_only_max_MPa"]',
     '"b1", "b2", "b3", "b4", "b5", "b6", "vm_skin_max_MPa", "vm_solid_only_max_MPa",\n              "skin_p99", "skin_p99_9", "skin_max_dTE_ge3", "solid_only_max"]'),
    ('for k in ("tip_mm", "mass_kg", "vm_skin_max_MPa", "vm_solid_only_max_MPa"):',
     'for k in ("tip_mm", "mass_kg", "vm_skin_max_MPa", "vm_solid_only_max_MPa", "skin_p99", "skin_p99_9",\n                          "skin_max_dTE_ge3", "solid_only_max"):'),
]

for old, new in EDITS:
    n = src.count(old)
    if n != 1:
        sys.exit("ABORT: anchor occurs %d times (need 1): %r" % (n, old[:70]))
shutil.copy(P, P + ".bak_stress")
for old, new in EDITS:
    src = src.replace(old, new)
open(P, "w").write(src)
py_compile.compile(P, doraise=True)
print("PATCH_OK  (backup: %s.bak_stress)" % P)
