#!/usr/bin/env python3
"""Guarded patch: update the Status block of README.md and the Phase B line of the development-log Status Log."""
import py_compile  # noqa: F401  (kept for symmetry with other patch scripts)
import shutil
import sys

EDITS = {
    "README.md": [(
        "Stage 6 (rib placement) not started; the project is not complete.",
        "Stage 6 (rib layout, rib count, rib thickness and skin-gauge sensitivity) is complete as a\n"
        "computational study: verification only, no validation, and no finer-mesh convergence check on\n"
        "the Stage 6 layouts yet. The project is not yet marked complete (open items are listed under\n"
        "Limitations of Stage 6)."
    )],
    "development-log.md": [(
        "- Phase B (rib placement search): not started\n",
        "- Stage 6 (rib placement, rib count, rib thickness, skin gauge): runs, analysis and documentation complete\n"
        "  (this entry, section \"Stage 6\"); 121 layouts solved at 1.5 mm ribs, 63 thin-rib runs, 12 skin-gauge runs;\n"
        "  finer-mesh check of the top layouts not run; local commit to follow, nothing pushed\n"
    )],
}

for path, edits in EDITS.items():
    src = open(path).read()
    for old, new in edits:
        n = src.count(old)
        if n != 1:
            sys.exit("ABORT: %s: anchor occurs %d times (need 1): %r" % (path, n, old[:60]))
for path, edits in EDITS.items():
    src = open(path).read()
    shutil.copy(path, path + ".bak_status6")
    for old, new in edits:
        src = src.replace(old, new)
    open(path, "w").write(src)
    print("PATCHED %s (backup %s.bak_status6)" % (path, path))
