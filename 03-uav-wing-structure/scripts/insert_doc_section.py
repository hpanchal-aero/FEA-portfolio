#!/usr/bin/env python3
"""
Guarded insertion of a markdown fragment before an anchor heading.
Usage: insert_doc_section.py <fragment.md> <target.md> "<anchor heading line>"
Aborts unless: the anchor occurs exactly once as a full line, the fragment's first line is a '## ' heading, and that
heading does not already exist in the target. Writes a .bak copy, inserts, prints the heading list.
"""
import re
import shutil
import sys

frag_p, tgt_p, anchor = sys.argv[1], sys.argv[2], sys.argv[3]
frag = open(frag_p).read().rstrip("\n") + "\n\n"
lines = open(tgt_p).read().split("\n")
first = frag.split("\n")[0]
if not first.startswith("## "):
    sys.exit("ABORT: fragment must start with a '## ' heading, got %r" % first)
if first in lines:
    sys.exit("ABORT: heading already present in target: %r" % first)
idx = [i for i, l in enumerate(lines) if l == anchor]
if len(idx) != 1:
    sys.exit("ABORT: anchor %r occurs %d times (need 1)" % (anchor, len(idx)))
shutil.copy(tgt_p, tgt_p + ".bak_stage6")
new = lines[:idx[0]] + frag.split("\n") + lines[idx[0]:]
open(tgt_p, "w").write("\n".join(new))
print("INSERTED before line %d of %s (backup %s.bak_stage6)" % (idx[0] + 1, tgt_p, tgt_p))
for n, l in enumerate(open(tgt_p).read().split("\n"), 1):
    if re.match(r"^#{1,3} ", l):
        print("%5d  %s" % (n, l))
