"""
Project 01 - Aerospace Mounting Bracket
Stage: Strip Gmsh skin (CPS6) elements from the .inp mesh export.

Purpose: Gmsh's Abaqus/.inp export writes standalone CPS6 (plane-stress
triangle) elements for each tagged Physical Surface, in addition to the
real C3D10 volume elements. These CPS6 blocks are only useful as
node-tagging metadata during mesh inspection - CalculiX interprets any
*ELEMENT block as a real element to analyze, and plane-stress elements
are only valid in the z=0 plane, causing a fatal error for any skin face
not on that plane.

This script keeps ONLY:
  - the *NODE block (all nodes, needed by both skin and volume elements)
  - the *ELEMENT, TYPE=C3D10 block (the actual volume mesh)
and discards all CPS6 *ELEMENT blocks. Node and element numbering is
unchanged, so this remains fully consistent with the previously
generated *NSET file and the facial-surface file.
"""

IN_FILE = "../mesh/verification_plate.inp"
OUT_FILE = "../mesh/verification_plate_clean.inp"

with open(IN_FILE) as f:
    lines = f.readlines()

out_lines = []
skip = False
kept_blocks = 0
dropped_blocks = 0

i = 0
while i < len(lines):
    line = lines[i]
    upper = line.strip().upper().replace(" ", "")

    if upper.startswith("*ELEMENT,TYPE=CPS6"):
        skip = True
        dropped_blocks += 1
        i += 1
        continue

    if upper.startswith("*ELEMENT,TYPE=C3D10"):
        skip = False
        kept_blocks += 1
        out_lines.append(line)
        i += 1
        continue

    if upper.startswith("*") and skip:
        # any new keyword ends the CPS6 block we were skipping
        skip = False

    if not skip:
        out_lines.append(line)

    i += 1

with open(OUT_FILE, "w") as f:
    f.writelines(out_lines)

print(f"Kept C3D10 blocks: {kept_blocks}")
print(f"Dropped CPS6 blocks: {dropped_blocks}")
print(f"Wrote: {OUT_FILE}")
