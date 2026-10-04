#!/usr/bin/env python3
"""
Convert a static skin .inp (from build_stage3_skin_inp.py) into a *BUCKLE input, mirroring the
syntax verified in the plate / panel buckling checks: boundary conditions in the model data,
then a single '*BUCKLE / <n modes>' step containing the loads (CalculiX solves the preload
itself). The buckling factor multiplies the WHOLE applied load (lift + gravity x n).

Usage: python3 scripts/make_buckle_inp.py <static.inp> <out_buckle.inp> [n_modes=6]
"""
import sys

src, out = sys.argv[1], sys.argv[2]
nmodes = int(sys.argv[3]) if len(sys.argv) > 3 else 6
txt = open(src).read()
marker = "*STEP\n*STATIC\n*BOUNDARY\nNROOT, 1, 6, 0.0\n"
if txt.count(marker) != 1:
    sys.exit("ERROR: expected step header not found exactly once -- nothing written.")
model, rest = txt.split(marker)
if "*NODE FILE" not in rest:
    sys.exit("ERROR: no *NODE FILE in step -- nothing written.")
loads = rest.split("*NODE FILE")[0]
if "*DLOAD" not in loads or "*CLOAD" not in loads:
    sys.exit("ERROR: loads block missing *DLOAD or *CLOAD -- nothing written.")
new = (model + "*BOUNDARY\nNROOT, 1, 6, 0.0\n"
       + f"*STEP\n*BUCKLE\n{nmodes}\n" + loads + "*NODE FILE\nU\n*END STEP\n")
open(out, "w").write(new)
print(f"wrote {out}: {len(new.splitlines())} lines (static input had {len(txt.splitlines())}); "
      f"{nmodes} modes; loads block {len(loads.splitlines())} lines")
