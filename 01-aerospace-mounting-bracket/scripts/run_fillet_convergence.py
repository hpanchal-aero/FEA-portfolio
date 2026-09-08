"""
Stage 3 convergence sweep driver. Generates mesh + NSETs for the
coarse (0.30mm) and fine (0.075mm) levels only - the medium (0.15mm)
level is the already-solved baseline and is not regenerated here.

Does NOT run the quality check or the solver - those remain separate,
inspectable steps per the established protocol.
"""

import sys
sys.path.insert(0, ".")
from fillet_convergence_lib import write_and_mesh

BASE_DIR = "01-aerospace-mounting-bracket/mesh/fillet_study/convergence"

LEVELS = {
    "level_intermediate2": 0.06,
}

for level_name, fillet_size in LEVELS.items():
    case_dir = f"{BASE_DIR}/{level_name}"
    write_and_mesh(case_dir, fillet_size)
