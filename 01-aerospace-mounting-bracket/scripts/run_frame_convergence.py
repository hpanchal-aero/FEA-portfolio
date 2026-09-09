"""
Stage 4a convergence sweep driver. The already-solved baseline
(CORNER_SIZE=0.5mm) becomes "level_medium" for consistency with
Stage 3's naming - not regenerated here, only referenced. Generates
a coarser and a finer level around it.
"""

import sys
sys.path.insert(0, ".")
from frame_convergence_lib import write_and_mesh

BASE_DIR = "01-aerospace-mounting-bracket/mesh/frame_study/convergence"

LEVELS = {
    "level_coarse": 1.0,
    "level_fine": 0.25,
}

for level_name, corner_size in LEVELS.items():
    case_dir = f"{BASE_DIR}/{level_name}"
    write_and_mesh(case_dir, corner_size)
