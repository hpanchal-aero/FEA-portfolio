"""
Targeted toe-convergence sweep driver. Baseline (toe implicitly
~0.5mm via the hypotenuse-face field alone) already solved - this
generates two further-refined levels for direct comparison.
"""

import sys
sys.path.insert(0, ".")
from frame_gusset_convergence_lib import write_and_mesh

BASE_DIR = "01-aerospace-mounting-bracket/mesh/frame_gusset_study/toe_convergence"

LEVELS = {
    "level_toe_intermediate": 0.10,
}

for level_name, toe_size in LEVELS.items():
    case_dir = f"{BASE_DIR}/{level_name}"
    write_and_mesh(case_dir, toe_size)
