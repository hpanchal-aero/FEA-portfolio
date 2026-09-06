"""
Project 01 - Aerospace Mounting Bracket
Stage 2: Hole-region mesh convergence library.

Reuses the exact locked geometry/material/BC/load setup - only
hole_size (local mesh refinement at the hole) varies between levels.
far_size and transition_width are held fixed at their validated
values (3.0mm, 1.5mm) so only the hole-region resolution changes.
"""

import os
import re
import subprocess
import sys

sys.path.insert(0, ".")
import hex_lib as lib

L, B, H = 60.0, 40.0, 4.0
F = 235.44
HOLE_X, HOLE_Y, HOLE_D = 30.0, 20.0, 5.0
HOLE_R = HOLE_D / 2.0
FAR_SIZE = 3.0
TRANSITION_WIDTH = 1.5

I_SECTION = B * H**3 / 12.0
C = H / 2.0


def sigma_theory(x):
    M = F * (L - x)
    return M * C / I_SECTION


def write_geo(case_dir, hole_size):
    geo_path = os.path.join(case_dir, "case_hole.geo")
    content = f"""// Hole-region convergence case: hole_size = {hole_size} mm
// far_size and transition_width held fixed at locked values.

SetFactory("OpenCASCADE");
L = {L}; B = {B}; H = {H};
hole_x = {HOLE_X}; hole_y = {HOLE_Y}; hole_r = {HOLE_R};

Box(1) = {{0, 0, 0, L, B, H}};
Cylinder(2) = {{hole_x, hole_y, -1.0, 0, 0, H+2.0, hole_r}};
BooleanDifference(3) = {{ Volume{{1}}; Delete; }}{{ Volume{{2}}; Delete; }};

eps = 1e-3;
root_surf()      = Surface In BoundingBox{{-eps,-eps,-eps, eps,B+eps,H+eps}};
tip_surf()       = Surface In BoundingBox{{L-eps,-eps,-eps, L+eps,B+eps,H+eps}};
bottom_surf()    = Surface In BoundingBox{{-eps,-eps,-eps, L+eps,B+eps,eps}};
top_surf()       = Surface In BoundingBox{{-eps,-eps,H-eps, L+eps,B+eps,H+eps}};
side_y0_surf()   = Surface In BoundingBox{{-eps,-eps,-eps, L+eps,eps,H+eps}};
side_y40_surf()  = Surface In BoundingBox{{-eps,B-eps,-eps, L+eps,B+eps,H+eps}};
hole_surf()      = Surface In BoundingBox{{hole_x-hole_r-eps, hole_y-hole_r-eps, -eps,
                                            hole_x+hole_r+eps, hole_y+hole_r+eps, H+eps}};

Physical Surface("root", 101)     = root_surf();
Physical Surface("tip", 102)      = tip_surf();
Physical Surface("bottom", 103)   = bottom_surf();
Physical Surface("top", 104)      = top_surf();
Physical Surface("side_y0", 105)  = side_y0_surf();
Physical Surface("side_y40", 106) = side_y40_surf();
Physical Surface("hole", 107)     = hole_surf();
Physical Volume("plate", 201)     = {{3}};

far_size = {FAR_SIZE};
hole_size = {hole_size};
transition_width = {TRANSITION_WIDTH};

Field[1] = Distance;
Field[1].SurfacesList = {{hole_surf()}};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hole_size;
Field[2].SizeMax = far_size;
Field[2].DistMin = hole_r;
Field[2].DistMax = hole_r + transition_width;
Background Field = 2;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthMax = far_size;
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 1;
Mesh.Algorithm3D = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh 3;
"""
    with open(geo_path, "w") as f:
        f.write(content)
    return geo_path


def mesh_only(case_dir):
    """Run gmsh CLI (not Python API) to generate the mesh, matching
    the exact command sequence used and verified in prior stages.

    Uses absolute paths throughout - passing relative paths together
    with cwd= to subprocess.run is a double-relative-path bug (paths
    get resolved against the new cwd, not the caller's cwd) - see
    chat log for the failure this caused.
    """
    case_dir = os.path.abspath(case_dir)
    geo_path = os.path.join(case_dir, "case_hole.geo")
    raw_inp = os.path.join(case_dir, "mesh_hole_raw.inp")
    result = subprocess.run(
        ["/usr/bin/gmsh", geo_path, "-3", "-o", raw_inp, "-format", "inp"],
        capture_output=True, text=True, cwd=case_dir
    )
    log_path = os.path.join(case_dir, "gmsh_log.txt")
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)
    # Extract final element count from log
    m = re.findall(r"(\d+) nodes (\d+) elements", result.stdout)
    n_nodes, n_elements = (int(m[-1][0]), int(m[-1][1])) if m else (None, None)
    return raw_inp, n_nodes, n_elements
