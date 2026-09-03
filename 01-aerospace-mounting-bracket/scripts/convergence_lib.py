"""
Project 01 - Aerospace Mounting Bracket
Shared library for the mesh convergence study. Every function here
reuses logic already hand-verified in the single baseline case
(generate_calculix_mesh.py, clean_mesh_for_ccx.py,
build_facial_surface.py, check_frd_displacement.py,
check_frd_stress.py) - the convergence driver does not reinvent
the pipeline, only parametrizes and loops it.
"""

import os
import re
import subprocess
import time

import gmsh

# ---- Fixed geometry (unchanged across all mesh densities) ----
L, B, H = 60.0, 40.0, 4.0
LOAD_Z = -235.44  # N, tip point load (see analysis deck)

# ---- Root extraction band (explicit, geometry-based, avoids the
# fixed-edge/free-surface singularity at x=0 - see chat discussion) ----
ROOT_X_MIN = 3.0
ROOT_X_MAX = 9.0

STL_PATH = os.path.abspath(
    "../geometry/exported/verification_plate.stl"
)

def write_geo(case_dir, elem_size):
    """Write a .geo file using Gmsh's OpenCASCADE kernel to construct
    the verification plate as an exact parametric box, rather than
    importing the OpenSCAD STL. This eliminates STL-faceting/discrete-
    surface-reconstruction artifacts that caused nonpositive-Jacobian
    tets near the box corners at finer mesh sizes (see chat log,
    mesh_0.8mm STL-based failure investigation).

    The .scad file remains the authoritative design/reference record
    for dimensions and design intent; this OCC box is used only for
    mesh generation of this simple verification geometry. Face
    identification uses Gmsh's deterministic bounding-box selector,
    not assumed surface numbering.
    """
    min_size = round(elem_size * 0.7, 4)
    geo_path = os.path.join(case_dir, "case.geo")
    eps = 1e-3
    content = f"""// Project 01 - Aerospace Mounting Bracket
// Convergence study case: element size = {elem_size} mm
// Geometry source: Gmsh OpenCASCADE exact box (NOT the OpenSCAD STL -
// see chat log for rationale: STL faceting caused nonpositive-Jacobian
// tets near box corners at fine mesh sizes).
// Reference design record: geometry/cad/01_verification_plate.scad
// (identical dimensions: L=60, B=40, H=4 mm)

SetFactory("OpenCASCADE");
Box(1) = {{0, 0, 0, {L}, {B}, {H}}};

eps = {eps};
root_surf()      = Surface In BoundingBox{{-eps,-eps,-eps, eps,{B}+eps,{H}+eps}};
tip_surf()       = Surface In BoundingBox{{{L}-eps,-eps,-eps, {L}+eps,{B}+eps,{H}+eps}};
bottom_surf()    = Surface In BoundingBox{{-eps,-eps,-eps, {L}+eps,{B}+eps,eps}};
top_surf()       = Surface In BoundingBox{{-eps,-eps,{H}-eps, {L}+eps,{B}+eps,{H}+eps}};
side_y0_surf()   = Surface In BoundingBox{{-eps,-eps,-eps, {L}+eps,eps,{H}+eps}};
side_y40_surf()  = Surface In BoundingBox{{-eps,{B}-eps,-eps, {L}+eps,{B}+eps,{H}+eps}};

Physical Surface("root", 101)     = root_surf();
Physical Surface("tip", 102)      = tip_surf();
Physical Surface("bottom", 103)   = bottom_surf();
Physical Surface("top", 104)      = top_surf();
Physical Surface("side_y0", 105)  = side_y0_surf();
Physical Surface("side_y40", 106) = side_y40_surf();
Physical Volume("plate", 201)     = {{1}};

Mesh.CharacteristicLengthMax = {elem_size};
Mesh.CharacteristicLengthMin = {min_size};
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

def mesh_and_export(case_dir, geo_path):
    """Runs the .geo through Gmsh's Python API in a single session,
    exports the .inp mesh and the NSET file from the SAME in-memory
    model (see generate_calculix_mesh.py rationale: guarantees
    consistent node numbering)."""
    raw_inp = os.path.join(case_dir, "mesh_raw.inp")
    nset_inp = os.path.join(case_dir, "nsets.inp")

    gmsh.initialize()
    gmsh.open(geo_path)
    gmsh.write(raw_inp)

    named_surfaces = ["root", "tip", "bottom", "top", "side_y0", "side_y40"]
    groups = gmsh.model.getPhysicalGroups(dim=2)
    name_to_tag = {gmsh.model.getPhysicalName(d, t): t for d, t in groups}

    with open(nset_inp, "w") as f:
        for name in named_surfaces:
            tag = name_to_tag[name]
            node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
            node_tags = sorted(int(n) for n in node_tags)
            f.write(f"*NSET,NSET={name}\n")
            for i in range(0, len(node_tags), 8):
                chunk = node_tags[i:i + 8]
                f.write(",".join(str(n) for n in chunk) + ",\n")

    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=3)
    n_elements = sum(len(t) for t in elem_tags)

    gmsh.finalize()
    return raw_inp, nset_inp, n_elements


def clean_mesh(case_dir, raw_inp):
    """Strip Gmsh's CPS6 skin elements, keep only *NODE and
    *ELEMENT,TYPE=C3D10 - see clean_mesh_for_ccx.py rationale."""
    clean_inp = os.path.join(case_dir, "mesh_clean.inp")
    with open(raw_inp) as f:
        lines = f.readlines()

    out_lines = []
    skip = False
    for line in lines:
        upper = line.strip().upper().replace(" ", "")
        if upper.startswith("*ELEMENT,TYPE=CPS6"):
            skip = True
            continue
        if upper.startswith("*ELEMENT,TYPE=C3D10"):
            skip = False
            out_lines.append(line)
            continue
        if upper.startswith("*") and skip:
            skip = False
        if not skip:
            out_lines.append(line)

    with open(clean_inp, "w") as f:
        f.writelines(out_lines)
    return clean_inp


def parse_clean_mesh(clean_inp):
    """Parse nodes and C3D10 elements from the cleaned .inp -
    identical parsing logic to build_facial_surface.py."""
    with open(clean_inp) as f:
        lines = f.readlines()

    nodes = {}
    i = 0
    while i < len(lines):
        if lines[i].strip().upper().startswith("*NODE") and \
           not lines[i].strip().upper().startswith("*NODE,"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = lines[i].strip().split(",")
                if len(parts) >= 4 and parts[0].strip().isdigit():
                    nid = int(parts[0])
                    nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
                i += 1
        else:
            i += 1

    elements = {}
    in_block = False
    buffer = []
    for line in lines:
        s = line.strip()
        if s.upper().startswith("*ELEMENT") and "TYPE=C3D10" in s.upper().replace(" ", ""):
            in_block = True
            continue
        if in_block:
            if s.startswith("*"):
                in_block = False
                continue
            buffer.append(s)

    tokens = []
    for line in buffer:
        tokens.extend([t.strip() for t in line.split(",") if t.strip()])
    idx = 0
    while idx < len(tokens):
        elem_id = int(tokens[idx])
        node_ids = [int(tokens[idx + k]) for k in range(1, 11)]
        elements[elem_id] = node_ids
        idx += 11

    return nodes, elements


def build_facial_surface(case_dir, nodes, elements, target_x=L, tol=1e-6):
    """Find C3D10 tet faces lying on the x=target_x plane - identical
    logic to build_facial_surface.py, cross-checked once already
    against Gmsh's own skin tagging in the baseline case."""
    face_defs = {
        "S1": (0, 1, 2), "S2": (0, 1, 3),
        "S3": (1, 2, 3), "S4": (0, 2, 3),
    }
    found = []
    for elem_id, node_ids in elements.items():
        corners = node_ids[0:4]
        coords = [nodes[n] for n in corners]
        for label, idx_triplet in face_defs.items():
            xs = [coords[i][0] for i in idx_triplet]
            if all(abs(x - target_x) < tol for x in xs):
                found.append((elem_id, label))

    path = os.path.join(case_dir, "tip_facial_surface.inp")
    with open(path, "w") as f:
        f.write("*SURFACE,NAME=TIP_FACE,TYPE=ELEMENT\n")
        for elem_id, label in found:
            f.write(f"{elem_id},{label}\n")
    return path, len(found)


def write_analysis_deck(case_dir):
    """Same physics as the hand-verified baseline deck - only the
    mesh/nset/surface include paths change per case."""
    path = os.path.join(case_dir, "analysis.inp")
    content = f"""** Project 01 - convergence study case
*INCLUDE, INPUT=mesh_clean.inp
*INCLUDE, INPUT=nsets.inp
*INCLUDE, INPUT=tip_facial_surface.inp

*NODE
50000, {L}, {B/2}, {H/2}
*NSET, NSET=REF_TIP
50000

*MATERIAL, NAME=AL7075T6
*ELASTIC
71700., 0.33

*SOLID SECTION, ELSET=plate, MATERIAL=AL7075T6

*COUPLING, REF NODE=50000, SURFACE=TIP_FACE, CONSTRAINT NAME=TIP_COUPLING
*DISTRIBUTING
1,3

*BOUNDARY
root, 1, 3

*STEP
*STATIC

*CLOAD
50000, 3, {LOAD_Z}

*NODE FILE
U
*EL FILE
S
*END STEP
"""
    with open(path, "w") as f:
        f.write(content)
    return path


def run_ccx(case_dir, analysis_inp):
    """Run CalculiX, capture full stdout to run_log.txt, extract
    solve time."""
    base = os.path.splitext(os.path.basename(analysis_inp))[0]
    log_path = os.path.join(case_dir, "run_log.txt")

    t0 = time.time()
    result = subprocess.run(
        ["ccx", base], cwd=case_dir,
        capture_output=True, text=True
    )
    wall_time = time.time() - t0

    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    if "*ERROR" in result.stdout or "*ERROR" in result.stderr:
        raise RuntimeError(
            f"CalculiX reported *ERROR in case {case_dir} - see {log_path}"
        )

    m = re.search(r"Total CalculiX Time:\s*([\d.]+)", result.stdout)
    solve_time = float(m.group(1)) if m else wall_time

    return solve_time


def parse_frd_block(frd_lines, marker, ncols):
    """Generic fixed-width .frd result block parser - validated
    approach from check_frd_displacement.py / check_frd_stress.py."""
    data = {}
    in_block = False
    for line in frd_lines:
        if marker in line and line.rstrip("\n").lstrip().startswith("-4"):
            in_block = True
            continue
        if in_block and line.rstrip("\n").lstrip().startswith("-3"):
            break
        if in_block and line[:3].strip() == "-1":
            nid = int(line[3:13])
            vals = []
            for c in range(ncols):
                start = 13 + c * 12
                end = start + 12
                vals.append(float(line[start:end]))
            data[nid] = vals
    return data


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = s
    return (0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2
            + 6 * (sxy ** 2 + syz ** 2 + szx ** 2))) ** 0.5


def parse_nset(nset_inp, name):
    nodes = set()
    with open(nset_inp) as f:
        lines = f.readlines()
    capture = False
    for line in lines:
        if line.strip().upper() == f"*NSET,NSET={name.upper()}":
            capture = True
            continue
        if capture:
            if line.strip().startswith("*"):
                break
            for tok in line.strip().rstrip(",").split(","):
                tok = tok.strip()
                if tok.isdigit():
                    nodes.add(int(tok))
    return nodes


def extract_results(case_dir, nset_inp, nodes_coords):
    """Extract tip Uz average and root-band max von Mises from this
    case's .frd file."""
    frd_path = os.path.join(case_dir, "analysis.frd")
    with open(frd_path) as f:
        frd_lines = f.readlines()

    disp = parse_frd_block(frd_lines, "DISP", 3)
    stress = parse_frd_block(frd_lines, "STRESS", 6)

    tip_nodes = parse_nset(nset_inp, "tip")
    tip_uz = [disp[n][2] for n in tip_nodes if n in disp]
    tip_uz_avg = sum(tip_uz) / len(tip_uz) if tip_uz else None

    root_band_nodes = [
        nid for nid, (x, y, z) in nodes_coords.items()
        if ROOT_X_MIN <= x <= ROOT_X_MAX
    ]
    root_vm = [von_mises(stress[n]) for n in root_band_nodes if n in stress]
    max_vm_root = max(root_vm) if root_vm else None

    return tip_uz_avg, max_vm_root
