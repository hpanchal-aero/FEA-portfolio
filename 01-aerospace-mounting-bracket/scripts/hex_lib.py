"""
Project 01 - Aerospace Mounting Bracket
Shared library for the STRUCTURED HEX convergence study
(verification plate ONLY - scoped methodology, see chat log).

Reuses the exact validated logic from the standalone hex mesh test:
  - case_hex.geo structure (OpenCASCADE box + Transfinite meshing)
  - generate_hex_nsets.py NSET generation (single Gmsh session)
  - build_hex_facial_surface.py facial surface search
  - analysis_hex.inp physics (material, BC, coupling, load)
with one correction: the reference node ID is now computed as
max_node_id + 1000 for each mesh, rather than a fixed guess - this
is what caused the standalone test's first solve failure (node 50000
collided with a real mesh node at 88386-node density).
"""

import os
import re
import subprocess
import time

import gmsh

L, B, H = 60.0, 40.0, 4.0
LOAD_Z = -235.44

ROOT_X_MIN = 3.0
ROOT_X_MAX = 9.0

FACE_DEFS = {
    "S1": (0,1,2,3), "S2": (4,5,6,7),
    "S3": (0,1,5,4), "S4": (1,2,6,5),
    "S5": (2,3,7,6), "S6": (3,0,4,7),
}

HEX_EDGES = [
    (0,1),(1,2),(2,3),(3,0),
    (4,5),(5,6),(6,7),(7,4),
    (0,4),(1,5),(2,6),(3,7),
]


def write_geo(case_dir, nx, ny, nz):
    """nx,ny,nz are POINT counts (elements = n-1) along L,B,H."""
    geo_path = os.path.join(case_dir, "case.geo")
    eps = 1e-3
    content = f"""// Project 01 - Aerospace Mounting Bracket
// Structured hex convergence case: nx={nx} ny={ny} nz={nz}
// Scoped to verification plate ONLY - see chat log.

SetFactory("OpenCASCADE");
L = {L}; B = {B}; H = {H};
Box(1) = {{0, 0, 0, L, B, H}};

eps = {eps};

xline1() = Curve In BoundingBox{{-eps,-eps,-eps, L+eps,eps,eps}};
xline2() = Curve In BoundingBox{{-eps,-eps,H-eps, L+eps,eps,H+eps}};
xline3() = Curve In BoundingBox{{-eps,B-eps,-eps, L+eps,B+eps,eps}};
xline4() = Curve In BoundingBox{{-eps,B-eps,H-eps, L+eps,B+eps,H+eps}};

yline1() = Curve In BoundingBox{{-eps,-eps,-eps, eps,B+eps,eps}};
yline2() = Curve In BoundingBox{{-eps,-eps,H-eps, eps,B+eps,H+eps}};
yline3() = Curve In BoundingBox{{L-eps,-eps,-eps, L+eps,B+eps,eps}};
yline4() = Curve In BoundingBox{{L-eps,-eps,H-eps, L+eps,B+eps,H+eps}};

zline1() = Curve In BoundingBox{{-eps,-eps,-eps, eps,eps,H+eps}};
zline2() = Curve In BoundingBox{{-eps,B-eps,-eps, eps,B+eps,H+eps}};
zline3() = Curve In BoundingBox{{L-eps,-eps,-eps, L+eps,eps,H+eps}};
zline4() = Curve In BoundingBox{{L-eps,B-eps,-eps, L+eps,B+eps,H+eps}};

Transfinite Curve {{xline1(), xline2(), xline3(), xline4()}} = {nx};
Transfinite Curve {{yline1(), yline2(), yline3(), yline4()}} = {ny};
Transfinite Curve {{zline1(), zline2(), zline3(), zline4()}} = {nz};

all_surfaces() = Surface In BoundingBox{{-eps,-eps,-eps, L+eps,B+eps,H+eps}};
Transfinite Surface {{all_surfaces()}};
Recombine Surface {{all_surfaces()}};
Transfinite Volume {{1}};

root_surf()      = Surface In BoundingBox{{-eps,-eps,-eps, eps,B+eps,H+eps}};
tip_surf()       = Surface In BoundingBox{{L-eps,-eps,-eps, L+eps,B+eps,H+eps}};
bottom_surf()    = Surface In BoundingBox{{-eps,-eps,-eps, L+eps,B+eps,eps}};
top_surf()       = Surface In BoundingBox{{-eps,-eps,H-eps, L+eps,B+eps,H+eps}};
side_y0_surf()   = Surface In BoundingBox{{-eps,-eps,-eps, L+eps,eps,H+eps}};
side_y40_surf()  = Surface In BoundingBox{{-eps,B-eps,-eps, L+eps,B+eps,H+eps}};

Physical Surface("root", 101)     = root_surf();
Physical Surface("tip", 102)      = tip_surf();
Physical Surface("bottom", 103)   = bottom_surf();
Physical Surface("top", 104)      = top_surf();
Physical Surface("side_y0", 105)  = side_y0_surf();
Physical Surface("side_y40", 106) = side_y40_surf();
Physical Volume("plate", 201)     = {{1}};

Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh.SecondOrderLinear = 1;

Mesh 3;
"""
    with open(geo_path, "w") as f:
        f.write(content)
    return geo_path


def mesh_and_export(case_dir, geo_path):
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
    clean_inp = os.path.join(case_dir, "mesh_clean.inp")
    with open(raw_inp) as f:
        lines = f.readlines()

    out_lines = []
    skip = False
    for line in lines:
        upper = line.strip().upper().replace(" ", "")
        if upper.startswith("*ELEMENT,TYPE=CPS8"):
            skip = True
            continue
        if upper.startswith("*ELEMENT,TYPE=C3D20"):
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
        if s.upper().startswith("*ELEMENT") and "TYPE=C3D20" in s.upper().replace(" ", ""):
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
        node_ids = [int(tokens[idx + k]) for k in range(1, 21)]
        elements[elem_id] = node_ids
        idx += 21

    return nodes, elements


def build_facial_surface(case_dir, nodes, elements, target_x=L, tol=1e-6):
    found = []
    for elem_id, node_ids in elements.items():
        corners = node_ids[0:8]
        coords = [nodes[n] for n in corners]
        for label, idx4 in FACE_DEFS.items():
            xs = [coords[i][0] for i in idx4]
            if all(abs(x - target_x) < tol for x in xs):
                found.append((elem_id, label))

    path = os.path.join(case_dir, "tip_facial_surface.inp")
    with open(path, "w") as f:
        f.write("*SURFACE,NAME=TIP_FACE,TYPE=ELEMENT\n")
        for elem_id, label in found:
            f.write(f"{elem_id},{label}\n")
    return path, len(found)


def write_analysis_deck(case_dir, ref_node_id):
    path = os.path.join(case_dir, "analysis.inp")
    content = f"""** Project 01 - structured hex convergence case
** Scoped to verification plate ONLY - see chat log.
*INCLUDE, INPUT=mesh_clean.inp
*INCLUDE, INPUT=nsets.inp
*INCLUDE, INPUT=tip_facial_surface.inp

*NODE
{ref_node_id}, {L}, {B/2}, {H/2}
*NSET, NSET=REF_TIP
{ref_node_id}

*MATERIAL, NAME=AL7075T6
*ELASTIC
71700., 0.33

*SOLID SECTION, ELSET=plate, MATERIAL=AL7075T6

*BOUNDARY
root, 1, 3

*COUPLING, REF NODE={ref_node_id}, SURFACE=TIP_FACE, CONSTRAINT NAME=TIP_COUPLING
*DISTRIBUTING
1,3

*STEP
*STATIC, SOLVER=ITERATIVE SCALING

*CLOAD
{ref_node_id}, 3, {LOAD_Z}

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
    base = os.path.splitext(os.path.basename(analysis_inp))[0]
    log_path = os.path.join(case_dir, "run_log.txt")

    t0 = time.time()
    result = subprocess.run(["ccx", base], cwd=case_dir,
                             capture_output=True, text=True)
    wall_time = time.time() - t0

    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    if "*ERROR" in result.stdout or "*ERROR" in result.stderr:
        raise RuntimeError(
            f"CalculiX reported *ERROR in case {case_dir} - see {log_path}"
        )

    m = re.search(r"Total CalculiX Time:\s*([\d.]+)", result.stdout)
    return float(m.group(1)) if m else wall_time


def parse_frd_block(frd_lines, marker, ncols):
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
    return (0.5 * ((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2
            + 6*(sxy**2 + syz**2 + szx**2))) ** 0.5


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
    frd_path = os.path.join(case_dir, "analysis.frd")
    with open(frd_path) as f:
        frd_lines = f.readlines()

    disp = parse_frd_block(frd_lines, "DISP", 3)
    stress = parse_frd_block(frd_lines, "STRESS", 6)

    tip_nodes = parse_nset(nset_inp, "tip")
    tip_uz = [disp[n][2] for n in tip_nodes if n in disp]
    tip_uz_avg = sum(tip_uz) / len(tip_uz) if tip_uz else None

    root_band_nodes = [nid for nid, (x, y, z) in nodes_coords.items()
                        if ROOT_X_MIN <= x <= ROOT_X_MAX]
    root_vm = [von_mises(stress[n]) for n in root_band_nodes if n in stress]
    max_vm_root = max(root_vm) if root_vm else None

    return tip_uz_avg, max_vm_root
