"""
Project 01 - Aerospace Mounting Bracket
Stage 3 (fillet study): parametrized geometry+mesh generation for the
convergence sweep. Refactors the validated logic from
build_fillet_case.py (fillet radius fix at r=0.95mm, robust
type+z-sign surface identification, sigmoid Distance+Threshold field)
into a reusable function, varying only FILLET_SIZE.

Locked/fixed for every level (NOT swept):
  - Geometry: 60x40mm, 6mm->4mm thickness step at x=30, r=0.95mm
  - FAR_SIZE=3.0, TRANSITION_WIDTH=1.5, Sigmoid interpolation
Varied per level:
  - FILLET_SIZE only
"""

import gmsh
import sys
import os

TOL = 1e-6

LENGTH = 60.0
WIDTH  = 40.0
D_ROOT = 6.0
D_TIP  = 4.0
X_TRANS = 30.0
FILLET_R = 0.95

HALF_ROOT = D_ROOT / 2.0
HALF_TIP  = D_TIP / 2.0

FAR_SIZE = 3.0
TRANSITION_WIDTH = 1.5


def write_and_mesh(case_dir, fillet_size):
    """
    Builds the Stage 3 fillet geometry, applies the sigmoid-blended
    Distance+Threshold refinement at the given fillet_size, meshes,
    and writes mesh_fillet_raw.inp + nsets_fillet.inp into case_dir.

    Returns (n_nodes, n_elements).
    """
    case_dir = os.path.abspath(case_dir)
    os.makedirs(case_dir, exist_ok=True)
    out_mesh = os.path.join(case_dir, "mesh_fillet_raw.inp")
    nset_out = os.path.join(case_dir, "nsets_fillet.inp")

    gmsh.initialize()
    gmsh.model.add("fillet_case")
    occ = gmsh.model.occ

    root_box = occ.addBox(0.0, 0.0, -HALF_ROOT, X_TRANS, WIDTH, D_ROOT)
    tip_box  = occ.addBox(X_TRANS, 0.0, -HALF_TIP, LENGTH - X_TRANS, WIDTH, D_TIP)
    occ.synchronize()

    out, out_map = occ.fuse([(3, root_box)], [(3, tip_box)])
    occ.synchronize()
    fused_vol = out[0][1]

    all_curves = gmsh.model.getEntities(1)
    top_edges, bottom_edges = [], []
    for dim, tag in all_curves:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        x_at_trans = abs(xmin - X_TRANS) < TOL and abs(xmax - X_TRANS) < TOL
        spans_width = (ymax - ymin) > (WIDTH - 1.0)
        if not (x_at_trans and spans_width):
            continue
        z_mid = 0.5 * (zmin + zmax)
        if (zmax - zmin) >= TOL:
            continue
        if abs(z_mid - HALF_TIP) < TOL:
            top_edges.append(tag)
        elif abs(z_mid + HALF_TIP) < TOL:
            bottom_edges.append(tag)

    if not top_edges or not bottom_edges:
        print(f"FATAL [{case_dir}]: could not identify notch edges.")
        gmsh.finalize()
        sys.exit(1)

    fillet_edges = top_edges + bottom_edges
    radii = [FILLET_R] * len(fillet_edges)

    filleted = occ.fillet([fused_vol], fillet_edges, radii, True)
    occ.synchronize()
    final_vol = filleted[0][1]

    surfaces = gmsh.model.getEntities(2)
    top_fillet_surf = None
    bottom_fillet_surf = None
    root_face_tag = None
    for dim, tag in surfaces:
        stype = gmsh.model.getType(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        if stype == "Cylinder":
            z_mid = 0.5 * (zmin + zmax)
            if z_mid > 0:
                top_fillet_surf = tag
            elif z_mid < 0:
                bottom_fillet_surf = tag
        if abs(xmin) < TOL and abs(xmax) < TOL:
            root_face_tag = tag

    if top_fillet_surf is None or bottom_fillet_surf is None or root_face_tag is None:
        print(f"FATAL [{case_dir}]: could not identify required surfaces.")
        gmsh.finalize()
        sys.exit(1)

    gmsh.model.addPhysicalGroup(3, [final_vol], name="plate")
    gmsh.model.addPhysicalGroup(2, [root_face_tag], name="root")
    gmsh.model.addPhysicalGroup(2, [top_fillet_surf], name="fillet_top")
    gmsh.model.addPhysicalGroup(2, [bottom_fillet_surf], name="fillet_bottom")

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [top_fillet_surf, bottom_fillet_surf])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", fillet_size)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", FAR_SIZE)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", TRANSITION_WIDTH)
    gmsh.model.mesh.field.setNumber(2, "Sigmoid", 1)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    gmsh.model.mesh.generate(3)
    gmsh.write(out_mesh)

    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(t) for t in elem_tags)

    named_surfaces = ["root", "fillet_top", "fillet_bottom"]
    groups = gmsh.model.getPhysicalGroups(dim=2)
    name_to_tag = {gmsh.model.getPhysicalName(d, t): t for d, t in groups}
    with open(nset_out, "w") as f:
        for name in named_surfaces:
            tag = name_to_tag[name]
            node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, tag)
            node_tags = sorted(int(nd) for nd in node_tags)
            f.write(f"*NSET,NSET={name}\n")
            for i in range(0, len(node_tags), 8):
                chunk = node_tags[i:i + 8]
                f.write(",".join(str(nd) for nd in chunk) + ",\n")

    gmsh.finalize()
    print(f"[{os.path.basename(case_dir)}] fillet_size={fillet_size}mm -> "
          f"{n_nodes} nodes, {n_elements} elements")
    return n_nodes, n_elements
