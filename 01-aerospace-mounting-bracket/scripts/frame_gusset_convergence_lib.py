"""
Stage 4b (right-angle frame with gusset) - targeted convergence study
at the gusset's tapered toe (x=38, z=2, spanning y) - the line where
the hypotenuse face meets Flange A's concave surface, identified as
the location of the ~-274 MPa peak in the baseline solve.

Geometry, material, BCs, load are UNCHANGED from the baseline
(build_frame_gusset_case.py) - only mesh refinement is varied, and
only in a NEW field targeted specifically at the toe edge, combined
with (not replacing) the existing hypotenuse-face field.

The baseline solve (already run) used only the hypotenuse-face
Distance field with SizeMin=0.5mm - this implicitly gave the toe
region ~0.5mm resolution by proximity, with no dedicated targeting.
This library adds a second Distance+Threshold field on TOE_SIZE,
anchored on the toe edge specifically, combined via a Min field.
"""

import gmsh
import sys
import os

TOL = 1e-6

L = 60.0
WIDTH = 40.0
THICK = 4.0
HALF_T = THICK / 2.0
GUSSET_LEG = 20.0

FAR_SIZE = 3.0
HYP_SIZE = 0.5          # unchanged from baseline - hypotenuse face field
TRANSITION_WIDTH = 4.0


def write_and_mesh(case_dir, toe_size):
    case_dir = os.path.abspath(case_dir)
    os.makedirs(case_dir, exist_ok=True)
    out_mesh = os.path.join(case_dir, "mesh_frame_gusset_raw.inp")
    nset_out = os.path.join(case_dir, "nsets_frame_gusset.inp")

    gmsh.initialize()
    gmsh.model.add("frame_gusset_case")
    occ = gmsh.model.occ

    flange_a = occ.addBox(0.0, 0.0, -HALF_T, L, WIDTH, THICK)
    flange_b = occ.addBox(L - HALF_T, 0.0, 0.0, THICK, WIDTH, L)
    occ.synchronize()

    corner_x = L - HALF_T
    corner_z = HALF_T

    p1 = occ.addPoint(corner_x, 0.0, corner_z)
    p2 = occ.addPoint(corner_x - GUSSET_LEG, 0.0, corner_z)
    p3 = occ.addPoint(corner_x, 0.0, corner_z + GUSSET_LEG)
    l1 = occ.addLine(p1, p2)
    l2 = occ.addLine(p2, p3)
    l3 = occ.addLine(p3, p1)
    wire = occ.addCurveLoop([l1, l2, l3])
    gusset_face = occ.addPlaneSurface([wire])
    occ.synchronize()

    gusset_extrude = occ.extrude([(2, gusset_face)], 0.0, WIDTH, 0.0)
    occ.synchronize()
    gusset_vol = None
    for dim, tag in gusset_extrude:
        if dim == 3:
            gusset_vol = tag
            break

    out, out_map = occ.fuse([(3, flange_a)], [(3, flange_b), (3, gusset_vol)])
    occ.synchronize()
    fused_vol = out[0][1]

    # ---- Identify hypotenuse face (as before) ----
    all_surfaces_search = gmsh.model.getEntities(2)
    hyp_face = None
    for dim, tag in all_surfaces_search:
        stype = gmsh.model.getType(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        spans_width = (ymax - ymin) > (WIDTH - 1.0)
        x_diagonal = abs(xmin - (corner_x - GUSSET_LEG)) < TOL and abs(xmax - corner_x) < TOL
        z_diagonal = abs(zmin - corner_z) < TOL and abs(zmax - (corner_z + GUSSET_LEG)) < TOL
        if stype == "Plane" and spans_width and x_diagonal and z_diagonal:
            hyp_face = tag
            break

    if hyp_face is None:
        print(f"FATAL [{case_dir}]: could not identify hypotenuse face.")
        gmsh.finalize()
        sys.exit(1)

    # ---- Identify the toe edge: x=(corner_x - GUSSET_LEG)=38, z=corner_z=2,
    # spanning full y-width - a straight, width-spanning, constant-x,
    # constant-z edge, same signature style as 4a's original corner edge ----
    all_curves = gmsh.model.getEntities(1)
    toe_edge = None
    for dim, tag in all_curves:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        x_at_toe = abs(xmin - (corner_x - GUSSET_LEG)) < TOL and abs(xmax - (corner_x - GUSSET_LEG)) < TOL
        z_at_toe = abs(zmin - corner_z) < TOL and abs(zmax - corner_z) < TOL
        spans_width = (ymax - ymin) > (WIDTH - 1.0)
        if x_at_toe and z_at_toe and spans_width:
            toe_edge = tag
            break

    if toe_edge is None:
        print(f"FATAL [{case_dir}]: could not identify the toe edge.")
        gmsh.finalize()
        sys.exit(1)

    # ---- Root / tip faces ----
    surfaces = gmsh.model.getEntities(2)
    root_face_tag = None
    tip_face_tag = None
    for dim, tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        if abs(xmin) < TOL and abs(xmax) < TOL:
            root_face_tag = tag
        if abs(zmin - L) < TOL and abs(zmax - L) < TOL:
            tip_face_tag = tag

    if root_face_tag is None or tip_face_tag is None:
        print(f"FATAL [{case_dir}]: could not identify root/tip face.")
        gmsh.finalize()
        sys.exit(1)

    gmsh.model.addPhysicalGroup(3, [fused_vol], name="frame_gusset")
    gmsh.model.addPhysicalGroup(2, [root_face_tag], name="root")
    gmsh.model.addPhysicalGroup(2, [tip_face_tag], name="tip")

    # ---- Field 1: hypotenuse face (unchanged from baseline) ----
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [hyp_face])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", HYP_SIZE)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", FAR_SIZE)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", TRANSITION_WIDTH)
    gmsh.model.mesh.field.setNumber(2, "Sigmoid", 1)

    # ---- Field 2: toe edge (NEW, swept parameter) ----
    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "CurvesList", [toe_edge])
    gmsh.model.mesh.field.setNumber(3, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 3)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", toe_size)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", FAR_SIZE)
    gmsh.model.mesh.field.setNumber(4, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(4, "DistMax", TRANSITION_WIDTH)
    gmsh.model.mesh.field.setNumber(4, "Sigmoid", 1)

    # ---- Combine: take the smaller (finer) size at every point ----
    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
    gmsh.model.mesh.field.setAsBackgroundMesh(5)

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

    named_surfaces = ["root", "tip"]
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
    print(f"[{os.path.basename(case_dir)}] toe_size={toe_size}mm -> "
          f"{n_nodes} nodes, {n_elements} elements")
    return n_nodes, n_elements
