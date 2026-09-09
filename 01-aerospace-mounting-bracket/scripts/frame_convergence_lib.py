"""
Stage 4a (right-angle frame) - parametrized geometry+mesh generation
for the convergence sweep. Refactors the validated logic from
build_frame_case.py, varying only CORNER_SIZE.

Locked/fixed for every level (NOT swept):
  - Geometry: 60mm flanges, 40x4mm section, sharp reentrant corner
  - FAR_SIZE=3.0, TRANSITION_WIDTH=4.0, Sigmoid interpolation
Varied per level:
  - CORNER_SIZE only

Convergence target: mid-span stations only (Flange A x=30, Flange B
s=30 i.e. z=30) - NOT the corner (s=60), which is a known,
non-convergent stress singularity by design (sharp reentrant corner,
no fillet in this sub-stage) and is explicitly excluded from the
convergence criterion.
"""

import gmsh
import sys
import os

TOL = 1e-6

L = 60.0
WIDTH = 40.0
THICK = 4.0
HALF_T = THICK / 2.0

FAR_SIZE = 3.0
TRANSITION_WIDTH = 4.0


def write_and_mesh(case_dir, corner_size):
    case_dir = os.path.abspath(case_dir)
    os.makedirs(case_dir, exist_ok=True)
    out_mesh = os.path.join(case_dir, "mesh_frame_raw.inp")
    nset_out = os.path.join(case_dir, "nsets_frame.inp")

    gmsh.initialize()
    gmsh.model.add("frame_case")
    occ = gmsh.model.occ

    flange_a = occ.addBox(0.0, 0.0, -HALF_T, L, WIDTH, THICK)
    flange_b = occ.addBox(L - HALF_T, 0.0, 0.0, THICK, WIDTH, L)
    occ.synchronize()

    out, out_map = occ.fuse([(3, flange_a)], [(3, flange_b)])
    occ.synchronize()
    fused_vol = out[0][1]

    all_curves = gmsh.model.getEntities(1)
    corner_edge = None
    for dim, tag in all_curves:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        x_at_corner = abs(xmin - (L - HALF_T)) < TOL and abs(xmax - (L - HALF_T)) < TOL
        z_at_corner = abs(zmin - HALF_T) < TOL and abs(zmax - HALF_T) < TOL
        spans_width = (ymax - ymin) > (WIDTH - 1.0)
        if x_at_corner and z_at_corner and spans_width:
            corner_edge = tag
            break

    if corner_edge is None:
        print(f"FATAL [{case_dir}]: could not identify the reentrant corner edge.")
        gmsh.finalize()
        sys.exit(1)

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
        print(f"FATAL [{case_dir}]: could not identify root and/or tip face.")
        gmsh.finalize()
        sys.exit(1)

    gmsh.model.addPhysicalGroup(3, [fused_vol], name="frame")
    gmsh.model.addPhysicalGroup(2, [root_face_tag], name="root")
    gmsh.model.addPhysicalGroup(2, [tip_face_tag], name="tip")

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [corner_edge])
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", corner_size)
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
    print(f"[{os.path.basename(case_dir)}] corner_size={corner_size}mm -> "
          f"{n_nodes} nodes, {n_elements} elements")
    return n_nodes, n_elements
