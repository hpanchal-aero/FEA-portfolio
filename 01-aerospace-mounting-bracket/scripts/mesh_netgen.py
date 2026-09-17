"""
Stage 5 -- alternative meshing path using Netgen instead of Gmsh, as a
diagnostic test of whether the persistent order-2 Jacobian failure is
specific to Gmsh's tetrahedralization.

Strategy:
  1. Load the verified BRep directly via netgen.occ.
  2. Generate a LINEAR (order-1) tetrahedral mesh only -- Netgen's own
     native algorithm, completely independent of Gmsh.
  3. Verify every element has positive linear volume (a necessary
     condition we can check directly).
  4. Promote to 10-node quadratic elements OURSELVES, by placing each
     midside node at the exact geometric midpoint of its parent edge.
     This guarantees (by the geometric identity proven earlier this
     session) that the resulting C3D10 elements are mathematically
     equivalent to their linear parents everywhere -- eliminating any
     ambiguity from Netgen's or Gmsh's own order-2 elevation/curving
     logic, which has been a repeated source of confusion this session.
  5. Write directly to a CalculiX-ready .inp in the documented
     Abaqus/CalculiX C3D10 node order: corners 1-4, then edges
     (1,2)(2,3)(1,3)(1,4)(2,4)(3,4).
"""

from netgen.occ import OCCGeometry
import numpy as np
import os

BREP_FILE = "mesh/full_bracket_study/full_bracket_geometry_final.brep"
INP_OUT = "simulation/stage5_netgen_solve.inp"

E_MOD = 71700.0
NU = 0.33
RHO = 2810e-9
FORCE_X = -235.44
REF_NODE_COORD = (60.0, 20.0, 63.0)

# --- Sizing: global coarse + local refinement at the 3 features, via
# point-based RestrictH (Netgen's local-density control mechanism) ---
MAXH_GLOBAL = 6.0
FEATURE_POINTS = [
    # (x, y, z, target_h)
    (20.0, 20.0, 0.0, 0.5),    # hole center
    (56.5, 20.0, 3.5, 0.5),    # R3 fillet region
    (38.0, 20.0, 2.0, 0.25),   # toe fillet
]

def main():
    if not os.path.exists(BREP_FILE):
        print(f"ERROR: {BREP_FILE} not found")
        return

    print("Loading BRep via netgen.occ...")
    geo = OCCGeometry(BREP_FILE)

    from netgen.meshing import Point3d
    for x, y, z, h in FEATURE_POINTS:
        geo.RestrictH(Point3d(x, y, z), h)

    print(f"Generating mesh (global maxh={MAXH_GLOBAL}, with local refinement)...")
    mesh = geo.GenerateMesh(maxh=MAXH_GLOBAL)

    n_points = mesh.Coordinates().shape[0] if hasattr(mesh, "Coordinates") else len(mesh.Points())
    print(f"Mesh generated.")

    # --- Extract points ---
    points = mesh.Points()
    n_pts = len(points)
    print(f"Points: {n_pts}")
    coords = np.zeros((n_pts + 1, 3))  # 1-indexed to match netgen point numbering
    for i, p in enumerate(points, start=1):
        pnt = p.p
        coords[i] = [pnt[0], pnt[1], pnt[2]]

    # --- Extract volume elements (linear tets, 4 corner nodes each) ---
    elements = mesh.Elements3D()
    n_elem = len(elements)
    print(f"Elements: {n_elem}")

    tets = []
    for el in elements:
        pts = [v.nr for v in el.vertices]  # 1-indexed netgen point numbers
        if len(pts) != 4:
            print(f"WARNING: non-tet element with {len(pts)} nodes found, skipping")
            continue
        tets.append(pts)
    tets = np.array(tets)
    print(f"Linear tets extracted: {len(tets)}")

    # --- Verify positive volume for every linear tet ---
    def signed_volume(p0, p1, p2, p3):
        return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

    # Netgen's tet vertex order is the mirror-image convention of the
    # standard positive-orientation formula above (confirmed: ALL 11,250
    # elements reported negative volume -- a uniform-sign result means a
    # convention mismatch, not real degeneracy). Swap nodes 2 and 3 in
    # the volume check AND in every element's actual connectivity below
    # to correct orientation to CalculiX's expected positive convention.
    tets = tets[:, [0, 1, 3, 2]]  # swap last two corner nodes, all elements

    vols = np.array([
        signed_volume(coords[t[0]], coords[t[1]], coords[t[2]], coords[t[3]])
        for t in tets
    ])
    n_bad = (vols <= 0).sum()
    print(f"Linear volume check: min={vols.min():.6f}, mean={vols.mean():.6f}, "
          f"negative/zero count={n_bad}")
    if n_bad > 0:
        print("ERROR: some elements have non-positive linear volume -- "
              "Netgen itself produced degenerate elements. Stopping.")
        return
    print("All linear tets have positive volume -- safe to promote to quadratic.")

    # --- Promote to quadratic: build midside nodes at exact edge midpoints ---
    print("\nPromoting to C3D10 (straight-sided, deterministic midside nodes)...")
    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]  # CalculiX C3D10 order

    midside_cache = {}  # (min_node, max_node) -> new node tag
    next_node_tag = n_pts + 1
    extra_coords = {}

    quad_elems = []
    for t in tets:
        full_conn = list(t)
        for (a, b) in edge_pairs:
            na, nb = t[a], t[b]
            key = (min(na, nb), max(na, nb))
            if key not in midside_cache:
                mid = (coords[na] + coords[nb]) / 2.0
                extra_coords[next_node_tag] = mid
                midside_cache[key] = next_node_tag
                next_node_tag += 1
            full_conn.append(midside_cache[key])
        quad_elems.append(full_conn)

    print(f"Unique midside nodes created: {len(midside_cache)}")
    total_nodes = n_pts + len(midside_cache)
    print(f"Total nodes (corners + midsides): {total_nodes}")
    print(f"Total elements: {len(quad_elems)}")

    # --- Identify root (x=0) and tip (z=60) nodes for BCs ---
    all_node_tags = list(range(1, n_pts + 1)) + list(extra_coords.keys())
    def get_coord(tag):
        return coords[tag] if tag <= n_pts else extra_coords[tag]

    root_nodes = [t for t in all_node_tags if abs(get_coord(t)[0] - 0.0) < 1e-3]
    tip_nodes = [t for t in all_node_tags if abs(get_coord(t)[2] - 60.0) < 1e-3]
    print(f"\nRoot nodes (x=0): {len(root_nodes)}")
    print(f"Tip nodes (z=60): {len(tip_nodes)}")

    if len(root_nodes) == 0 or len(tip_nodes) == 0:
        print("ERROR: root or tip node set empty -- check geometry orientation.")
        return

    # --- Write .inp ---
    os.makedirs(os.path.dirname(INP_OUT), exist_ok=True)
    ref_tag = total_nodes + 1
    with open(INP_OUT, "w") as f:
        f.write("*HEADING\n")
        f.write("Stage 5 -- Netgen-generated mesh, deterministic straight C3D10 promotion\n")
        f.write("Al 7075-T6, F=235.44N in -x at Flange B tip, root encastre\n")

        f.write("*NODE\n")
        for tag in range(1, n_pts + 1):
            c = coords[tag]
            f.write(f"{tag}, {c[0]!r}, {c[1]!r}, {c[2]!r}\n")
        for tag, c in extra_coords.items():
            f.write(f"{tag}, {c[0]!r}, {c[1]!r}, {c[2]!r}\n")
        f.write(f"{ref_tag}, {REF_NODE_COORD[0]}, {REF_NODE_COORD[1]}, {REF_NODE_COORD[2]}\n")

        f.write("*ELEMENT, TYPE=C3D10, ELSET=EBRACKET\n")
        for i, conn in enumerate(quad_elems, start=1):
            f.write(f"{i}, " + ", ".join(str(n) for n in conn) + "\n")

        f.write("*NSET, NSET=NROOT\n")
        for i in range(0, len(root_nodes), 10):
            f.write(", ".join(str(n) for n in root_nodes[i:i+10]) + "\n")
        f.write("*NSET, NSET=NTIP\n")
        for i in range(0, len(tip_nodes), 10):
            f.write(", ".join(str(n) for n in tip_nodes[i:i+10]) + "\n")

        f.write("*SURFACE, NAME=STIP, TYPE=NODE\nNTIP\n")

        f.write("*MATERIAL, NAME=AL7075T6\n*ELASTIC\n")
        f.write(f"{E_MOD}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=EBRACKET, MATERIAL=AL7075T6\n")

        f.write("*STEP\n*STATIC\n")
        f.write("*BOUNDARY\nNROOT, 1, 3, 0.0\n")
        f.write(f"*COUPLING, CONSTRAINT NAME=CTIP, REF NODE={ref_tag}, SURFACE=STIP\n")
        f.write("*KINEMATIC\n1, 3\n")
        f.write(f"*CLOAD\n{ref_tag}, 1, {FORCE_X}\n")
        f.write("*NODE FILE\nU, RF\n*EL FILE\nS, E\n*END STEP\n")

    print(f"\nWrote: {INP_OUT}")
    print(f"Reference node: {ref_tag}")
    print(f"Total DOF estimate: {3 * total_nodes}")

if __name__ == "__main__":
    main()
