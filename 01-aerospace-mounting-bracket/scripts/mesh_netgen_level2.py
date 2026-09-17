"""
Stage 5 -- Netgen mesh, Level 2 of 3-level toe-fillet convergence study.
Identical methodology to Level 1 (mesh_netgen.py): straight-sided C3D10
via deterministic midpoint promotion from a verified-positive-volume
Netgen linear mesh. Only the local feature sizes change.

Level 1: hole=0.50, R3=0.50, toe=0.25mm
Level 2: hole=0.35, R3=0.35, toe=0.15mm  <-- this script
Level 3: hole=0.25, R3=0.25, toe=0.10mm
"""

from netgen.occ import OCCGeometry
from netgen.meshing import Point3d
import numpy as np
import os

BREP_FILE = "mesh/full_bracket_study/full_bracket_geometry_final.brep"
INP_OUT = "simulation/stage5_netgen_level2_solve.inp"

E_MOD = 71700.0
NU = 0.33
RHO = 2810e-9
FORCE_X = -235.44
REF_NODE_COORD = (60.0, 20.0, 63.0)

MAXH_GLOBAL = 6.0
FEATURE_POINTS = [
    (20.0, 20.0, 0.0, 0.35),   # hole center
    (56.5, 20.0, 3.5, 0.35),   # R3 fillet region
    (38.0, 20.0, 2.0, 0.15),   # toe fillet
]

def main():
    if not os.path.exists(BREP_FILE):
        print(f"ERROR: {BREP_FILE} not found")
        return

    print("Loading BRep via netgen.occ...")
    geo = OCCGeometry(BREP_FILE)

    for x, y, z, h in FEATURE_POINTS:
        geo.RestrictH(Point3d(x, y, z), h)

    print(f"Generating mesh (global maxh={MAXH_GLOBAL}, with local refinement)...")
    mesh = geo.GenerateMesh(maxh=MAXH_GLOBAL)
    print("Mesh generated.")

    points = mesh.Points()
    n_pts = len(points)
    print(f"Points: {n_pts}")
    coords = np.zeros((n_pts + 1, 3))
    for i, p in enumerate(points, start=1):
        pnt = p.p
        coords[i] = [pnt[0], pnt[1], pnt[2]]

    elements = mesh.Elements3D()
    n_elem = len(elements)
    print(f"Elements: {n_elem}")

    tets = []
    for el in elements:
        pts = [v.nr for v in el.vertices]
        if len(pts) != 4:
            print(f"WARNING: non-tet element with {len(pts)} nodes, skipping")
            continue
        tets.append(pts)
    tets = np.array(tets)
    tets = tets[:, [0, 1, 3, 2]]  # Netgen->CalculiX orientation fix, confirmed in Level 1
    print(f"Linear tets extracted: {len(tets)}")

    def signed_volume(p0, p1, p2, p3):
        return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

    vols = np.array([
        signed_volume(coords[t[0]], coords[t[1]], coords[t[2]], coords[t[3]])
        for t in tets
    ])
    n_bad = (vols <= 0).sum()
    print(f"Linear volume check: min={vols.min():.6f}, mean={vols.mean():.6f}, "
          f"negative/zero count={n_bad}")
    if n_bad > 0:
        print("ERROR: non-positive volume elements found. Stopping.")
        return
    print("All linear tets have positive volume -- safe to promote to quadratic.")

    print("\nPromoting to C3D10 (straight-sided, deterministic midside nodes)...")
    edge_pairs = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]

    midside_cache = {}
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

    total_nodes = n_pts + len(midside_cache)
    print(f"Unique midside nodes created: {len(midside_cache)}")
    print(f"Total nodes: {total_nodes}")
    print(f"Total elements: {len(quad_elems)}")

    all_node_tags = list(range(1, n_pts + 1)) + list(extra_coords.keys())
    def get_coord(tag):
        return coords[tag] if tag <= n_pts else extra_coords[tag]

    root_nodes = [t for t in all_node_tags if abs(get_coord(t)[0] - 0.0) < 1e-3]
    tip_nodes = [t for t in all_node_tags if abs(get_coord(t)[2] - 60.0) < 1e-3]
    print(f"\nRoot nodes (x=0): {len(root_nodes)}")
    print(f"Tip nodes (z=60): {len(tip_nodes)}")

    if len(root_nodes) == 0 or len(tip_nodes) == 0:
        print("ERROR: root or tip node set empty.")
        return

    est_equations = 3 * total_nodes
    ceiling = 1_180_000
    print(f"\nEquation estimate: {est_equations} ({100*est_equations/ceiling:.1f}% of ceiling)")

    os.makedirs(os.path.dirname(INP_OUT), exist_ok=True)
    ref_tag = total_nodes + 1
    with open(INP_OUT, "w") as f:
        f.write("*HEADING\n")
        f.write("Stage 5 -- Netgen mesh Level 2, deterministic straight C3D10\n")
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

if __name__ == "__main__":
    main()
