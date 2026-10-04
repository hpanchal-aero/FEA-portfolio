#!/usr/bin/env python3
"""
Diagnostic: classify boundary faces of the root rib in the Stage 2 solid mesh.

Question: is the rib's outer rim a distinct meshed surface (boundary triangles
spanning z=-0.75 to z=+0.75), or only an edge where the two flat faces meet?

Method: true boundary faces = corner-triples appearing in exactly ONE C3D10.
Corner nodes only (first 4 nodes of each C3D10); midside nodes not classified.

Run with: python3 scripts/diagnose_rib_rim.py simulation/stage2_ribs_h2p0.inp
(inside `conda activate ccx223`)
"""
import sys
import numpy as np

TOL = 1e-3  # mm


def parse_inp(path):
    nodes = {}
    elems = []
    state = None
    buf = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                # flush any partial element
                buf = []
                key = line.upper()
                if key.startswith("*NODE") and not key.startswith("*NODE PRINT") \
                        and not key.startswith("*NODE FILE"):
                    state = "node"
                elif key.startswith("*ELEMENT") and "C3D10" in key:
                    state = "elem"
                else:
                    state = None
                continue
            if state == "node":
                p = [t.strip() for t in line.split(",") if t.strip()]
                nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            elif state == "elem":
                buf += [int(t) for t in line.split(",") if t.strip()]
                if len(buf) >= 11:
                    elems.append(buf[:11])
                    buf = []
    return nodes, elems


def main():
    path = sys.argv[1]
    nodes, elems = parse_inp(path)
    print(f"nodes parsed   : {len(nodes)}")
    print(f"C3D10 elements : {len(elems)}")
    if not nodes or not elems:
        print("ERROR: parse produced nothing; check *NODE/*ELEMENT keywords")
        sys.exit(1)

    ids = np.array(sorted(nodes))
    xyz = np.array([nodes[i] for i in ids])
    idx = {n: k for k, n in enumerate(ids)}

    E = np.array(elems)[:, 1:5]  # 4 corner node IDs
    Ek = np.vectorize(idx.get)(E)  # to row indices

    tri = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    faces = np.concatenate([np.sort(Ek[:, list(t)], axis=1) for t in tri])
    uniq, cnt = np.unique(faces, axis=0, return_counts=True)
    print(f"count histogram (faces shared by n elements): "
          f"{dict(zip(*np.unique(cnt, return_counts=True)))}")
    bf = uniq[cnt == 1]
    print(f"boundary faces : {len(bf)}")

    zmin = xyz[:, 2].min()
    zlo, zhi = zmin, zmin + 1.5
    print(f"mesh z_min = {zmin:.6f}  -> root rib slab z in [{zlo:.4f}, {zhi:.4f}]")

    z = xyz[bf, 2]  # (nf,3)
    on_lo = np.abs(z - zlo) < TOL
    on_hi = np.abs(z - zhi) < TOL
    all_lo = on_lo.all(axis=1)
    all_hi = on_hi.all(axis=1)
    span = on_lo.any(axis=1) & on_hi.any(axis=1)
    inside = ((z > zlo - TOL) & (z < zhi + TOL)).all(axis=1)
    strictly_between = ((z > zlo + 0.1) & (z < zhi - 0.1)).any(axis=1)

    print("\nBoundary faces with all 3 nodes inside the root-rib slab:")
    print(f"  entirely on z_lo plane (flat face) : {int((all_lo).sum())}")
    print(f"  entirely on z_hi plane (flat face) : {int((all_hi).sum())}")
    print(f"  SPANNING lo<->hi (rim candidates)  : {int(span.sum())}")
    print(f"  any node strictly between planes   : {int(strictly_between.sum())}")
    print(f"  slab-contained total               : {int(inside.sum())}")

    rim = bf[span & inside]
    if len(rim):
        rn = np.unique(rim)
        p = xyz[rn]
        cen = xyz[rim].mean(axis=1)
        up = cen[:, 1] > 0
        print(f"\nRim faces: {len(rim)}; unique nodes {len(rn)}")
        print(f"  x range {p[:,0].min():.3f} .. {p[:,0].max():.3f}")
        print(f"  y range {p[:,1].min():.3f} .. {p[:,1].max():.3f}")
        print(f"  rim faces with centroid y>0 (upper half): {int(up.sum())}")
        # rim edge length distribution (edges along the outline, i.e. lo-lo or hi-hi)
        print("\nSample of 6 rim faces (node coords):")
        for f in rim[:6]:
            print("  ", [tuple(np.round(xyz[k], 3)) for k in f])
    else:
        print("\nNo spanning rim faces found.")


if __name__ == "__main__":
    main()
