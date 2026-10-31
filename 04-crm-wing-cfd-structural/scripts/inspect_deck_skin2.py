#!/usr/bin/env python3
"""Read-only: classify skin nodes by shell element orientation (|n_z| >= 0.7) and compare with the dataset-loaded nodes.
usage: inspect_deck_skin2.py <deck.inp> [nz_threshold]
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import inspect_deck_loads as idl


def parse_shells(path):
    elems = []
    groups = []
    section = None
    cur = None
    is_shell = False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                parts = [p.strip() for p in line.split(",")]
                section = parts[0].upper()
                is_shell = False
                cur = None
                if section == "*ELEMENT":
                    typ = ""
                    name = ""
                    for p in parts[1:]:
                        u = p.upper()
                        if u.startswith("TYPE="):
                            typ = p.split("=", 1)[1].upper()
                        elif u.startswith("ELSET="):
                            name = p.split("=", 1)[1]
                    is_shell = typ.startswith("S4")
                    cur = name
                continue
            if section == "*ELEMENT" and is_shell:
                tok = [t.strip() for t in line.split(",") if t.strip() != ""]
                if len(tok) >= 5:
                    elems.append([int(t) for t in tok[1:5]])
                    groups.append(cur)
    return np.array(elems, dtype=np.int64), groups


def main(path, thr):
    nodes, block_of, cload, _, boundary = idl.parse(path)
    loaded = set(i for i in cload if i in nodes)
    constrained = set(k for k, _, _ in boundary if isinstance(k, int))
    E, groups = parse_shells(path)
    print("deck:", path)
    print("S4 elements:", len(E))
    X = np.array([[nodes[n] for n in e] for e in E])  # (ne, 4, 3)
    d1 = X[:, 2, :] - X[:, 0, :]
    d2 = X[:, 3, :] - X[:, 1, :]
    nv = np.cross(d1, d2)
    nl = np.linalg.norm(nv, axis=1)
    ok = nl > 0
    nz = np.zeros(len(E))
    nz[ok] = np.abs(nv[ok, 2]) / nl[ok]
    print("degenerate (zero-area) elements:", int((~ok).sum()))
    h, edges = np.histogram(nz, bins=10, range=(0.0, 1.0))
    print("histogram of |n_z| over shell elements:")
    for a, b, c in zip(edges[:-1], edges[1:], h):
        print("  [%.1f, %.1f): %d" % (a, b, c))

    skin_el = nz >= thr
    print("\nskin elements (|n_z| >= %.2f): %d of %d" % (thr, skin_el.sum(), len(E)))
    skin_nodes = set(int(n) for n in E[skin_el].ravel())
    web_nodes = set(int(n) for n in E[~skin_el].ravel())
    print("skin nodes:", len(skin_nodes), "  web/other-element nodes:", len(web_nodes), "  nodes in both:", len(skin_nodes & web_nodes))

    out = set(n for n in skin_nodes if nodes[n][1] >= 2.94)
    print("\n--- outboard check (y >= 2.94 m) ---")
    print("skin nodes y >= 2.94:", len(out), "  dataset-loaded:", len(loaded))
    print("  loaded but not skin:", len(loaded - skin_nodes), "   skin (y>=2.94) but not loaded:", len(out - loaded))
    miss = sorted(loaded - skin_nodes)[:5]
    for n in miss:
        print("   loaded-not-skin example node", n, nodes[n])
    extra = sorted(out - loaded)[:5]
    for n in extra:
        print("   skin-not-loaded example node", n, nodes[n])

    stub = set(n for n in skin_nodes if nodes[n][1] < 2.94)
    print("\n--- stub (y < 2.94 m) ---")
    print("stub skin nodes:", len(stub))
    if stub:
        ys = np.array([nodes[n][1] for n in stub])
        print("  y range %.3f..%.3f" % (ys.min(), ys.max()))
        print("  constrained (any DOF):", len(stub & constrained), "  at y < 0.01 m:", int((ys < 0.01).sum()))
        free = stub - constrained
        print("  free stub skin nodes:", len(free))
        yy = np.array([nodes[n][1] for n in free])
        print("  free stub y range %.3f..%.3f" % (yy.min(), yy.max()))
        for lo, hi in ((0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 2.94)):
            print("    free stub skin nodes with y in [%.1f, %.1f): %d" % (lo, hi, int(((yy >= lo) & (yy < hi)).sum())))
    print("all NALL nodes with y < 2.94 m:", sum(1 for n in nodes if block_of[n] == "NALL" and nodes[n][1] < 2.94))


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) == 3 else 0.7)
