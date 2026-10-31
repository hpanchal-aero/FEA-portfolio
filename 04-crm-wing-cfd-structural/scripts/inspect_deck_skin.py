#!/usr/bin/env python3
"""Read-only: identify the skin nodes of the structural deck.
usage: inspect_deck_skin.py <deck.inp>
A group (ELSET) is treated as skin-like if >= 90% of its nodes with y >= 3.1 m are dataset-loaded (*CLOAD) nodes.
"""
import os
import sys
import collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import inspect_deck_loads as idl


def parse_elements(path):
    groups = collections.OrderedDict()
    section = None
    cur = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                parts = [p.strip() for p in line.split(",")]
                section = parts[0].upper()
                cur = None
                if section == "*ELEMENT":
                    typ = None
                    name = None
                    for p in parts[1:]:
                        u = p.upper()
                        if u.startswith("TYPE="):
                            typ = p.split("=", 1)[1]
                        elif u.startswith("ELSET="):
                            name = p.split("=", 1)[1]
                    cur = name if name else "ELEM_%d" % len(groups)
                    groups[cur] = {"type": typ, "elems": 0, "nodes": set()}
                continue
            if section == "*ELEMENT" and cur is not None:
                tok = [t.strip() for t in line.split(",") if t.strip() != ""]
                groups[cur]["elems"] += 1
                for t in tok[1:]:
                    groups[cur]["nodes"].add(int(t))
    return groups


def main(path):
    nodes, block_of, cload, _, boundary = idl.parse(path)
    loaded = set(i for i in cload if i in nodes)
    groups = parse_elements(path)
    constrained = set(k for k, _, _ in boundary if isinstance(k, int))
    print("deck:", path)
    print("element groups:", len(groups), " types:", dict(collections.Counter(g["type"] for g in groups.values())))

    rows = []
    for name, g in groups.items():
        if not (g["type"] or "").upper().startswith("S"):
            continue
        ns = [n for n in g["nodes"] if n in nodes]
        out = [n for n in ns if nodes[n][1] >= 3.1]
        frac = (sum(1 for n in out if n in loaded) / len(out)) if out else float("nan")
        rows.append((name, g["type"], g["elems"], len(ns), len(out), frac))
    fr = np.array([r[5] for r in rows])
    print("\nshell groups:", len(rows))
    bins = [(0.0, 0.05), (0.05, 0.5), (0.5, 0.9), (0.9, 1.0001)]
    for a, b in bins:
        k = int(np.sum((fr >= a) & (fr < b)))
        print("  loaded fraction in [%.2f, %.2f): %d groups" % (a, b, k))
    print("  groups with no nodes at y >= 3.1 m (fraction undefined):", int(np.sum(np.isnan(fr))))
    print("\nlargest 15 shell groups by node count: name type elems nodes nodes(y>=3.1) loaded-fraction")
    for r in sorted(rows, key=lambda r: -r[3])[:15]:
        print("  %-8s %-5s %6d %6d %6d   %.3f" % r)

    skin = set()
    nskin_groups = 0
    for r in rows:
        if r[5] == r[5] and r[5] >= 0.9:
            skin |= groups[r[0]]["nodes"]
            nskin_groups += 1
    skin = set(n for n in skin if n in nodes)
    print("\nskin-like groups (fraction >= 0.9):", nskin_groups, "  union nodes:", len(skin))
    skin_out = set(n for n in skin if nodes[n][1] >= 2.94)
    print("skin nodes with y >= 2.94 m:", len(skin_out), "  dataset-loaded nodes:", len(loaded))
    print("  loaded but not in skin union:", len(loaded - skin), "   skin (y>=2.94) but not loaded:", len(skin_out - loaded))
    stub = set(n for n in skin if nodes[n][1] < 2.94)
    print("\nstub skin nodes (y < 2.94 m):", len(stub))
    if stub:
        ys = np.array([nodes[n][1] for n in stub])
        print("  y range %.3f..%.3f" % (ys.min(), ys.max()))
        print("  of which constrained (any DOF):", len(stub & constrained))
        print("  of which at y < 0.01 m (clamp plane):", int(np.sum(ys < 0.01)))
        free = [n for n in stub if n not in constrained]
        print("  free stub skin nodes (not constrained):", len(free))
        if free:
            yf = np.array([nodes[n][1] for n in free])
            print("  free stub y range %.3f..%.3f" % (yf.min(), yf.max()))
    all_stub_nodes = [n for n in nodes if block_of[n] == "NALL" and nodes[n][1] < 3.1]
    print("all NALL nodes with y < 3.1 m:", len(all_stub_nodes))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
