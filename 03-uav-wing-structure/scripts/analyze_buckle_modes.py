#!/usr/bin/env python3
"""
Where are the buckling modes?  Parses every DISP block of a *BUCKLE .frd (one per mode), uses the
.frd node coordinates (skin-only original IDs are absent from the .frd; expanded nodes carry the
skin), and for each mode reports the region with |U| > 50% of the peak: count, centroid, z/x
extent (2nd-98th percentile), nearest bay between ribs, upper/lower skin split, proximity to
leading/trailing edge, and the dominant displacement component.

Usage (system or conda Python, numpy only):
  python3 scripts/analyze_buckle_modes.py <job> [<job> ...]
"""
import sys
import numpy as np

RIBS = np.array([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])


def read_modes(path):
    ids, xyz, blocks = [], [], []
    on_nodes, cur = False, None
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                on_nodes = True
                continue
            if on_nodes:
                if line.startswith(" -3"):
                    on_nodes = False
                elif line.startswith(" -1"):
                    ids.append(int(line[3:13]))
                    xyz.append([float(line[13 + 12 * k:25 + 12 * k]) for k in range(3)])
                continue
            if line.startswith(" -4"):
                cur = [] if line[5:13].strip() == "DISP" else None
                if cur is not None:
                    blocks.append(cur)
                continue
            if cur is not None:
                if line.startswith(" -3"):
                    cur = None
                elif line.startswith(" -1"):
                    cur.append((int(line[3:13]),
                                [float(line[13 + 12 * k:25 + 12 * k]) for k in range(3)]))
    return np.array(ids), np.array(xyz), blocks


def main():
    for job in sys.argv[1:]:
        ids, xyz, blocks = read_modes(job + ".frd")
        order = np.argsort(ids)
        ids, xyz = ids[order], xyz[order]
        print(f"\n######## {job}: {len(ids)} frd nodes, {len(blocks)} displacement blocks (modes)")
        chord = 300.0 - 0.15 * xyz[:, 2]
        dTE = 0.7 * chord - xyz[:, 0]
        dLE = xyz[:, 0] + 0.3 * chord
        for k, blk in enumerate(blocks, start=1):
            bid = np.array([b[0] for b in blk])
            U = np.array([b[1] for b in blk])
            pos = np.searchsorted(ids, bid)
            ok = (pos < len(ids)) & (ids[np.minimum(pos, len(ids) - 1)] == bid)
            bid, U, pos = bid[ok], U[ok], pos[ok]
            P = xyz[pos]
            mag = np.linalg.norm(U, axis=1)
            mmax = mag.max()
            hot = mag > 0.5 * mmax
            Ph, Uh = P[hot], U[hot]
            zlo, zhi = np.percentile(Ph[:, 2], [2, 98])
            xlo, xhi = np.percentile(Ph[:, 0], [2, 98])
            zc = Ph[:, 2].mean()
            bay = int(np.clip(np.searchsorted(RIBS, zc) - 1, 0, 4))
            kmax = int(np.argmax(mag))
            comp = np.abs(Uh).mean(axis=0)
            comp = comp / comp.sum()
            print(f"mode {k}: peak |U| {mmax:.3e} (normalised eigenvector) at "
                  f"({P[kmax,0]:.1f},{P[kmax,1]:.1f},{P[kmax,2]:.1f})")
            print(f"   hot nodes (>50% peak): {int(hot.sum())}; z range {zlo:.0f}..{zhi:.0f} "
                  f"(centroid {zc:.0f}, bay between ribs {RIBS[bay]:.0f}-{RIBS[bay+1]:.0f}); "
                  f"x range {xlo:.0f}..{xhi:.0f}")
            print(f"   upper skin (y>0): {100*(Ph[:,1] > 0).mean():.0f}%; within 20 mm of TE: "
                  f"{100*(dTE[pos][hot] < 20).mean():.0f}%; within 20 mm of LE: "
                  f"{100*(dLE[pos][hot] < 20).mean():.0f}%; within 3 mm of a rib plane: "
                  f"{100*(np.abs(Ph[:,2][:, None] - RIBS[None, :]).min(axis=1) < 3).mean():.0f}%; "
                  f"|U| components x/y/z: {comp[0]:.2f}/{comp[1]:.2f}/{comp[2]:.2f}")


if __name__ == "__main__":
    main()
