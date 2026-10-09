#!/usr/bin/env python3
"""
Compare two raw Netgen meshes (.npz written by mesh_stage3_skin_raw.py): per array, shape and exact equality.
Usage: /usr/bin/python3 scripts/compare_raw_meshes.py <a.npz> <b.npz>
"""
import sys
import numpy as np


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
    ka, kb = sorted(a.files), sorted(b.files)
    print(f"arrays: {ka}; same names: {ka == kb}")
    allsame = ka == kb
    for k in ka:
        if k not in kb:
            continue
        x, y = a[k], b[k]
        same = x.shape == y.shape and np.array_equal(x, y)
        allsame = allsame and same
        extra = ""
        if not same:
            extra = (f"; max |diff| {np.abs(x - y).max():.3g}" if x.shape == y.shape
                     else f"; shapes {x.shape} vs {y.shape}")
        print(f"  {k}: shape {x.shape}, identical: {same}{extra}")
    print("RESULT: meshes identical" if allsame else "RESULT: meshes DIFFER")


if __name__ == "__main__":
    main()
