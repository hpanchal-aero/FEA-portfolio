#!/usr/bin/env python3
"""
Parse the BUCKLING FACTOR OUTPUT table from a CalculiX .dat for a skin-model *BUCKLE job and
report the critical factor against the requirement (>= 1.5 at limit load). The factor multiplies
the entire applied load (lift + gravity x n). No closed-form comparison is made here.

Usage: python3 scripts/extract_stage3_buckle.py <job_path_without_ext> [required_factor=1.5]
"""
import re
import sys

ROW_RE = re.compile(r"^\s*(\d+)\s+([-+]?\d*\.?\d+[EeDd][+-]?\d+)\s*$")
LC_RE = re.compile(r"LOADCASE\s+LIFT_N=([-+0-9.eE]+)\s+NG=([-+0-9.eE]+)")


def main():
    job = sys.argv[1]
    req = float(sys.argv[2]) if len(sys.argv) > 2 else 1.5
    lines = open(job + ".dat").readlines()
    hi = None
    for i, line in enumerate(lines):
        if "BUCKLINGFACTOROUTPUT" in "".join(line.upper().split()):
            hi = i
            break
    if hi is None:
        print("ERROR: 'BUCKLING FACTOR OUTPUT' header not found in .dat"); sys.exit(1)
    rows = []
    for line in lines[hi:hi + 60]:
        m = ROW_RE.match(line)
        if m:
            rows.append((int(m.group(1)), float(m.group(2).replace("D", "E").replace("d", "E"))))
        elif rows and line.strip() == "":
            break
    if not rows:
        print("ERROR: no buckling factor rows parsed")
        for line in lines[hi:hi + 15]:
            print(f"  {line!r}")
        sys.exit(1)
    lc = None
    with open(job + ".inp") as f:
        for _ in range(12):
            m = LC_RE.search(f.readline())
            if m:
                lc = (float(m.group(1)), float(m.group(2)))
                break
    print(f"Job: {job}")
    if lc:
        print(f"Applied load = limit load: lift {lc[0]:.2f} N per half-wing, n = {lc[1]}")
    for mode, fac in rows:
        print(f"  mode {mode}: buckling factor = {fac:.4f}")
    crit = min(f for _, f in rows)
    print(f"\nCritical buckling factor: {crit:.4f}  (requirement >= {req}) -> "
          f"{'PASS' if crit >= req else 'FAIL'}; margin vs requirement {100 * (crit / req - 1):+.1f}%")


if __name__ == "__main__":
    main()
