"""
Parses the BUCKLING FACTOR OUTPUT table from a CalculiX .ccx.log for a
*BUCKLE step and compares the lowest (critical) factor against the
closed-form target of 1.0 (since the reference load applied IS the
closed-form P_cr).

Usage: python3 scripts/buckle_verify_extract.py <job_path_without_ext>
"""
import re
import sys

ROW_RE = re.compile(r"^\s*(\d+)\s+([-+]?\d*\.?\d+[EeDd][+-]?\d+)\s*$")


def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: buckle_verify_extract.py <job> [closed_form_target_factor]"); sys.exit(1)
    job = sys.argv[1]
    # The closed-form comparison is ONLY meaningful when the applied reference load IS the
    # closed-form critical load (flat-plate verification, target 1.0). Pass the target explicitly.
    target = float(sys.argv[2]) if len(sys.argv) == 3 else None
    lines = open(job + ".dat").readlines()

    header_i = None
    for i, line in enumerate(lines):
        if "BUCKLINGFACTOROUTPUT" in "".join(line.upper().split()):
            header_i = i
            break
    if header_i is None:
        print("ERROR: 'BUCKLING FACTOR OUTPUT' header not found in log")
        sys.exit(1)

    rows = []
    for line in lines[header_i:header_i + 40]:
        m = ROW_RE.match(line)
        if m:
            rows.append((int(m.group(1)), float(m.group(2).replace("D", "E").replace("d", "E"))))
        elif rows and line.strip() == "":
            break  # blank line after the table ends it

    if not rows:
        print("ERROR: no buckling factor rows parsed after header")
        print("Lines examined:")
        for line in lines[header_i:header_i + 15]:
            print(f"  {line!r}")
        sys.exit(1)

    print(f"Job: {job}")
    for mode, factor in rows:
        print(f"  mode {mode}: buckling factor = {factor:.6f}")

    crit_mode, crit_factor = min(rows, key=lambda r: r[1])
    print(f"\nCritical mode: {crit_mode}, factor={crit_factor:.6f}")
    if target is not None:
        dev = 100 * (crit_factor - target) / target
        print(f"Closed-form target factor: {target:.6f}; deviation: {dev:+.2f}%")
    else:
        print("(no closed-form target given: no theory comparison printed)")


if __name__ == "__main__":
    main()
