#!/usr/bin/env python3
"""Compare CalculiX eigenfrequencies (.dat) with the NASTRAN SOL 103 reference (coarse mesh).

usage: post_modal.py <job.dat> [tolerance_percent=5] [n_check=5]

Mode matching is by order only: a mode swap or a spurious soft mode shifts the list, so eigenvalues
below 10 (rad/s)^2 are flagged as suspect.
"""
import sys
from pathlib import Path

REF = [1.629275, 2.631756, 3.816032, 4.683773, 6.315011,
       8.112136, 10.59417, 12.84292, 13.97464, 18.96142]

path = Path(sys.argv[1])
tol = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
n_check = int(sys.argv[3]) if len(sys.argv) > 3 else 5

rows, started = [], False
for line in path.read_text().splitlines():
    if "E I G E N V A L U E" in line:
        started = True
        continue
    if not started:
        continue
    t = line.split()
    if len(t) == 5 and t[0].isdigit():
        try:
            rows.append((int(t[0]), float(t[1]), float(t[3])))
        except ValueError:
            pass
    elif rows and not line.strip():
        break
if not rows:
    sys.exit("ABORT: no eigenvalue table found")

print(f"{'mode':>4} {'eigenvalue':>12} {'f_ccx [Hz]':>11} {'f_nastran':>10} {'ratio':>7} {'err %':>7}  flag")
n_suspect = 0
worst = 0.0
for k, ev, f in rows:
    ref = REF[k - 1] if k <= len(REF) else None
    flag = "SUSPECT soft mode" if ev < 10.0 else ""
    n_suspect += ev < 10.0
    if ref:
        ratio = f / ref
        err = 100.0 * (ratio - 1.0)
        if k <= n_check:
            worst = max(worst, abs(err))
        print(f"{k:>4} {ev:>12.5e} {f:>11.5f} {ref:>10.5f} {ratio:>7.4f} {err:>+7.2f}  {flag}")
    else:
        print(f"{k:>4} {ev:>12.5e} {f:>11.5f} {'-':>10} {'-':>7} {'-':>7}  {flag}")

print()
if n_suspect:
    print(f"RESULT: INVALID - {n_suspect} suspect soft mode(s); the comparison by mode order is not meaningful")
else:
    verdict = "PASS" if worst <= tol else "FAIL"
    print(f"RESULT: {verdict} - worst error over first {n_check} modes {worst:.2f} % (tolerance {tol:.1f} %)")
