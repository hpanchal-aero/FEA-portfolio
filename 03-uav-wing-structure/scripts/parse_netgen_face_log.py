#!/usr/bin/env python3
"""
Summarise a Netgen verbose log written by `diag_netgen_faces.py verbose`: which faces were started and finished
('Face k / n ...' and 'Surface meshing done'), where ERROR lines sit, and where the volume-meshing progress lines
('faces = N trials = N elements = N') begin.
Usage: python3 scripts/parse_netgen_face_log.py <log>
"""
import re
import sys

FACE = re.compile(r"Face (\d+) / (\d+)")
TRIALS = re.compile(r"faces = (\d+) trials = (\d+) elements = (\d+)")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    lines = open(sys.argv[1], errors="replace").read().splitlines()
    started, done_set, errors = [], set(), []
    cur, nfaces = None, None
    trials = []
    done_order = []
    for ln, line in enumerate(lines, 1):
        m = FACE.search(line)
        if m:
            cur, nfaces = int(m.group(1)), int(m.group(2))
            started.append((ln, cur))
            continue
        if "Surface meshing done" in line:
            if cur is not None:
                done_set.add(cur)
                done_order.append(cur)
            continue
        if "ERROR" in line:
            errors.append((ln, cur))
            continue
        t = TRIALS.search(line)
        if t:
            trials.append((ln, int(t.group(1)), int(t.group(2)), int(t.group(3))))
    ks = [k for _, k in started]
    print(f"{sys.argv[1]}: {len(lines)} lines; faces declared {nfaces}; 'Face k/n' lines {len(started)} "
          f"(distinct {len(set(ks))}); 'Surface meshing done' {len(done_order)} (distinct faces {len(done_set)})")
    print(f"face order strictly increasing: {all(a < b for a, b in zip(ks, ks[1:]))}; "
          f"last face started: {ks[-1] if ks else None} at line {started[-1][0] if started else None}")
    print(f"ERROR lines (line, face in progress): {errors[:10]}")
    for ln, k in errors[:2]:
        print(f"--- context around line {ln}:")
        for j in range(max(0, ln - 8), min(len(lines), ln + 5)):
            print(f"{j + 1:6d}: {lines[j][:140]}")
    not_done = [k for _, k in started if k not in done_set]
    print(f"faces started but not done: {len(not_done)}; first 20: {not_done[:20]}")
    if nfaces:
        never = sorted(set(range(1, nfaces + 1)) - set(ks))
        print(f"faces never started: {len(never)}; first 10: {never[:10]}; last 10: {never[-10:]}")
    if trials:
        print(f"volume-mesher progress lines: {len(trials)}; first (line, faces, trials, elements) {trials[0]}; "
              f"last {trials[-1]}")
        if started:
            print(f"first progress line {trials[0][0]} vs last face line {started[-1][0]}")
    else:
        print("no volume-mesher progress lines")


if __name__ == "__main__":
    main()
