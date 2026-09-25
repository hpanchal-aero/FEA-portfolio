"""
Project 02, Configuration B -- mesh feasibility (counts only, no .inp,
no solve) for randomly generated pocket assemblies, at multiple mesh
sizes.

Usage: /usr/bin/python3 scripts/mesh_config_b_random_feasibility.py <tag|--all> [maxh ...]
  default maxh list: 3.0 2.0 1.8

Reads mesh/config_b_study/config_b_random_assembly_<tag>.step, generates
a linear-tet mesh at each size, reports element/node counts and % of the
documented ceilings (230,000 elements / 1,180,000 equations). No .inp
written. Failures (Netgen exceptions) are caught per-case so one bad
geometry doesn't stop the batch.
"""
import csv
import glob
import os
import sys
import time

from netgen.occ import OCCGeometry

ELEM_CEILING = 230_000
EQ_CEILING = 1_180_000

OUT_DIR = "mesh/config_b_study"
RESULTS_CSV = "results/config_b_random_mesh_feasibility.csv"


def mesh_one(tag, maxh):
    step_file = os.path.join(OUT_DIR, f"config_b_random_assembly_{tag}.step")
    if not os.path.isfile(step_file):
        return {"tag": tag, "maxh": maxh, "status": "step_missing"}
    try:
        geo = OCCGeometry(step_file)
        n_solids = len(list(geo.shape.solids))
        n_faces = len(list(geo.shape.faces))
        t0 = time.time()
        mesh = geo.GenerateMesh(maxh=maxh)
        dt = time.time() - t0
        n_pts = len(mesh.Points())
        n_tet = len(mesh.Elements3D())
        # rough C3D10 node estimate consistent with mesh_config_b.py's
        # actual promoted counts is not computed here (no promotion done);
        # this reports LINEAR mesh counts only, as a fast feasibility proxy.
        return {
            "tag": tag, "maxh": maxh, "status": "ok",
            "n_solids": n_solids, "n_faces": n_faces,
            "linear_nodes": n_pts, "tets": n_tet, "gen_time_s": round(dt, 2),
        }
    except Exception as e:
        return {"tag": tag, "maxh": maxh, "status": f"FAILED: {type(e).__name__}: {e}"}


def main():
    if len(sys.argv) < 2:
        print("Usage: mesh_config_b_random_feasibility.py <tag|--all> [maxh ...]")
        sys.exit(1)
    sizes = [float(a) for a in sys.argv[2:]] or [3.0, 2.0, 1.8]

    if sys.argv[1] == "--all":
        steps = sorted(glob.glob(os.path.join(OUT_DIR, "config_b_random_assembly_*.step")))
        tags = [os.path.basename(s)[len("config_b_random_assembly_"):-len(".step")] for s in steps]
    else:
        tags = [sys.argv[1]]

    print(f"Processing {len(tags)} tag(s) at sizes {sizes}")
    os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
    fieldnames = ["tag", "maxh", "status", "n_solids", "n_faces",
                  "linear_nodes", "tets", "gen_time_s"]
    n_ok, n_fail = 0, 0
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for i, tag in enumerate(tags, 1):
            for maxh in sizes:
                r = mesh_one(tag, maxh)
                w.writerow(r)
                f.flush()
                if r["status"] == "ok":
                    n_ok += 1
                    print(f"[{i}/{len(tags)}] {tag} maxh={maxh}: tets={r['tets']} "
                          f"nodes={r['linear_nodes']} time={r['gen_time_s']}s")
                else:
                    n_fail += 1
                    print(f"[{i}/{len(tags)}] {tag} maxh={maxh}: {r['status']}")

    print(f"\nDone: {n_ok} ok, {n_fail} failed/missing out of {n_ok + n_fail} (tag x size combos)")
    print(f"Saved: {RESULTS_CSV}")


if __name__ == "__main__":
    main()
