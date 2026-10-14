#!/usr/bin/env python3
"""
Project 03, Stage 6: resumable two-pass layout study driver.

  pass1 : per layout -> geometry (build_stage3_solid) -> imprint (build_stage3_skin_geom) -> Netgen mesh
          -> CalculiX input (build_stage3_skin_inp) -> Stage 4 pressure load (apply_stage4_pressure_load).
          All geometry/mesh gates run here, before any solver time is spent.
  pass2 : per layout -> static solve -> extract (tip deflection, mass, stresses) -> buckling solve ->
          buckling factors -> cleanup of the large solver files.

Monitoring (second terminal):
  tail -f results/stage6/driver.log     concise: step start/end, gate numbers, heartbeat, ETA
  tail -f results/stage6/live.log       every output line of every child process, prefixed [layout step]
Full per-step logs: results/stage6/logs/<key>_<step>.log ; per-layout records: results/stage6/records/<key>.json

Usage (system python, project directory):
  scripts/run_layout_study.py --stage check
  scripts/run_layout_study.py --stage pass1 --layouts L000
  scripts/run_layout_study.py --stage pass2 --layouts L000
  scripts/run_layout_study.py --stage pass1 --layouts all
  scripts/run_layout_study.py --stage summary
Options: --rib-t <mm> (default 1.5), --retry-failed, --max-consecutive-fail N (default 3)

Provisional gate thresholds (chosen from the baseline imprint distribution; the minima are recorded per layout):
  smallest free imprint piece >= 5 mm^2 ; smallest shared piece >= 0.5 mm^2.
Solid 4 mm / skin 8 mm mesh, NETGEN_KW="closeedgefac=None", skin t = 1.0 mm, limit load 372.78 N, n = 3.8.
"""
import argparse
import csv
import fcntl
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(PROJECT)
sys.path.insert(0, os.path.join(PROJECT, "scripts"))
import variant_tag  # noqa: E402

SYS_PY = "/usr/bin/python3"
CONDA_ENV = "/home/harsh/miniconda3/envs/ccx223"
CONDA_PY = "/usr/bin/python3"   # the ccx223 env has no python of its own; system python has numpy/scipy/pyvista
CONDA_SH = "/home/harsh/miniconda3/etc/profile.d/conda.sh"
LAYOUT_CSV = "results/stage6_layouts.csv"
RES = "results/stage6"
LOGS = RES + "/logs"
RECS = RES + "/records"
STATIC_DIR = RES + "/static"

RSS_LIMIT_GB = 8
MEM_NEED_GB = 7.0
MIN_FREE_AREA = 5.0
MIN_SHARED_AREA = 0.5
MOMENT_TOL_PCT = 0.5
SPAR_V2_VOLUME = 89910.4455          # mm^3, spar v2 alone (Stage 3)
RHO_G_MM3 = 2.78e-3                  # 2780 kg/m^3
DEFAULT_RIB_T = 1.5
SKIN_T = 1.0          # skin thickness [mm]; changed by --skin-t (stage "skin" only)
SKIN_SUFFIX = ""     # "" at 1.0 mm, "_sk0p8" otherwise (record keys, summary names)


def skin_p():
    return ("%g" % SKIN_T).replace(".", "p")
TIMEOUT = dict(geom=600, imprint=900, mesh=150, inp=600, load=900, solve=1200, buckle=1800, extract=600, small=300)

BASE_TIP = 4.197838931623932         # Stage 4 baseline, mm
BASE_BUCKLE = 3.1938                 # Stage 4 baseline (4 decimals)
BASE_MASS = 1.6352093547873905       # kg


class Log:
    def __init__(self):
        os.makedirs(LOGS, exist_ok=True)
        os.makedirs(RECS, exist_ok=True)
        os.makedirs(STATIC_DIR, exist_ok=True)
        self.drv = open(RES + "/driver.log", "a", buffering=1)
        self.live = open(RES + "/live.log", "a", buffering=1)
        self.lock = threading.Lock()

    def say(self, msg):
        line = "%s  %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
        with self.lock:
            self.drv.write(line + "\n")
            self.live.write(line + "\n")
        print(line, flush=True)

    def raw(self, prefix, text):
        with self.lock:
            self.live.write(prefix + text + "\n")


LOG = None


def free_disk_gb():
    """Smaller of the WSL virtual-disk free space and the real Windows C: free space."""
    v = shutil.disk_usage(PROJECT).free / 1e9
    try:
        v = min(v, shutil.disk_usage("/mnt/c").free / 1e9)
    except OSError:
        pass
    return v


def mem_avail_gb():
    for ln in open("/proc/meminfo"):
        if ln.startswith("MemAvailable:"):
            return int(ln.split()[1]) / 1048576.0
    return 0.0


def kill_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    for _ in range(10):
        if proc.poll() is not None:
            break
        time.sleep(0.5)
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run(cmd, label, c, timeout, env=None):
    """Run a child process in its own process group; tee every line to the step log and live.log."""
    t0 = time.time()
    logpath = "%s/%s_%s.log" % (LOGS, c.key, label)
    shown = cmd if isinstance(cmd, str) else " ".join(cmd)
    LOG.say("[%s] %s: START  %s" % (c.key, label, shown))
    proc = subprocess.Popen(cmd, cwd=PROJECT, env=env if env is not None else c.env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, errors="replace", bufsize=1,
                            start_new_session=True)
    lines = []
    f = open(logpath, "w")

    def reader():
        for ln in proc.stdout:
            ln = ln.rstrip("\n")
            lines.append(ln)
            f.write(ln + "\n")
            f.flush()
            LOG.raw("[%s %s] " % (c.key, label), ln)

    th = threading.Thread(target=reader, daemon=True)
    th.start()
    timed_out = False
    last_hb = 0.0
    try:
        while True:
            try:
                proc.wait(timeout=5)
                break
            except subprocess.TimeoutExpired:
                el = time.time() - t0
                if el > timeout:
                    timed_out = True
                    LOG.say("[%s] %s: TIMEOUT after %.0f s -> killing the process group" % (c.key, label, el))
                    kill_group(proc)
                    proc.wait()
                    break
                if el - last_hb >= 60:
                    last_hb = el
                    LOG.say("[%s] %s: running %.0f s (limit %d s), free memory %.1f GB"
                            % (c.key, label, el, timeout, mem_avail_gb()))
    except BaseException:
        kill_group(proc)
        raise
    th.join(timeout=10)
    f.close()
    wall = time.time() - t0
    LOG.say("[%s] %s: END rc=%s wall=%.1f s%s" % (c.key, label, proc.returncode, wall, "  TIMED OUT" if timed_out else ""))
    return dict(rc=proc.returncode, wall=wall, text="\n".join(lines), timed_out=timed_out, log=logpath)


def grab(pattern, text, cast=float, group=1):
    m = re.search(pattern, text)
    return cast(m.group(group)) if m else None


def compute_tag(L, rib_t):
    keys = ("RIB_SPACING", "RIB_STATIONS", "LAYOUT_ID", "RIB_T")
    saved = {k: os.environ.get(k) for k in keys}
    try:
        os.environ.pop("RIB_SPACING", None)
        os.environ["RIB_STATIONS"] = L["stations"]
        os.environ["LAYOUT_ID"] = L["id"]
        os.environ["RIB_T"] = "%g" % rib_t
        return variant_tag.configuration(1000.0)[3]
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def layout_env(L, rib_t):
    e = os.environ.copy()
    for k in ("RIB_SPACING", "NETGEN_KW", "RIB_STATIONS", "LAYOUT_ID", "RIB_T"):
        e.pop(k, None)
    e["RIB_STATIONS"] = L["stations"]
    e["LAYOUT_ID"] = L["id"]
    e["RIB_T"] = "%g" % rib_t
    return e


def names(tag):
    n = {}
    n["step"] = "mesh/stage3_solid/stage3_spar_ribs_v2%s.step" % tag
    n["json"] = "mesh/stage3_skin/stage3_skin_faces%s.json" % tag
    n["vol"] = "mesh/stage3_skin/stage3_solid_volume%s.txt" % tag
    n["raw"] = "mesh/stage3_skin/raw_s4_k8%s_closeedgefacnone.npz" % tag
    n["static_job"] = "simulation/stage3_skin_s4_k8%s_closeedgefacnone_t%s_n3p8" % (tag, skin_p())
    n["job"] = "simulation/stage4_s4_k8%s_closeedgefacnone_t%s_n3p8" % (tag, skin_p())
    n["viewer"] = "mesh/stage3_skin/stage3_skin_s4_k8%s_closeedgefacnone_t%s_n3p8_viewer.npz" % (tag, skin_p())
    n["static_inp"] = n["static_job"] + ".inp"
    n["job_inp"] = n["job"] + ".inp"
    return n


class Ctx:
    def __init__(self, L, rib_t):
        self.L = L
        self.rib_t = rib_t
        self.tag = compute_tag(L, rib_t)
        self.key = self.tag.lstrip("_") + SKIN_SUFFIX
        self.env = layout_env(L, rib_t)
        self.N = names(self.tag)
        self.rec = None


# ------------------------------------------------------------------ records
def rec_path(key):
    return "%s/%s.json" % (RECS, key)


def save_rec(c):
    p = rec_path(c.key)
    with open(p + ".tmp", "w") as f:
        json.dump(c.rec, f, indent=2)
    os.replace(p + ".tmp", p)


def load_rec(c):
    p = rec_path(c.key)
    if os.path.exists(p):
        d = json.load(open(p))
        if d.get("stations") == c.L["stations"] and d.get("rib_t") == c.rib_t and d.get("skin_t", 1.0) == SKIN_T:
            return d
        LOG.say("[%s] stored record does not match the layout table -> moved to .stale" % c.key)
        os.replace(p, p + ".stale")
    return dict(id=c.L["id"], key=c.key, family=c.L["family"], stations=c.L["stations"], rib_t=c.rib_t,
                tag=c.tag, skin_t=SKIN_T, steps={}, metrics={}, failed_step=None, failed_reason=None, complete=False)


# ------------------------------------------------------------------ pass 1 steps
def step_geom(c):
    t0 = time.time()
    r = run([SYS_PY, "-u", "scripts/build_stage3_solid.py", "ribs"], "geom", c, TIMEOUT["geom"])
    txt = r["text"]
    ok = (r["rc"] == 0 and not r["timed_out"] and "PASS" in txt and "FAIL" not in txt
          and os.path.exists(c.N["step"]) and os.path.getmtime(c.N["step"]) >= t0 - 1)
    reason = "" if ok else "rc=%s timed_out=%s has_PASS=%s has_FAIL=%s step_file=%s; see %s" % (
        r["rc"], r["timed_out"], "PASS" in txt, "FAIL" in txt, os.path.exists(c.N["step"]), r["log"])
    return ok, {"wall_s": round(r["wall"], 1)}, reason


def step_imprint(c):
    t0 = time.time()
    r = run([SYS_PY, "-u", "scripts/build_stage3_skin_geom.py"], "imprint", c, TIMEOUT["imprint"])
    txt = r["text"]
    info = {"wall_s": round(r["wall"], 1)}
    ok = (r["rc"] == 0 and not r["timed_out"] and "OVERALL VALIDATION: PASS" in txt
          and os.path.exists(c.N["json"]) and os.path.getmtime(c.N["json"]) >= t0 - 1)
    if not ok:
        return False, info, "rc=%s timed_out=%s validation_pass=%s; see %s" % (
            r["rc"], r["timed_out"], "OVERALL VALIDATION: PASS" in txt, r["log"])
    m = re.search(r"Skin pieces: (\d+)\s+\(shared with solid: (\d+), free: (\d+)\)", txt)
    if m:
        info.update(pieces=int(m.group(1)), shared=int(m.group(2)), free=int(m.group(3)))
    m2 = re.search(r"Imported solid volume: ([0-9.]+) mm\^3 \(expect ([0-9.]+)\)", txt)
    if m2:
        info.update(volume=float(m2.group(1)), volume_expected=float(m2.group(2)))
    pieces = json.load(open(c.N["json"]))
    fr = sorted((p for p in pieces if not p["shared"]), key=lambda p: p["area"])
    sh = sorted((p for p in pieces if p["shared"]), key=lambda p: p["area"])
    info.update(min_free_area=fr[0]["area"], min_shared_area=sh[0]["area"], n_json=len(pieces))
    LOG.say("[%s] imprint: pieces=%s shared=%s free=%s; min free piece %.4f mm2 (z=%.1f), min shared piece %.4f mm2 (z=%.1f)"
            % (c.key, info.get("pieces"), info.get("shared"), info.get("free"), fr[0]["area"], fr[0]["cz"],
               sh[0]["area"], sh[0]["cz"]))
    if fr[0]["area"] < MIN_FREE_AREA:
        return False, info, "sliver gate: smallest free piece %.4f mm2 < %.1f at cz=%.2f cx=%.2f" % (
            fr[0]["area"], MIN_FREE_AREA, fr[0]["cz"], fr[0]["cx"])
    if sh[0]["area"] < MIN_SHARED_AREA:
        return False, info, "sliver gate: smallest shared piece %.4f mm2 < %.1f at cz=%.2f cx=%.2f" % (
            sh[0]["area"], MIN_SHARED_AREA, sh[0]["cz"], sh[0]["cx"])
    return True, info, ""


def step_mesh(c):
    info = {"attempts": []}
    txt = ""
    for attempt in (1, 2):
        if os.path.exists(c.N["raw"]):
            os.remove(c.N["raw"])
        env = dict(c.env)
        env["NETGEN_KW"] = "closeedgefac=None"
        r = run([SYS_PY, "-u", "scripts/mesh_stage3_skin_raw.py", "4.0", "8.0"], "mesh-try%d" % attempt, c,
                TIMEOUT["mesh"], env=env)
        txt = r["text"]
        has_err = "ERROR" in txt
        good = (r["rc"] == 0 and not r["timed_out"] and "GenerateMesh OK" in txt and not has_err
                and os.path.exists(c.N["raw"]))
        info["attempts"].append(dict(rc=r["rc"], timed_out=r["timed_out"], error_line=has_err,
                                     wall_s=round(r["wall"], 1)))
        if good:
            break
        LOG.say("[%s] mesh attempt %d FAILED (rc=%s timed_out=%s ERROR line=%s)" % (
            c.key, attempt, r["rc"], r["timed_out"], has_err))
    else:
        return False, info, "Netgen failed twice (see logs/%s_mesh-try*.log)" % c.key
    if len(info["attempts"]) > 1:
        LOG.say("[%s] NOTE: mesh needed %d attempts (recorded)" % (c.key, len(info["attempts"])))
    m = re.search(r"Points (\d+); tets (\d+); surface tris (\d+)", txt)
    if m:
        info.update(points=int(m.group(1)), tets=int(m.group(2)), tris=int(m.group(3)))
    orphan = grab(r"without a surface tri: (\d+)", txt, int)
    info["tet_faces_without_tri"] = orphan
    mf = re.search(r"free skin (\d+) -> maxh [0-9.]+; shared skin (\d+) \+ non-skin (\d+)", txt)
    imp = c.rec["steps"].get("imprint", {})
    if mf:
        info.update(faces_free=int(mf.group(1)), faces_shared=int(mf.group(2)), faces_nonskin=int(mf.group(3)))
        if imp and (int(mf.group(1)) != imp.get("free") or int(mf.group(2)) != imp.get("shared")):
            return False, info, "mesh face classes (free %s, shared %s) differ from the imprint (free %s, shared %s)" % (
                mf.group(1), mf.group(2), imp.get("free"), imp.get("shared"))
    else:
        return False, info, "could not parse the face-class line of the mesh log"
    if orphan != 0:
        return False, info, "tet-boundary faces without a surface tri = %s (must be 0)" % orphan
    LOG.say("[%s] mesh: %s tets, %s surface tris, %s points" % (c.key, info.get("tets"), info.get("tris"), info.get("points")))
    return True, info, ""


def step_inp(c):
    t0 = time.time()
    r = run([CONDA_PY, "-u", "scripts/build_stage3_skin_inp.py", c.N["raw"], "%g" % SKIN_T, "372.78", "3.8"], "inp", c,
            TIMEOUT["inp"])
    txt = r["text"]
    info = {"wall_s": round(r["wall"], 1)}
    if r["rc"] != 0 or r["timed_out"] or not os.path.exists(c.N["static_inp"]) or os.path.getmtime(c.N["static_inp"]) < t0 - 1:
        return False, info, "rc=%s timed_out=%s static_inp_exists=%s; see %s" % (
            r["rc"], r["timed_out"], os.path.exists(c.N["static_inp"]), r["log"])
    m = re.search(r"Matched (\d+)/(\d+) skin pieces", txt)
    mt = re.search(r"Tets (\d+): min volume ([0-9.eE+-]+) \(bad (\d+)\)", txt)
    me = re.search(r"Elements: (\d+) C3D10 \+ (\d+) S6 = (\d+) \(([0-9.]+)% of ceiling\)", txt)
    mq = re.search(r"Equation estimate: (\d+) \(([0-9.]+)% of ceiling", txt)
    mr = re.search(r"Root face z_min=([-0-9.]+)", txt)
    mn = re.search(r"of which NOT solid nodes: (\d+)", txt)
    ml = re.search(r"load moment ([0-9.]+) vs ([0-9.]+) N\*mm \(([-+0-9.]+)%\)", txt)
    if not (m and mt and me and mq and mr and mn and ml):
        return False, info, "could not parse all gate lines of the input builder (see %s)" % r["log"]
    info.update(matched=int(m.group(1)), matched_of=int(m.group(2)), tets=int(mt.group(1)), tet_min_volume=float(mt.group(2)),
                tets_bad=int(mt.group(3)), c3d10=int(me.group(1)), s6=int(me.group(2)), elements=int(me.group(3)),
                elem_pct=float(me.group(4)), eq_estimate=int(mq.group(1)), eq_pct=float(mq.group(2)),
                root_z=float(mr.group(1)), load_nodes_not_solid=int(mn.group(1)), load_moment_err_pct=float(ml.group(3)))
    LOG.say("[%s] input: %s C3D10 + %s S6, equation estimate %s (%.1f%% of ceiling), load moment %+.4f%%" % (
        c.key, info["c3d10"], info["s6"], info["eq_estimate"], info["eq_pct"], info["load_moment_err_pct"]))
    if info["matched"] != info["matched_of"]:
        return False, info, "skin pieces matched %d/%d" % (info["matched"], info["matched_of"])
    if info["tets_bad"] != 0:
        return False, info, "%d bad tets" % info["tets_bad"]
    if abs(info["root_z"] + c.rib_t / 2.0) > 1e-6:
        return False, info, "root face z_min %.6f != %.6f" % (info["root_z"], -c.rib_t / 2.0)
    if info["load_nodes_not_solid"] != 0:
        return False, info, "%d load nodes are not solid nodes" % info["load_nodes_not_solid"]
    if abs(info["load_moment_err_pct"]) > MOMENT_TOL_PCT:
        return False, info, "load moment error %+.4f%% beyond %.1f%%" % (info["load_moment_err_pct"], MOMENT_TOL_PCT)
    if info["eq_pct"] > 100.0:
        return False, info, "equation estimate %.1f%% of the ceiling" % info["eq_pct"]
    return True, info, ""


def step_load(c):
    t0 = time.time()
    r = run([SYS_PY, "-u", "scripts/apply_stage4_pressure_load.py", c.N["static_inp"], c.N["job_inp"], "full"],
            "load", c, TIMEOUT["load"])
    txt = r["text"]
    info = {"wall_s": round(r["wall"], 1)}
    ok = (r["rc"] == 0 and not r["timed_out"] and "GATE FAIL" not in txt and "ABORT" not in txt
          and ("written %s" % c.N["job_inp"]) in txt and os.path.exists(c.N["job_inp"])
          and os.path.getmtime(c.N["job_inp"]) >= t0 - 1)
    info["scale"] = grab(r"-> scale to target = ([0-9.]+)", txt)
    info["root_moment_err_pct"] = grab(r"root moment sum\(Fy\*z\) = .* err ([-+0-9.]+)%", txt)
    info["torque_err_pct"] = grab(r"analytic \(CP at 25% chord\) .* err ([-+0-9.]+)%", txt)
    if not ok:
        return False, info, "rc=%s timed_out=%s gate_fail=%s abort=%s written=%s; see %s" % (
            r["rc"], r["timed_out"], "GATE FAIL" in txt, "ABORT" in txt, os.path.exists(c.N["job_inp"]), r["log"])
    LOG.say("[%s] load: scale %s, root-moment err %s%%, torque err %s%%" % (
        c.key, info["scale"], info["root_moment_err_pct"], info["torque_err_pct"]))
    for p in (c.N["static_inp"], c.N["viewer"]):
        if os.path.exists(p):
            os.remove(p)
    return True, info, ""


# ------------------------------------------------------------------ pass 2 steps
def wait_for_memory(c, label):
    t0 = time.time()
    n = 0
    while mem_avail_gb() < MEM_NEED_GB:
        if time.time() - t0 > 1800:
            return False
        if n % 10 == 0:
            LOG.say("[%s] %s: waiting for memory (%.1f GB free, need %.1f GB)" % (c.key, label, mem_avail_gb(), MEM_NEED_GB))
        n += 1
        time.sleep(30)
    return True


def solve(c, job, label, timeout):
    if not wait_for_memory(c, label):
        return False, {"deferred": True}, "deferred: free memory stayed below %.1f GB for 30 min" % MEM_NEED_GB
    cmd = ["bash", "-c", "source %s && conda activate ccx223 && exec %s -u scripts/run_ccx_watchdog.py %s %d"
           % (CONDA_SH, SYS_PY, job, RSS_LIMIT_GB)]
    r = run(cmd, label, c, timeout)
    txt = r["text"]
    warns = len(re.findall(r"\*WARNING", txt, re.I))
    if os.path.exists(job + ".ccx.log"):
        warns = max(warns, len(re.findall(r"\*WARNING", open(job + ".ccx.log", errors="replace").read(), re.I)))
    peak = next((ln.strip() for ln in txt.splitlines() if "peak" in ln.lower()), "")
    info = {"wall_s": round(r["wall"], 1), "warnings": warns, "peak_line": peak}
    ok = (r["rc"] == 0 and not r["timed_out"] and "SOLVE OK" in txt and os.path.exists(job + ".dat") and warns == 0)
    if not ok:
        return False, info, "rc=%s timed_out=%s solve_ok=%s dat=%s warnings=%d; see %s" % (
            r["rc"], r["timed_out"], "SOLVE OK" in txt, os.path.exists(job + ".dat"), warns, r["log"])
    return True, info, ""


def step_static_solve(c):
    if not os.path.exists(c.N["job_inp"]):
        return False, {}, "missing %s (pass 1 incomplete)" % c.N["job_inp"]
    return solve(c, c.N["job"], "solve-static", TIMEOUT["solve"])


def step_static_extract(c):
    job = c.N["job"]
    r = run([SYS_PY, "-u", "scripts/extract_stage3_skin.py", job], "extract-static", c, TIMEOUT["extract"])
    out = "results/%s_results.json" % os.path.basename(job)
    if r["rc"] != 0 or r["timed_out"] or not os.path.exists(out):
        return False, {"wall_s": round(r["wall"], 1)}, "rc=%s timed_out=%s json=%s; see %s" % (
            r["rc"], r["timed_out"], os.path.exists(out), r["log"])
    d = json.load(open(out))
    shutil.move(out, "%s/%s_static.json" % (STATIC_DIR, c.key))
    mass = (d["W_solid_N"] + d["W_skin_N"]) / 9.81 / d["n"]
    vals = [d[k] for k in ("uy_spar_axis_mm", "vm_skin_max_MPa", "vm_solid_only_max_MPa")]
    if any((v is None) or math.isnan(v) or math.isinf(v) for v in vals):
        return False, {}, "non-finite value in the extractor results"
    sr = run([SYS_PY, "scripts/stress_metrics_cli.py", c.N["job"]], "stress-metrics", c, TIMEOUT["extract"])
    sm = None
    for line in (sr["text"] or "").splitlines():
        if line.startswith("STRESS_JSON "):
            sm = json.loads(line[len("STRESS_JSON "):])
    if sr["rc"] != 0 or sr["timed_out"] or sm is None:
        return False, {"wall_s": round(r["wall"], 1)}, "stress metrics failed rc=%s; see %s" % (sr["rc"], sr["log"])
    sm = {k: sm[k] for k in ("skin_p99", "skin_p99_9", "skin_max_dTE_ge3", "solid_only_max")}
    c.rec["metrics"].update(tip_mm=d["uy_spar_axis_mm"], mass_kg=mass, solid_mass_kg=d["W_solid_N"] / 9.81 / d["n"],
                            skin_mass_kg=d["W_skin_N"] / 9.81 / d["n"], vm_skin_max_MPa=d["vm_skin_max_MPa"],
                            vm_solid_only_max_MPa=d["vm_solid_only_max_MPa"],
                            vm_skin_max_z_gt_10_MPa=d["vm_skin_max_z_gt_10_MPa"],
                            bracket_lo_mm=d["bracket_lo_mm"], bracket_hi_mm=d["bracket_hi_mm"])
    c.rec["metrics"].update(sm)
    LOG.say("[%s] static: tip %.6f mm, mass %.6f kg, skin max VM %.1f MPa (TE-vertex peaks included)" % (
        c.key, d["uy_spar_axis_mm"], mass, d["vm_skin_max_MPa"]))
    return True, {"wall_s": round(r["wall"], 1)}, ""


def step_buckle_solve(c):
    bjob = c.N["job"] + "_buckle"
    r = run([SYS_PY, "scripts/make_buckle_inp.py", c.N["job_inp"], bjob + ".inp", "6"], "buckle-make", c, TIMEOUT["small"])
    if r["rc"] != 0 or not os.path.exists(bjob + ".inp"):
        return False, {}, "make_buckle_inp failed rc=%s; see %s" % (r["rc"], r["log"])
    return solve(c, bjob, "solve-buckle", TIMEOUT["buckle"])


def step_buckle_extract(c):
    bjob = c.N["job"] + "_buckle"
    r = run([SYS_PY, "scripts/extract_stage3_buckle.py", bjob], "extract-buckle", c, TIMEOUT["small"])
    fac = [(int(a), float(b)) for a, b in re.findall(r"mode (\d+): buckling factor = ([0-9.eE+-]+)", r["text"])]
    crit = grab(r"Critical buckling factor: ([0-9.eE+-]+)", r["text"])
    if r["rc"] != 0 or len(fac) != 6 or crit is None:
        return False, {}, "parsed %d factors, critical=%s rc=%s; see %s" % (len(fac), crit, r["rc"], r["log"])
    vals = [v for _, v in sorted(fac)]
    if min(vals) <= 0 or abs(crit - min(vals)) > 1e-9:
        return False, {}, "buckling factors inconsistent: %s critical %s" % (vals, crit)
    c.rec["metrics"].update(buckling_crit=crit, buckling_factors=vals)
    LOG.say("[%s] buckling: critical %.4f ; modes %s" % (c.key, crit, ", ".join("%.4f" % v for v in vals)))
    return True, {}, ""


def step_cleanup(c):
    n = 0
    for job in (c.N["job"], c.N["job"] + "_buckle"):
        for ext in (".inp", ".frd", ".dat", ".12d", ".cvg"):
            p = job + ext
            if os.path.exists(p):
                os.remove(p)
                n += 1
    LOG.say("[%s] cleanup: removed %d large solver files (.sta and .ccx.log kept)" % (c.key, n))
    return True, {"removed": n}, ""


PASS1 = [("geom", step_geom), ("imprint", step_imprint), ("mesh", step_mesh), ("inp", step_inp), ("load", step_load)]
PASS2 = [("static_solve", step_static_solve), ("static_extract", step_static_extract),
         ("buckle_solve", step_buckle_solve), ("buckle_extract", step_buckle_extract), ("cleanup", step_cleanup)]
SKIN_STEPS = PASS1[3:] + PASS2      # inp, load, then all of pass 2 (mesh/imprint/geometry are reused)


def step_done(c, name):
    s = c.rec["steps"].get(name)
    if not s or s.get("status") != "ok":
        return False
    N = c.N
    if name == "geom":
        return os.path.exists(N["step"])
    if name == "imprint":
        return os.path.exists(N["json"])
    if name == "mesh":
        return os.path.exists(N["raw"])
    if name == "inp":
        return os.path.exists(N["static_inp"]) or os.path.exists(N["job_inp"])
    if name == "load":
        return os.path.exists(N["job_inp"])
    return True


def baseline_check(c):
    m = c.rec["metrics"]
    d_tip = abs(m["tip_mm"] - BASE_TIP) / BASE_TIP
    d_buc = abs(m["buckling_crit"] - BASE_BUCKLE)
    d_mas = abs(m["mass_kg"] - BASE_MASS) / BASE_MASS
    ok = d_tip < 1e-6 and d_buc < 1e-4 and d_mas < 1e-6
    LOG.say("[%s] BASELINE CHECK vs Stage 4: tip rel diff %.2e, buckling abs diff %.2e, mass rel diff %.2e -> %s" % (
        c.key, d_tip, d_buc, d_mas, "PASS" if ok else "FAIL"))
    c.rec["baseline_check"] = dict(passed=ok, tip_rel=d_tip, buckling_abs=d_buc, mass_rel=d_mas)
    return ok


def process(c, steps, retry_failed):
    """Run the steps of one layout. Returns 'ok', 'failed', 'deferred' or 'skipped'."""
    if c.rec.get("failed_step") and not retry_failed:
        LOG.say("[%s] skipped: previously failed at %s (%s); use --retry-failed" % (
            c.key, c.rec["failed_step"], c.rec.get("failed_reason")))
        return "skipped"
    for name, fn in steps:
        if step_done(c, name):
            LOG.say("[%s] %s: already done, skipping" % (c.key, name))
            continue
        LOG.say("[%s] ===== step %s =====" % (c.key, name))
        t0 = time.time()
        ok, info, reason = fn(c)
        info = dict(info)
        if info.get("deferred"):
            LOG.say("[%s] %s: DEFERRED (%s)" % (c.key, name, reason))
            return "deferred"
        info.update(status="ok" if ok else "failed", reason=reason, step_wall_s=round(time.time() - t0, 1))
        c.rec["steps"][name] = info
        if ok:
            c.rec["failed_step"], c.rec["failed_reason"] = None, None
            save_rec(c)
        else:
            c.rec["failed_step"], c.rec["failed_reason"] = name, reason
            save_rec(c)
            LOG.say("[%s] *** FAILED at %s: %s" % (c.key, name, reason))
            return "failed"
    if steps is PASS2 or steps is SKIN_STEPS:
        c.rec["complete"] = True
        save_rec(c)
    return "ok"


# ------------------------------------------------------------------ summary
def write_summary(rows_layouts, rib_t):
    fields = ["id", "key", "family", "stations", "rib_t", "status", "failed_step", "est_rib_mass_g", "rib_mass_g", "tets",
              "s6", "eq_estimate", "min_free_area", "min_shared_area", "mesh_attempts", "tip_mm", "mass_kg",
              "buckling_crit", "b1", "b2", "b3", "b4", "b5", "b6", "vm_skin_max_MPa", "vm_solid_only_max_MPa",
              "skin_p99", "skin_p99_9", "skin_max_dTE_ge3", "solid_only_max"]
    out = "results/stage6_results.csv" if rib_t == DEFAULT_RIB_T else "results/stage6_results_rt%s.csv" % ("%g" % rib_t).replace(".", "p")
    if SKIN_T != 1.0:
        out = out.replace(".csv", "_sk%s.csv" % skin_p())
    counts = {}
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for L in rows_layouts:
            try:
                tag = compute_tag(L, rib_t)
            except SystemExit as e:
                w.writerow(dict(id=L["id"], family=L["family"], stations=L["stations"], rib_t=rib_t, status="rejected: %s" % e))
                continue
            key = tag.lstrip("_") + SKIN_SUFFIX
            row = dict(id=L["id"], key=key, family=L["family"], stations=L["stations"], rib_t=rib_t,
                       est_rib_mass_g=L["est_rib_mass_g"], status="pending")
            p = rec_path(key)
            if os.path.exists(p):
                d = json.load(open(p))
                st = d["steps"]
                m = d.get("metrics", {})
                if d.get("complete"):
                    row["status"] = "complete"
                elif d.get("failed_step"):
                    row["status"] = "failed"
                    row["failed_step"] = d["failed_step"]
                elif st.get("load", {}).get("status") == "ok":
                    row["status"] = "pass1_done"
                row["tets"] = st.get("inp", {}).get("tets")
                row["s6"] = st.get("inp", {}).get("s6")
                row["eq_estimate"] = st.get("inp", {}).get("eq_estimate")
                row["min_free_area"] = st.get("imprint", {}).get("min_free_area")
                row["min_shared_area"] = st.get("imprint", {}).get("min_shared_area")
                row["mesh_attempts"] = len(st.get("mesh", {}).get("attempts", [])) or None
                vol = st.get("imprint", {}).get("volume_expected") or st.get("imprint", {}).get("volume")
                if vol:
                    row["rib_mass_g"] = round((vol - SPAR_V2_VOLUME) * RHO_G_MM3, 3)
                for k in ("tip_mm", "mass_kg", "vm_skin_max_MPa", "vm_solid_only_max_MPa", "skin_p99", "skin_p99_9",
                          "skin_max_dTE_ge3", "solid_only_max"):
                    row[k] = m.get(k)
                row["buckling_crit"] = m.get("buckling_crit")
                for i, v in enumerate(m.get("buckling_factors", []), 1):
                    row["b%d" % i] = v
            counts[row["status"]] = counts.get(row["status"], 0) + 1
            w.writerow(row)
    LOG.say("summary written to %s : %s" % (out, ", ".join("%s=%d" % kv for kv in sorted(counts.items()))))


# ------------------------------------------------------------------ check and main
def do_check(layouts, rib_t):
    ok = True

    def chk(label, cond, extra=""):
        nonlocal ok
        LOG.say("check  %-46s %s %s" % (label, "OK" if cond else "MISSING/FAIL", extra))
        ok = ok and cond

    for p in (SYS_PY, CONDA_ENV + "/bin/ccx", CONDA_SH):
        chk(p, os.path.exists(p))
    for s in ("build_stage3_solid", "build_stage3_skin_geom", "mesh_stage3_skin_raw", "build_stage3_skin_inp",
              "apply_stage4_pressure_load", "make_buckle_inp", "run_ccx_watchdog", "extract_stage3_skin",
              "extract_stage3_buckle", "variant_tag"):
        chk("scripts/%s.py" % s, os.path.exists("scripts/%s.py" % s))
    chk("layout table %s" % LAYOUT_CSV, os.path.exists(LAYOUT_CSV), "(%d layouts)" % len(layouts))
    free = free_disk_gb()
    chk("real free disk space >= 10 GB", free >= 10.0, "(%.1f GB)" % free)
    LOG.say("check  free memory now: %.1f GB (a solve needs %.1f GB)" % (mem_avail_gb(), MEM_NEED_GB))
    bad = 0
    for L in layouts:
        try:
            compute_tag(L, rib_t)
        except SystemExit as e:
            bad += 1
            LOG.say("check  layout %s REJECTED by variant_tag: %s" % (L["id"], e))
    chk("all layouts accepted by variant_tag", bad == 0, "(%d rejected)" % bad)
    r = subprocess.run(["bash", "-c", "source %s && conda activate ccx223 && which ccx && %s -c 'import scipy, numpy, pyvista; print(\"scipy\", scipy.__version__)'"
                        % (CONDA_SH, SYS_PY)], capture_output=True, text=True)
    LOG.say("check  conda activation: rc=%d | %s" % (r.returncode, " | ".join((r.stdout + r.stderr).split("\n")[:4])))
    ok = ok and r.returncode == 0
    LOG.say("CHECK %s" % ("PASSED" if ok else "FAILED"))
    return ok


def main():
    global LOG, SKIN_T, SKIN_SUFFIX
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["check", "pass1", "pass2", "summary", "skin"])
    ap.add_argument("--layouts", default="all")
    ap.add_argument("--rib-t", type=float, default=DEFAULT_RIB_T)
    ap.add_argument("--skin-t", type=float, default=1.0)
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--max-consecutive-fail", type=int, default=3)
    a = ap.parse_args()
    if a.stage == "skin" and a.skin_t == 1.0:
        sys.exit("ERROR: --stage skin needs --skin-t other than 1.0 (1.0 mm is the existing pass 1/pass 2 study).")
    if a.stage in ("pass1", "pass2") and a.skin_t != 1.0:
        sys.exit("ERROR: --skin-t other than 1.0 is only for --stage skin (or summary).")
    if a.skin_t <= 0:
        sys.exit("ERROR: --skin-t must be positive.")
    SKIN_T = a.skin_t
    SKIN_SUFFIX = "" if SKIN_T == 1.0 else "_sk" + skin_p()
    os.makedirs(RES, exist_ok=True)
    lock = open(RES + "/driver.lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit("ERROR: another driver instance holds %s/driver.lock -- STOP." % RES)
    LOG = Log()
    allrows = list(csv.DictReader(open(LAYOUT_CSV)))
    if a.layouts == "all":
        sel = allrows
    else:
        want = [s.strip() for s in a.layouts.split(",")]
        byid = {r["id"]: r for r in allrows}
        unknown = [w for w in want if w not in byid]
        if unknown:
            sys.exit("ERROR: unknown layout ids %s" % unknown)
        sel = [byid[w] for w in want]
    LOG.say("=== driver start: stage=%s layouts=%d rib_t=%g skin_t=%g retry_failed=%s ===" % (a.stage, len(sel), a.rib_t, SKIN_T, a.retry_failed))
    if a.stage == "check":
        sys.exit(0 if do_check(sel, a.rib_t) else 1)
    if a.stage == "summary":
        write_summary(allrows, a.rib_t)
        return
    free = free_disk_gb()
    if free < (10.0 if a.stage == "pass1" else 8.0):
        sys.exit("ERROR: only %.1f GB free disk -- STOP." % free)
    steps = PASS1 if a.stage == "pass1" else (SKIN_STEPS if a.stage == "skin" else PASS2)
    t_start = time.time()
    done, consecutive, results = 0, 0, {}
    try:
        for i, L in enumerate(sel, 1):
            try:
                c = Ctx(L, a.rib_t)
            except SystemExit as e:
                LOG.say("[%s] REJECTED by variant_tag: %s" % (L["id"], e))
                continue
            c.rec = load_rec(c)
            if a.stage == "skin":
                miss = [k for k in ("raw", "json", "vol") if not os.path.exists(c.N[k])]
                if miss:
                    LOG.say("[%s] skin stage skipped: missing reusable files %s" % (c.key, miss))
                    continue
            if a.stage == "pass2" and not step_done(c, "load"):
                LOG.say("[%s] pass 2 skipped: pass 1 not complete for this layout" % c.key)
                continue
            LOG.say("----- layout %s (%d of %d)  stations %s  rib_t %g -----" % (c.key, i, len(sel), L["stations"], a.rib_t))
            res = process(c, steps, a.retry_failed)
            results[res] = results.get(res, 0) + 1
            if res == "ok":
                consecutive = 0
                done += 1
                if a.stage == "pass2" and L["id"] == "L000" and a.rib_t == DEFAULT_RIB_T:
                    if not baseline_check(c):
                        save_rec(c)
                        LOG.say("*** BASELINE MISMATCH: pass 2 aborted. Everything after this point is suspect. ***")
                        write_summary(allrows, a.rib_t)
                        sys.exit(3)
                    save_rec(c)
            elif res == "failed":
                consecutive += 1
                if consecutive >= a.max_consecutive_fail:
                    LOG.say("*** %d consecutive failures: stopping (probably systematic). ***" % consecutive)
                    write_summary(allrows, a.rib_t)
                    sys.exit(2)
            el = time.time() - t_start
            if done:
                LOG.say("progress: %d/%d layouts processed, %d ok; elapsed %.0f min; avg %.1f min per ok layout; ETA %.0f min"
                        % (i, len(sel), done, el / 60, el / 60 / done, (len(sel) - i) * el / 60 / max(done, 1)))
    except KeyboardInterrupt:
        LOG.say("*** interrupted by the user; children killed; records are saved and the run can be resumed ***")
        write_summary(allrows, a.rib_t)
        sys.exit(130)
    LOG.say("=== %s finished: %s ; elapsed %.1f min ===" % (a.stage, results, (time.time() - t_start) / 60))
    write_summary(allrows, a.rib_t)


if __name__ == "__main__":
    main()
