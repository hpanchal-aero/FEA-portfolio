"""
Run CalculiX 2.23 on a job with an RSS-polling memory watchdog.

Usage:
    /usr/bin/python3 scripts/run_ccx_watchdog.py <job_path_without_ext> [rss_limit_gb]

Example:
    /usr/bin/python3 scripts/run_ccx_watchdog.py simulation/config_a_h3p0 8

The watchdog polls the solver's actual resident memory (VmRSS from
/proc/<pid>/status) every 0.2 s and sends SIGKILL if it exceeds the limit.
It deliberately does NOT use RLIMIT_AS: that caps virtual address space, not
physical memory, and caused spurious allocation failures in Project 01.
"""

import os
import subprocess
import sys
import threading
import time

CCX_ENV = "/home/harsh/miniconda3/envs/ccx223"
CCX = os.path.join(CCX_ENV, "bin", "ccx")
POLL_S = 0.2


def rss_mb(pid):
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except (FileNotFoundError, ProcessLookupError):
        return None
    return None


def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: run_ccx_watchdog.py <job_path_without_ext> [rss_limit_gb]")
        sys.exit(1)
    job = sys.argv[1]
    limit_gb = float(sys.argv[2]) if len(sys.argv) == 3 else 8.0
    limit_mb = limit_gb * 1024.0

    inp = job + ".inp"
    log = job + ".ccx.log"
    if not os.path.isfile(inp):
        print(f"ERROR: {inp} not found.")
        sys.exit(1)
    if not os.path.isfile(CCX):
        print(f"ERROR: solver binary not found at {CCX}")
        sys.exit(1)

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.path.join(CCX_ENV, "lib") + ":" + env.get("LD_LIBRARY_PATH", "")

    print(f"Job: {job}   RSS limit: {limit_gb:.1f} GB   log: {log}")
    t0 = time.time()
    with open(log, "w") as lf:
        proc = subprocess.Popen([CCX, "-i", job], stdout=lf, stderr=subprocess.STDOUT, env=env)

    state = {"peak": 0.0, "killed": False}

    def watch():
        while proc.poll() is None:
            r = rss_mb(proc.pid)
            if r is not None:
                if r > state["peak"]:
                    state["peak"] = r
                if r > limit_mb:
                    state["killed"] = True
                    os.kill(proc.pid, 9)
                    return
            time.sleep(POLL_S)

    th = threading.Thread(target=watch, daemon=True)
    th.start()
    rc = proc.wait()
    th.join()
    wall = time.time() - t0

    print(f"\nExit code: {rc}   wall time: {wall:.1f} s   peak RSS: {state['peak']:.0f} MB "
          f"({state['peak'] / 1024:.2f} GB)")
    if state["killed"]:
        print(f"WATCHDOG KILL: RSS exceeded {limit_gb:.1f} GB. Solve aborted.")
        sys.exit(2)

    with open(log) as lf:
        lines = lf.read().splitlines()
    print("\n--- last 15 lines of solver log ---")
    for ln in lines[-15:]:
        print(ln)
    print("-----------------------------------")

    finished = any("Job finished" in ln for ln in lines)
    frd = job + ".frd"
    frd_ok = os.path.isfile(frd) and os.path.getsize(frd) > 0
    print(f"'Job finished' in log: {finished}")
    print(f".frd present and non-empty: {frd_ok}")
    if rc != 0 or not finished or not frd_ok:
        print("SOLVE FAILED -- do not proceed; send me the log.")
        sys.exit(3)
    print("SOLVE OK")


if __name__ == "__main__":
    main()
