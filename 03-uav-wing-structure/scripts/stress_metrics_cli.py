#!/usr/bin/env python3
"""
Print the Stage-4 skin/solid stress statistics of one solved job as a single JSON line.
Usage:  /usr/bin/python3 scripts/stress_metrics_cli.py <job path without extension>
Output: STRESS_JSON {"skin_p99": ..., "skin_p99_9": ..., "skin_max_dTE_ge3": ..., "solid_only_max": ...}
Thin wrapper around collect_stage3_results.stress_stats so the driver can call it like any other step.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_stage3_results import stress_stats  # noqa: E402

KEYS = ("skin_p99", "skin_p99_9", "skin_max_dTE_ge3", "solid_only_max")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: stress_metrics_cli.py <job>")
    st, _ = stress_stats(sys.argv[1])
    out = {}
    for k in KEYS:
        if k not in st:
            sys.exit("stress_stats did not return key %s (keys: %s)" % (k, sorted(st)))
        v = float(st[k])
        if math.isnan(v) or math.isinf(v):
            sys.exit("non-finite value for %s" % k)
        out[k] = v
    print("STRESS_JSON " + json.dumps(out))


if __name__ == "__main__":
    main()
