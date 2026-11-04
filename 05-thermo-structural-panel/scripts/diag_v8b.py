#!/usr/bin/python3
"""V8b diagnostic: how does the structural run's temperature evolve over the step?

Usage: /usr/bin/python3 scripts/diag_v8b.py
Uses existing results only (no solver run). Hypothesis tested: the structural step ramps linearly from the initial
temperature to the LAST frame of the thermal results, T(t) = T0 + (t / T_END) (T_last - T0).
Prints a table for the centre-column top node (run B) and the largest deviation from this ramp rule over
all run-B times (centre column) and all run-A times (all nodes).
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
import run_v8b as v8b

T0, TEND, DT = v8b.v8a.T0, v8b.v8a.T_END, v8b.DT_T


def main():
    root = os.path.dirname(HERE)
    sim = os.path.join(root, "simulation")
    nodes = ccxio.read_nodes(os.path.join(root, "mesh", v8b.MESH_S + ".inp"))
    th = v8b.temp_blocks(os.path.join(sim, "v8b_thermal", "v8b_thermal.dat"), "COL_CEN")
    sb = v8b.temp_blocks(os.path.join(sim, "v8b_struct_b", "v8b_struct_b.dat"), "COL_CEN")
    sa = v8b.temp_blocks(os.path.join(sim, "v8b_struct_a", "v8b_struct_a.dat"), "NALL")
    ref = v8b.temp_blocks(os.path.join(sim, "v8a2_n20_dt0p5", "v8a2_n20_dt0p5.dat"), "NALL")
    ids, T_last = th[TEND]
    top = int(np.argmax([nodes[int(i)][2] for i in ids]))
    TT = dict(th)
    TT[0.0] = (ids, np.full(len(ids), T0))

    def thermal_at(t):
        if t >= TEND - 1e-9:
            return T_last
        k0 = int(math.floor(t / DT + 1e-9))
        t0, t1 = round(k0 * DT, 6), round((k0 + 1) * DT, 6)
        f = (t - t0) / DT
        return (1 - f) * TT[t0][1] + f * TT[t1][1]

    print("centre-column top node id %d, z = %.5f m" % (ids[top], nodes[int(ids[top])][2]))
    print("    t[s]   struct run B [C]   thermal history, interpolated [C]   ramp-to-last-frame rule [C]")
    for t in (0.25, 1.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 29.75, 30.0):
        s = sb[round(t, 6)][1][top]
        ramp = T0 + (t / TEND) * (T_last - T0)[top]
        print("  %6.2f %18.4f %34.4f %27.4f" % (t, s, thermal_at(t)[top], ramp))
    dev_b = max(float(np.abs(blk[1] - (T0 + (t / TEND) * (T_last - T0))).max()) for t, blk in sb.items())
    ramp_dev_thermal = max(float(np.abs(thermal_at(t) - (T0 + (t / TEND) * (T_last - T0))).max()) for t in sb)
    print("run B: max |struct - ramp rule| over all %d times, centre column: %.3e K" % (len(sb), dev_b))
    print("for scale: max |thermal history - ramp rule| over the same times: %.3e K" % ramp_dev_thermal)
    T30 = ref[TEND][1]
    dev_a = max(float(np.abs(blk[1] - (T0 + (t / TEND) * (T30 - T0))).max()) for t, blk in sa.items())
    print("run A: max |struct - ramp rule| over t = %s s, all %d nodes: %.3e K" % (sorted(sa), len(T30), dev_a))


if __name__ == "__main__":
    main()
