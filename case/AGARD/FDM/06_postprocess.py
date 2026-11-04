#!/usr/bin/env python3
"""Stage 06 - damping rate, frequency and per-cycle amplitude from
modal_diagnostics.csv.

The file is written by fastDynamicFvMesh once per completed time step, but the
column layout is read by name, not position:

    Time, Force_1, Force_2, PressureForce_*, ..., Disp_1, Vel_1, Acc_1, ...,
    Disp_2, Vel_2, Acc_2, AppliedDisp_2, StructuralScale_2

Method: one row per time (the last row of each step), modal-1 displacement is
mean-removed, cycles are bounded by upward zero crossings, and

    sigma = -ln(A_{n+1} / A_n) / T_n+1

with A = peak-to-peak amplitude inside a cycle. sigma > 0 means decaying.
"""
import argparse
import csv
import json
import math
import os
import sys


def load(path):
    with open(path) as f:
        header = [h.strip() for h in f.readline().strip().split(",")]
    idx = {n: i for i, n in enumerate(header)}
    if "Disp_1" not in idx or "Time" not in idx:
        sys.exit("[06] %s lacks Time/Disp_1 columns" % path)
    rows, last_t = [], None
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.strip().split(",")
            if len(p) < len(header):
                continue
            t = float(p[0])
            if last_t is not None and t == last_t:
                rows[-1] = p
            else:
                rows.append(p)
                last_t = t
    return idx, rows


def analyse(path):
    idx, rows = load(path)
    if len(rows) < 5:
        return None
    T = [float(r[0]) for r in rows]
    d1 = [float(r[idx["Disp_1"]]) for r in rows]
    d2 = [float(r[idx["Disp_2"]]) for r in rows] if "Disp_2" in idx else [0.0] * len(rows)
    n = len(T)
    tail = d1[int(0.6 * n):] or d1
    off = sum(tail) / len(tail)
    s = [x - off for x in d1]

    cross = []
    for i in range(1, n):
        if s[i - 1] < 0.0 <= s[i]:
            w = -s[i - 1] / (s[i] - s[i - 1])
            cross.append(T[i - 1] + w * (T[i] - T[i - 1]))

    cycles = []
    for k in range(len(cross) - 1):
        t0, t1 = cross[k], cross[k + 1]
        i0 = next(i for i in range(n) if T[i] >= t0)
        i1 = next(i for i in range(n) if T[i] >= t1)
        pp1 = max(d1[i0:i1 + 1]) - min(d1[i0:i1 + 1])
        pp2 = max(d2[i0:i1 + 1]) - min(d2[i0:i1 + 1])
        cycles.append({"t_start": t0, "period": t1 - t0, "f_Hz": 1.0 / (t1 - t0),
                       "pp_1": pp1, "pp_2": pp2})

    for k in range(1, len(cycles)):
        c, p = cycles[k], cycles[k - 1]
        c["sigma"] = -math.log(c["pp_1"] / p["pp_1"]) / c["period"] if p["pp_1"] > 0 else None

    sig = [c.get("sigma") for c in cycles if c.get("sigma") is not None]
    return {"cycles": cycles,
            "sigma_mean": sum(sig) / len(sig) if sig else None,
            "sigma_last": sig[-1] if sig else None,
            "f_mean": sum(c["f_Hz"] for c in cycles) / len(cycles) if cycles else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-root", required=True)
    ap.add_argument("--csv", default=None, help="default: <case>/modal_diagnostics.csv")
    ap.add_argument("--json", default=None, help="write the summary as JSON")
    args = ap.parse_args()

    path = args.csv or os.path.join(args.case_root, "modal_diagnostics.csv")
    if not os.path.isfile(path):
        sys.exit("[06] not found: %s" % path)

    res = analyse(path)
    if res is None or len(res["cycles"]) < 2:
        # Not an error: a short smoke run covers less than one structural
        # period, so no damping rate exists yet. Report what is available.
        idx, rows = load(path)
        d1 = [float(r[idx["Disp_1"]]) for r in rows]
        d2 = [float(r[idx["Disp_2"]]) for r in rows] if "Disp_2" in idx else [0.0] * len(rows)
        print("[06] %s" % path)
        print("[06]   time span      : %.4g s over %d steps" % (float(rows[-1][0]) - float(rows[0][0]), len(rows)))
        print("[06]   Disp_1 range   : %.3e .. %.3e" % (min(d1), max(d1)))
        print("[06]   Disp_2 range   : %.3e .. %.3e" % (min(d2), max(d2)))
        ncyc = len(res["cycles"]) if res else 0
        print("[06]   cycles found   : %d  (need >= 2 for sigma)" % ncyc)
        print("[06]   INSUFFICIENT DATA: run at least ~3 structural periods "
              "(endTime >= 0.3 for this wing) to obtain sigma/f.")
        return 0

    print("[06] %s" % path)
    print("[06]   cyc  t_start   T[s]     f[Hz]    Disp_1 p2p    Disp_2 p2p   sigma[1/s]")
    for k, c in enumerate(res["cycles"], 1):
        sg = "" if c.get("sigma") is None else "%+8.2f" % c["sigma"]
        print("[06]   %2d  %8.4f  %7.4f  %6.2f  %12.4e  %12.4e  %s"
              % (k, c["t_start"], c["period"], c["f_Hz"], c["pp_1"], c["pp_2"], sg))

    if res["sigma_mean"] is not None:
        print("[06]   sigma(mean of cycles 2..) = %+.2f 1/s   f = %.2f Hz"
              % (res["sigma_mean"], res["f_mean"]))
        verdict = "decaying" if res["sigma_mean"] > 0 else "growing"
        print("[06]   verdict: %s" % verdict)

    if args.json:
        with open(args.json, "w") as f:
            json.dump(res, f, indent=2)
        print("[06]   written: %s" % args.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
