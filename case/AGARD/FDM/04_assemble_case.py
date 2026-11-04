#!/usr/bin/env python3
"""Stage 04 - render the myHisa case dictionaries from the templates.

Input : templates/myHisa/{system,constant,0}
        <case-root>/constant/polyMesh  (stage 01)
        <case-root>/mode              (stage 02)
Output: <case-root>/system/{controlDict,fvSchemes,fvSolution,decomposeParDict}
        <case-root>/constant/{dynamicMeshDict,g,thermophysicalProperties,
                              turbulenceProperties}
        <case-root>/0/{p,rho,T,U,k,omega,nut,alphat}

Reference state baked into the templates: M = 0.96, q = 2936 Pa
(p_inf = 4551.7832 Pa, T_inf = 257.7964 K, U_inf = 309 m/s). At fixed Mach the
free-stream pressure and density scale linearly with q, so only p and rho are
rescaled; T, U and the turbulence quantities stay as they are.

The initial field is a uniform free-stream state. Stage 03 replaces it with a
converged steady solution.
"""
import argparse
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tpl  # noqa: E402  (templates live in templates/*.zip)

# reference free-stream state of the templates (M=0.96, q=2936 Pa)
Q_REF = 2936.0
P_INF_REF = 4551.7832
T_INF_REF = 257.7964
U_INF_REF = "(309 0 0)"
AVAILABLE_MODES = 4          # data/modes_4m.zip


def read(path):
    with open(path) as f:
        return f.read()


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def sub_once(text, pattern, repl, what):
    new, n = re.subn(pattern, repl, text, count=1)
    if n != 1:
        sys.exit("[04] could not patch %s (pattern: %s)" % (what, pattern))
    return new


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-root", required=True)
    ap.add_argument("--mach", type=float, default=0.96)
    ap.add_argument("--pressure", type=float, default=2936.0, help="free-stream q [Pa]")
    ap.add_argument("--modes", type=int, default=2, help="modes to superpose (nActiveModes)")
    ap.add_argument("--amr", dest="amr", action="store_true", default=True)
    ap.add_argument("--no-amr", dest="amr", action="store_false")
    ap.add_argument("--end-time", type=float, default=0.02, help="default: 20 steps (smoke)")
    ap.add_argument("--delta-t", type=float, default=1e-3)
    ap.add_argument("--write-interval", type=int, default=50)
    ap.add_argument("--nprocs", type=int, default=8)
    args = ap.parse_args()

    cr = args.case_root
    if not os.path.isdir(os.path.join(cr, "constant", "polyMesh")):
        sys.exit("[04] mesh missing - run stage 01 first")
    if not os.path.isdir(os.path.join(cr, "mode")):
        sys.exit("[04] mode/ missing - run stage 02 first")
    if not os.path.isfile(tpl.MYHISA_ZIP):
        sys.exit("[04] template archive not found: %s" % tpl.MYHISA_ZIP)

    if args.mach != 0.96:
        print("[04] WARNING: templates were exported at M=0.96; M=%.3f is applied to "
              "the dictionaries but the free-stream state (T, U) is NOT recomputed."
              % args.mach)

    s = args.pressure / Q_REF
    p_inf = P_INF_REF * s
    u_inf = args.mach * (1.4 * 287.1 * T_INF_REF) ** 0.5

    # ---- system ------------------------------------------------------------
    ctrl = tpl.text(tpl.MYHISA_ZIP, "myHisa", "system/controlDict")
    ctrl = sub_once(ctrl, r"endTime\s+\S+;", "endTime         %g;" % args.end_time, "endTime")
    ctrl = sub_once(ctrl, r"deltaT\s+\S+;", "deltaT          %g;" % args.delta_t, "deltaT")
    ctrl = sub_once(ctrl, r"writeInterval\s+\d+;",
                    "writeInterval   %d;" % args.write_interval, "writeInterval")
    ctrl = sub_once(ctrl, r"startFrom\s+\S+;", "startFrom       startTime;", "startFrom")
    ctrl = sub_once(ctrl, r"startTime\s+\S+;", "startTime       0;", "startTime")
    write(os.path.join(cr, "system", "controlDict"), ctrl)

    for f in ("fvSchemes", "fvSolution"):
        write(os.path.join(cr, "system", f), tpl.text(tpl.MYHISA_ZIP, "myHisa", "system/" + f))

    dec = tpl.text(tpl.MYHISA_ZIP, "myHisa", "system/decomposeParDict")
    dec = sub_once(dec, r"numberOfSubdomains\s+\d+;",
                   "numberOfSubdomains %d;" % args.nprocs, "numberOfSubdomains")
    write(os.path.join(cr, "system", "decomposeParDict"), dec)

    # ---- constant ----------------------------------------------------------
    for f in ("g", "thermophysicalProperties", "turbulenceProperties"):
        tpl.copy_out(tpl.MYHISA_ZIP, "myHisa", "constant/" + f, os.path.join(cr, "constant", f))

    dmd = tpl.text(tpl.MYHISA_ZIP, "myHisa", "constant/dynamicMeshDict")
    dmd = sub_once(dmd, r"pRef\s+\S+;", "pRef               %.10g;" % p_inf, "pRef")
    # modalMass/modalDamp must cover ALL modes present in mode/, not just the
    # ones superposed: nActiveModes truncates the superposition, it does not
    # shorten these lists.
    dmd = sub_once(dmd, r"modalMass\s+\([^)]*\);",
                   "modalMass          (%s);" % " ".join(["1"] * AVAILABLE_MODES),
                   "modalMass")
    dmd = sub_once(dmd, r"modalDamp\s+\([^)]*\);",
                   "modalDamp          (%s);" % " ".join(["0"] * AVAILABLE_MODES),
                   "modalDamp")

    if args.modes < AVAILABLE_MODES:
        if re.search(r"nActiveModes", dmd):
            dmd = re.sub(r"nActiveModes\s+\d+;", "nActiveModes       %d;" % args.modes, dmd)
        else:
            dmd = sub_once(dmd, r"(modalDamp\s+\([^)]*\);)",
                           r"\1\n    nActiveModes       %d;" % args.modes, "insert nActiveModes")
        print("[04] nActiveModes = %d (mode/ holds %d)" % (args.modes, AVAILABLE_MODES))

    if not args.amr:
        dmd = sub_once(dmd, r"refineInterval\s+\d+;",
                       "refineInterval     100000000;", "refineInterval (AMR off)")
        dmd = sub_once(dmd, r"meshRefinementSupport\s+\S+;",
                       "meshRefinementSupport        false;", "meshRefinementSupport")
        print("[04] AMR disabled")
    else:
        print("[04] AMR enabled (refineInterval 5)")

    write(os.path.join(cr, "constant", "dynamicMeshDict"), dmd)

    # ---- 0/ ----------------------------------------------------------------
    zero = os.path.join(cr, "0")
    for name in tpl.names(tpl.MYHISA_ZIP, "myHisa", "0/"):
        t = tpl.text(tpl.MYHISA_ZIP, "myHisa", "0/" + name)

        if name in ("p", "rho"):
            val = p_inf if name == "p" else p_inf / (287.1 * T_INF_REF)
            t = re.sub(r"uniform\s+[-+0-9.eE]+(?:[ \t]+[-+0-9.eE]+)*",
                       "uniform %.8g" % val, t)
        elif name == "T":
            t = sub_once(t, r"internalField\s+uniform\s+\S+\s*;",
                         "internalField   uniform %g;" % T_INF_REF, "T internalField")
        elif name == "U":
            t = sub_once(t, r"internalField\s+uniform\s+\([^)]*\)\s*;",
                         "internalField   uniform (%g 0 0);" % u_inf, "U internalField")

        write(os.path.join(zero, name), t)

    print("[04] case rendered: %s" % cr)
    print("[04]   q = %g Pa (scale %.4f), p_inf = %.6g Pa, modes = %d, endTime = %g"
          % (args.pressure, s, p_inf, args.modes, args.end_time))
    return 0


if __name__ == "__main__":
    sys.exit(main())
