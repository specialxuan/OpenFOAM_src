#!/usr/bin/env python3
"""Stage 03 - converged steady RANS initial field.

The flutter case needs a steady solution as its initial field. This stage
builds a throw-away steady case in <case-root>/.steady (static mesh, no FSI,
no modes), solves it with myHisa, and copies the converged fields into
<case-root>/0/.

Templates: templates/steady  (dynamicFvMesh staticFvMesh, dualTime steadyState,
residual/pseudo-time tolerance in fvSolution).
Reference state: M=0.96, q=2936 Pa (same as the transient templates).

Cost: ~12.5 min on 32 ranks, ~25 min on 8 ranks.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tpl  # noqa: E402  (templates live in templates/*.zip)

Q_REF = 2936.0
P_INF_REF = 4551.7832
T_INF_REF = 257.7964
FIELDS = ("p", "U", "T", "rho", "k", "omega", "nut", "alphat")


def check(cond, msg):
    if not cond:
        sys.exit("[03] " + msg)


def read(p):
    with open(p) as f:
        return f.read()


def write(p, s):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        f.write(s)


def sub_once(text, pattern, repl, what):
    new, n = re.subn(pattern, repl, text, count=1)
    if n != 1:
        sys.exit("[03] could not patch %s" % what)
    return new


def build_steady_case(case_root, pressure, nprocs):
    st = os.path.join(case_root, ".steady")
    shutil.rmtree(st, ignore_errors=True)
    os.makedirs(os.path.join(st, "constant"))
    os.makedirs(os.path.join(st, "system"))

    # static mesh: reuse the case mesh (symlink keeps it cheap)
    shutil.copytree(os.path.join(case_root, "constant", "polyMesh"),
                    os.path.join(st, "constant", "polyMesh"))

    s = pressure / Q_REF
    p_inf = P_INF_REF * s

    for f in ("g", "thermophysicalProperties", "turbulenceProperties"):
        tpl.copy_out(tpl.STEADY_ZIP, "steady", "constant/" + f,
                     os.path.join(st, "constant", f))

    dmd = tpl.text(tpl.STEADY_ZIP, "steady", "constant/dynamicMeshDict")
    dmd = sub_once(dmd, r"pRef\s+\S+;", "pRef               %.10g;" % p_inf, "pRef") \
        if re.search(r"pRef\s+\S+;", dmd) else dmd
    write(os.path.join(st, "constant", "dynamicMeshDict"), dmd)

    for f in ("controlDict", "fvSchemes", "fvSolution", "decomposeParDict"):
        t = tpl.text(tpl.STEADY_ZIP, "steady", "system/" + f)
        if f == "decomposeParDict":
            t = sub_once(t, r"numberOfSubdomains\s+\d+;",
                         "numberOfSubdomains %d;" % nprocs, "numberOfSubdomains")
        write(os.path.join(st, "system", f), t)

    # initial field: uniform free-stream, p/rho scaled with q
    for name in FIELDS:
        t = tpl.text(tpl.MYHISA_ZIP, "myHisa", "0/" + name)
        if name in ("p", "rho"):
            val = p_inf if name == "p" else p_inf / (287.1 * T_INF_REF)
            t = re.sub(r"uniform\s+[-+0-9.eE]+(?:[ \t]+[-+0-9.eE]+)*",
                       "uniform %.8g" % val, t)
        if name == "U":
            t = re.sub(r"internalField\s+uniform\s+\([^)]*\)\s*;",
                       "internalField   uniform (309.023 0 0);", t)
        write(os.path.join(st, "0", name), t)

    # the static mesh dict must not reference FSI patches lacking fields
    return st


def run(cmd, cwd, timeout):
    print("[03] + %s" % " ".join(cmd))
    r = subprocess.run(cmd, cwd=cwd, timeout=timeout,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-4000:])
        sys.exit("[03] command failed (%d): %s" % (r.returncode, " ".join(cmd)))
    return r.stdout


def latest_time(st):
    times = []
    for d in os.listdir(st):
        try:
            v = float(d)
        except ValueError:
            continue
        if v > 0:
            times.append((v, d))
    return max(times)[1] if times else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-root", required=True)
    ap.add_argument("--mach", type=float, default=0.96)
    ap.add_argument("--pressure", type=float, default=2936.0)
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=7200, help="seconds")
    ap.add_argument("--keep", action="store_true", help="keep .steady/ afterwards")
    args = ap.parse_args()

    cr = args.case_root
    check(os.path.isdir(os.path.join(cr, "constant", "polyMesh")),
          "mesh missing - run stage 01 first")
    check(os.path.isfile(tpl.STEADY_ZIP), "templates/steady.zip not found")
    if args.mach != 0.96:
        print("[03] WARNING: steady templates were exported at M=0.96")

    st = build_steady_case(cr, args.pressure, args.nprocs)
    print("[03] steady case: %s" % st)

    run(["decomposePar", "-case", st, "-force"], cr, args.timeout)
    run(["mpirun", "--allow-run-as-root", "-np", str(args.nprocs),
         "myHisa", "-parallel", "-case", st], cr, args.timeout)

    # a parallel run writes processor*/<time> only; bring the result back to the
    # case root before picking the latest time
    run(["reconstructPar", "-case", st, "-latestTime"], cr, args.timeout)

    t = latest_time(st)
    if t is None:
        # fall back: newest time present in processor0
        p0 = os.path.join(st, "processor0")
        cand = []
        if os.path.isdir(p0):
            for d in os.listdir(p0):
                try:
                    v = float(d)
                except ValueError:
                    continue
                if v > 0:
                    cand.append((v, d))
        check(cand, "steady run produced no time directory")
        t = max(cand)[1]
        print("[03] NOTE: root has no time dir; using processor0/%s" % t)
    src = os.path.join(st, t)
    print("[03] converged time: %s" % t)

    # stage 04 renders 0/ later; create it here so the steady field can land
    zero = os.path.join(cr, "0")
    os.makedirs(zero, exist_ok=True)
    for name in FIELDS:
        f = os.path.join(src, name)
        if not os.path.isfile(f):
            sys.exit("[03] steady result lacks field %s" % name)
        shutil.copy2(f, os.path.join(zero, name))
    print("[03] copied %d fields into %s" % (len(FIELDS), zero))

    if not args.keep:
        shutil.rmtree(st, ignore_errors=True)
        print("[03] removed .steady/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
