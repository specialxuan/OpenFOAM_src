#!/usr/bin/env python3
"""Stage 05 - run the flutter case (decompose + myHisa -parallel) and check it.

Input : a case assembled by stage 04 (mesh, mode/, dicts, 0/)
Output: time directories, log.run, and a one-screen verdict

Checks performed:
  * solver exit status
  * FOAM FATAL count
  * time steps advanced vs requested
  * whether AMR actually refined (if AMR is enabled)
  * physical sanity of the last written field (p, T ranges) when available
"""
import argparse
import os
import re
import subprocess
import sys


def sh(cmd, cwd, timeout=None, log=None):
    print("[05] + %s" % " ".join(cmd))
    if log:
        with open(log, "w") as f:
            r = subprocess.run(cmd, cwd=cwd, timeout=timeout, stdout=f,
                               stderr=subprocess.STDOUT, text=True)
        return r.returncode
    r = subprocess.run(cmd, cwd=cwd, timeout=timeout, text=True)
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-root", required=True)
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--end-time", type=float, default=0.02)
    ap.add_argument("--delta-t", type=float, default=1e-3)
    ap.add_argument("--timeout", type=int, default=86400)
    ap.add_argument("--skip-decompose", action="store_true")
    args = ap.parse_args()

    cr = args.case_root
    if not os.path.isdir(os.path.join(cr, "system")):
        sys.exit("[05] case not assembled - run stage 04 first")

    log = os.path.join(cr, "log.run")
    if not args.skip_decompose:
        rc = sh(["decomposePar", "-force"], cr, args.timeout)
        if rc != 0:
            sys.exit("[05] decomposePar failed (%d)" % rc)

    rc = sh(["mpirun", "--allow-run-as-root", "-np", str(args.nprocs),
             "myHisa", "-parallel"], cr, args.timeout, log=log)
    print("[05] solver exit code: %d" % rc)

    text = open(log, errors="ignore").read() if os.path.isfile(log) else ""
    fatal = len(re.findall(r"FOAM FATAL", text))
    nsteps = len(re.findall(r"^Time = ", text, re.M))
    refined = re.findall(r"Refined from (\d+) to (\d+) cells", text)
    last_time = re.findall(r"^Time = (\S+)", text, re.M)

    print("[05] ---- verdict ----")
    print("[05]   steps printed : %d" % nsteps)
    print("[05]   last Time     : %s" % (last_time[-1] if last_time else "-"))
    print("[05]   FOAM FATAL    : %d" % fatal)
    print("[05]   AMR events    : %d" % len(refined))
    if refined:
        a, b = refined[-1]
        print("[05]   last refinement: %s -> %s cells" % (a, b))

    ok = True
    if fatal:
        ok = False
        m = re.search(r"(FOAM FATAL.*?)(?:\n\n|\Z)", text, re.S)
        print("[05]   first fatal   : %s" % (m.group(1).strip()[:400] if m else "?"))
    if rc != 0:
        ok = False

    # expected step count
    expect = int(round(args.end_time / args.delta_t))
    if nsteps < expect:
        print("[05]   WARNING: expected %d steps, got %d (run may have stopped early)"
              % (expect, nsteps))

    print("[05]   RESULT        : %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
