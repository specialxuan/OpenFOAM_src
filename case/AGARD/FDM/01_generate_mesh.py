#!/usr/bin/env python3
"""Stage 01 - unpack the Fluid8 mesh into a case directory.

Input : <repo>/case/AGARD/FDM/data/fluid8_mesh.zip  (contains shared/polyMesh/*)
Output: <case-root>/constant/polyMesh/{points,faces,owner,neighbour,boundary,sets}

The archive is deliberately kept as a zip in git (LFS-friendly) so that the case
can be generated without any external workspace.
"""
import argparse
import os
import shutil
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ZIP = os.path.join(HERE, "data", "fluid8_mesh.zip")
REQUIRED = ("points", "faces", "owner", "neighbour", "boundary")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-root", required=True, help="case directory to create/populate")
    ap.add_argument("--mesh-zip", default=DEFAULT_ZIP, help="default: data/fluid8_mesh.zip")
    ap.add_argument("--force", action="store_true", help="overwrite an existing polyMesh")
    args = ap.parse_args()

    dst = os.path.join(args.case_root, "constant", "polyMesh")
    if os.path.isdir(dst) and not args.force:
        have = set(os.listdir(dst))
        missing = [f for f in REQUIRED if f not in have]
        if not missing:
            print("[01] polyMesh already present, skipping: %s" % dst)
            return 0
        print("[01] incomplete polyMesh (missing %s), rebuilding" % ",".join(missing))

    if not os.path.isfile(args.mesh_zip):
        sys.exit("[01] mesh archive not found: %s" % args.mesh_zip)

    tmp = os.path.join(args.case_root, ".stage01_tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    with zipfile.ZipFile(args.mesh_zip) as z:
        z.extractall(tmp)

    src = None
    for root, dirs, files in os.walk(tmp):
        if os.path.basename(root) == "polyMesh" and "points" in files:
            src = root
            break
    if src is None:
        sys.exit("[01] archive does not contain a polyMesh directory")

    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copytree(src, dst)
    shutil.rmtree(tmp, ignore_errors=True)

    have = sorted(os.listdir(dst))
    print("[01] polyMesh written: %s" % dst)
    print("[01]   files: %s" % " ".join(have))
    for f in REQUIRED:
        if f not in have:
            sys.exit("[01] missing required file: %s" % f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
