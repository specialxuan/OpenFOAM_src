#!/usr/bin/env python3
"""Stage 02 - unpack the modal set and gate it against the mesh.

Input : <repo>/case/AGARD/FDM/data/modes_4m.zip  (FluidNodeCoor + Disp1..4)
        <case-root>/constant/polyMesh/points      (from stage 01)
Output: <case-root>/mode/{FluidNodeCoor,Disp1..4,FluidPara}.csv

FluidPara.csv is generated: two header integers (19, 11 - unchanged from the
original export) followed by one initial modal velocity per mode.

Gate: the number of nodes in FluidNodeCoor.csv must equal the number of points
in the mesh. A mismatch means the modal set and the mesh come from different
meshes and must never be mixed (the library aborts on it anyway).
"""
import argparse
import os
import shutil
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ZIP = os.path.join(HERE, "data", "modes_4m.zip")
PARA_HEAD = ("19,", "11,")


def count_csv_rows(csv_path):
    """Read the modal CSV header: '<x>, <nNodes>, <nModes>,' (one node per data row)."""
    with open(csv_path) as f:
        head = f.readline()
    parts = [p.strip() for p in head.split(",") if p.strip()]
    if len(parts) < 3:
        return -1, -1
    return int(float(parts[1])), int(float(parts[2]))


def count_points(polymesh_points):
    """Number of points in an ascii/binary OpenFOAM points file."""
    with open(polymesh_points) as f:
        in_header = True
        for line in f:
            s = line.strip()
            if not s or s.startswith("//") or s.startswith("/*"):
                continue
            if in_header:
                if s == "}":
                    in_header = False
                continue
            # first token after the FoamFile header is the point count
            if s[0].isdigit():
                return int(s)
    return -1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-root", required=True)
    ap.add_argument("--modes-zip", default=DEFAULT_ZIP)
    ap.add_argument("--n-modes", type=int, default=4,
                    help="number of modal files kept in mode/ (default 4)")
    ap.add_argument("--v0", type=float, default=1e-4,
                    help="initial modal velocity written to FluidPara.csv")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    dst = os.path.join(args.case_root, "mode")
    if os.path.isdir(dst) and not args.force:
        have = set(os.listdir(dst))
        need = {"FluidNodeCoor.csv", "FluidPara.csv"} | \
               {"FluidNodeDisp%d.csv" % i for i in range(1, args.n_modes + 1)}
        if need <= have:
            print("[02] mode/ already present, skipping: %s" % dst)
            return 0

    if not os.path.isfile(args.modes_zip):
        sys.exit("[02] mode archive not found: %s" % args.modes_zip)

    tmp = os.path.join(args.case_root, ".stage02_tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    with zipfile.ZipFile(args.modes_zip) as z:
        z.extractall(tmp)

    src = None
    for root, dirs, files in os.walk(tmp):
        if "FluidNodeCoor.csv" in files:
            src = root
            break
    if src is None:
        sys.exit("[02] archive does not contain FluidNodeCoor.csv")

    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst)
    shutil.copy2(os.path.join(src, "FluidNodeCoor.csv"), dst)
    for i in range(1, args.n_modes + 1):
        f = os.path.join(src, "FluidNodeDisp%d.csv" % i)
        if not os.path.isfile(f):
            sys.exit("[02] archive lacks FluidNodeDisp%d.csv (asked for %d modes)"
                     % (i, args.n_modes))
        shutil.copy2(f, dst)
    with open(os.path.join(dst, "FluidPara.csv"), "w") as f:
        for h in PARA_HEAD:
            f.write(h + "\n")
        for _ in range(args.n_modes):
            f.write("%g,\n" % args.v0)
    shutil.rmtree(tmp, ignore_errors=True)

    # ---- gate: mode nodes vs mesh points -----------------------------------
    n_modes_nodes, n_available = count_csv_rows(os.path.join(dst, "FluidNodeCoor.csv"))
    mesh_points = os.path.join(args.case_root, "constant", "polyMesh", "points")
    if not os.path.isfile(mesh_points):
        sys.exit("[02] mesh not found - run stage 01 first (%s)" % mesh_points)
    n_mesh = count_points(mesh_points)
    print("[02] mode/ written: %s  (%d files)" % (dst, len(os.listdir(dst))))
    print("[02]   modal nodes = %d (header), available modes = %d, mesh points = %d"
          % (n_modes_nodes, n_available, n_mesh))
    if n_modes_nodes != n_mesh:
        sys.exit("[02] GATE FAILED: modal set and mesh do not match")
    if n_available < args.n_modes:
        sys.exit("[02] GATE FAILED: archive holds %d modes, %d requested"
                 % (n_available, args.n_modes))
    print("[02]   gate OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
